from fastapi import FastAPI, HTTPException
from pydantic import BaseModel
from google.cloud import storage
import tempfile
import os
from typing import List
import shutil
from app import carriers

app = FastAPI()

class CarrierRequest(BaseModel):
    geno_bucket: str  # GCS bucket containing genotype files
    geno_prefix: str  # Common prefix for all genotype files
    snplist_path: str  # GCS path to SNP list file
    key_file_path: str  # GCS path to key file
    ancestry_labels: List[str]  # List of ancestry labels to process
    output_bucket: str  # GCS bucket for output files

def download_from_gcs(bucket_name: str, blob_path: str, local_path: str):
    """Download a file from GCS to local storage."""
    storage_client = storage.Client()
    bucket = storage_client.bucket(bucket_name)
    blob = bucket.blob(blob_path)
    blob.download_to_filename(local_path)

def upload_to_gcs(bucket_name: str, local_path: str, blob_path: str):
    """Upload a file from local storage to GCS."""
    storage_client = storage.Client()
    bucket = storage_client.bucket(bucket_name)
    blob = bucket.blob(blob_path)
    blob.upload_from_filename(local_path)

@app.post("/process_carriers")
async def process_carriers(request: CarrierRequest):
    """
    Process carrier information from genotype files stored in GCS.
    Returns paths to the generated files in GCS.
    """
    try:
        # Create temporary working directory
        with tempfile.TemporaryDirectory() as temp_dir:
            # Set up directory structure
            geno_dir = os.path.join(temp_dir, "genotypes")
            output_dir = os.path.join(temp_dir, "outputs")
            os.makedirs(geno_dir, exist_ok=True)
            os.makedirs(output_dir, exist_ok=True)

            # Download SNP list and key file
            snplist_local = os.path.join(temp_dir, "snps.csv")
            key_file_local = os.path.join(temp_dir, "key.csv")
            
            bucket, blob_path = request.snplist_path.replace("gs://", "").split("/", 1)
            download_from_gcs(bucket, blob_path, snplist_local)
            
            bucket, blob_path = request.key_file_path.replace("gs://", "").split("/", 1)
            download_from_gcs(bucket, blob_path, key_file_local)

            # Process each ancestry label
            results_by_label = {}
            for label in request.ancestry_labels:
                # Create ancestry-specific directories
                label_dir = os.path.join(geno_dir, label)
                os.makedirs(label_dir, exist_ok=True)

                # Download genotype files (.pgen, .pvar, .psam)
                for ext in ['.pgen', '.pvar', '.psam']:
                    gcs_path = f"{request.geno_prefix}/{label}/{label}_release9_vwb{ext}"
                    local_path = os.path.join(label_dir, f"{label}_release9_vwb{ext}")
                    bucket, blob_path = gcs_path.replace("gs://", "").split("/", 1)
                    download_from_gcs(request.geno_bucket, blob_path, local_path)

                # Process carriers for this ancestry
                results = carriers.extract_carriers(
                    geno_path=os.path.join(label_dir, f"{label}_release9_vwb"),
                    snplist_path=snplist_local,
                    out_path=os.path.join(output_dir, label),
                    return_dfs=True
                )
                results_by_label[label] = results

            # Combine results
            combined_results = carriers.combine_carrier_files(
                results_by_label=results_by_label,
                key_file=key_file_local,
                output_dir=output_dir
            )

            # Upload results to GCS
            gcs_results = {}
            for key, local_path in combined_results.items():
                filename = os.path.basename(local_path)
                gcs_path = f"carriers_output/{filename}"
                upload_to_gcs(request.output_bucket, local_path, gcs_path)
                gcs_results[key] = f"gs://{request.output_bucket}/{gcs_path}"

            return {
                "status": "success",
                "output_files": gcs_results
            }

    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e)) 