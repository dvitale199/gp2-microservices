from fastapi import FastAPI, HTTPException, Depends
from pydantic import BaseModel
from google.cloud import storage
import tempfile
import os
from typing import List
import shutil
from app import carriers
from app.security import get_api_key

app = FastAPI()

class CarrierRequest(BaseModel):
    raw_geno_path: str  # Base GCS path (e.g., "gs://gp2_carriers/api_test/raw_genotypes")
    snplist_path: str  # GCS path to SNP list file
    key_file_path: str  # GCS path to key file
    output_path: str  # Full GCS path including desired prefix (e.g., "gs://bucket/path/prefix_name")

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

ANCESTRY_LABELS = ['AAC', 'AFR', 'AJ', 'AMR', 'CAH', 'CAS', 'EAS', 'EUR', 'FIN', 'MDE', 'SAS']

@app.post("/process_carriers")
async def process_carriers(
    request: CarrierRequest,
    api_key: str = Depends(get_api_key)
):
    """
    Process carrier information from genotype files stored in GCS.
    Returns paths to the generated files in GCS.
    """
    try:
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
            base_bucket = request.raw_geno_path.replace("gs://", "").split("/")[0]
            base_prefix = "/".join(request.raw_geno_path.replace("gs://", "").split("/")[1:])
            
            # Store individual label results
            label_results = {}
            
            for label in ANCESTRY_LABELS:
                # Create ancestry-specific directories
                label_dir = os.path.join(geno_dir, label)
                os.makedirs(label_dir, exist_ok=True)

                # Download genotype files (.pgen, .pvar, .psam)
                for ext in ['.pgen', '.pvar', '.psam']:
                    gcs_path = f"{base_prefix}/{label}/{label}_release9_vwb{ext}"
                    local_path = os.path.join(label_dir, f"{label}_release9_vwb{ext}")
                    download_from_gcs(base_bucket, gcs_path, local_path)

                # Process carriers for this ancestry
                results = carriers.extract_carriers(
                    geno_path=os.path.join(label_dir, f"{label}_release9_vwb"),
                    snplist_path=snplist_local,
                    out_path=os.path.join(output_dir, label),
                    return_dfs=True
                )
                results_by_label[label] = results
                
                # Upload individual label results to GCS
                label_gcs_results = {}
                output_bucket = request.output_path.replace("gs://", "").split("/")[0]
                output_prefix = "/".join(request.output_path.replace("gs://", "").split("/")[1:])
                
                for key, local_path in results.items():
                    if key.endswith('_df'):  # Skip DataFrame objects
                        continue
                    filename = os.path.basename(local_path)
                    gcs_path = f"{output_prefix}/{label}/{filename}"
                    upload_to_gcs(output_bucket, local_path, gcs_path)
                    label_gcs_results[key] = f"gs://{output_bucket}/{gcs_path}"
                
                label_results[label] = label_gcs_results

            # Combine results with specified output path
            combined_results = carriers.combine_carrier_files(
                results_by_label=results_by_label,
                key_file=key_file_local,
                out_path=request.output_path,
                temp_dir=output_dir
            )

            # Upload combined results to GCS
            gcs_results = {}
            output_bucket = request.output_path.replace("gs://", "").split("/")[0]
            output_prefix = "/".join(request.output_path.replace("gs://", "").split("/")[1:])
            
            for key, local_path in combined_results.items():
                filename = os.path.basename(local_path)
                gcs_path = f"{output_prefix}/{filename}"
                upload_to_gcs(output_bucket, local_path, gcs_path)
                gcs_results[key] = f"gs://{output_bucket}/{gcs_path}"

            return {
                "status": "success",
                "processed_labels": ANCESTRY_LABELS,
                "label_outputs": label_results,
                "combined_outputs": gcs_results
            }

    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e)) 