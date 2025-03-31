from fastapi import FastAPI, HTTPException, Depends
from pydantic import BaseModel
from google.cloud import storage
import tempfile
import json
import os
from typing import List
import shutil
from browser_prep import prep_browser_files
from security import get_api_key

app = FastAPI()

class BrowserRequest(BaseModel):
    # master_key_path: str # WORK IN PROGRESS - depends on restructuring of keys
    gt_out_path: str  # GCS path to GT outputs
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


@app.post("/prep_browser")
async def prep_browser(
    request: BrowserRequest,
    api_key: str = Depends(get_api_key)
):
    """
    Process Master Keys and GenoTools outputs stored in GCS.
    Returns paths to the generated files in GCS.
    """
    try:
        with tempfile.TemporaryDirectory() as temp_dir:
            # Set up directory structure
            gt_dir = os.path.join(temp_dir, "genotools")
            output_dir = os.path.join(temp_dir, "outputs")
            os.makedirs(gt_dir, exist_ok=True)
            os.makedirs(output_dir, exist_ok=True)

            # Download SNP list and key file
            gt_out_local = os.path.join(temp_dir, "gt_out.json")
            
            bucket, blob_path = request.gt_out_path.replace("gs://", "").split("/", 1)
            download_from_gcs(bucket, blob_path, gt_out_local)

            # Combine results with specified output path
            final_files = prep_browser_files(
                gt_output = gt_out_local,
                temp_dir = output_dir
            )

            # Upload combined results to GCS
            final_gcs_paths = []
            output_bucket = request.output_path.replace("gs://", "").split("/")[0]
            output_prefix = "/".join(request.output_path.replace("gs://", "").split("/")[1:])
            
            for filename in final_files:
                gcs_path = f"{output_prefix}/{filename}"
                upload_to_gcs(output_bucket, f"{output_dir}/{filename}", gcs_path)
                final_gcs_paths.append(f"gs://{output_bucket}/{gcs_path}")

            return {
                "status": "success",
                "final_ouput_paths": final_gcs_paths
            }

    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e)) 