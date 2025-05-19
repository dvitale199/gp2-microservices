import requests
import json
from pprint import pprint
from google.cloud import secretmanager

def get_secret(project_id: str, secret_id: str, version_id: str = "latest") -> str:
    """Retrieve a secret from Secret Manager."""
    client = secretmanager.SecretManagerServiceClient()
    name = f"projects/{project_id}/secrets/{secret_id}/versions/{version_id}"
    response = client.access_secret_version(request={"name": name})
    return response.payload.data.decode("UTF-8")

# Get API key from Secret Manager
api_key = get_secret("gp2-release-terra", "browser-api-key")

### only use this for debugging
# print("Raw API key:", api_key)

payload = {
    "release_num": 10,
    "master_key_path": "gs://gp2_release10/genotools_output/master_key_release10_final_vwb.csv",
    "gt_out_path": "gs://gp2_release10/genotools_output/GP2_r10_final_post_genotools.json",
    "output_path": "gs://gt_app_utils/testing"
}

headers = {
    "X-API-Key": api_key
}

# for local testing - uvicorn main:app --host 127.0.0.1 --port 8000 --reload
url = 'http://127.0.0.1:8000/prep_browser' 
response = requests.post(
    url,
    json=payload,
    headers=headers
)
pprint(response.json())