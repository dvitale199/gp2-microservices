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
api_key = get_secret("gp2-release-terra", "carriers-api-key")

## only use this for debugging
# print("Raw API key:", api_key)
# print("Character analysis:")
# for char in api_key:
#     print(f"'{char}': ASCII value {ord(char)}")

payload = {
    "raw_geno_path": "gs://gp2tier2_vwb/release9_18122024/raw_genotypes",
    "snplist_path": "gs://gp2_carriers/snps_out.csv",
    "key_file_path": "gs://gp2_carriers/nba_app_key.csv",
    "output_path": "gs://gp2_carriers/release9_carriers"
}

headers = {
    "X-API-Key": api_key
}
url = "https://carriers-api-776926281950.europe-west4.run.app/process_carriers"
# for local testing
# url = 'http://127.0.0.1:8000/process_carriers' 
response = requests.post(
    url,
    json=payload,
    headers=headers
)
pprint(response.json())