import requests
import json
from pprint import pprint

payload = {
    "raw_geno_path": "gs://gp2tier2_vwb/release9_18122024/raw_genotypes",
    "snplist_path": "gs://gp2_carriers/api_test/snps_out.csv",
    "key_file_path": "gs://gp2_carriers/api_test/nba_app_key.csv",
    "output_path": "gs://gp2_carriers/api_test/output/release0_carriers"
}

response = requests.post("http://localhost:8000/process_carriers", json=payload)
pprint(response.json())