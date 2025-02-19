import requests
import json
from pprint import pprint

payload = {
    "raw_geno_path": "gs://gp2tier2_vwb/release9_18122024/raw_genotypes",
    "snplist_path": "gs://gp2_carriers/snps_out.csv",
    "key_file_path": "gs://gp2_carriers/nba_app_key.csv",
    "output_path": "gs://gp2_carriers/release9_carriers"
}

response = requests.post("https://carriers-api-776926281950.europe-west4.run.app/process_carriers", json=payload)
pprint(response.json())