import os
import sys
import pandas as pd
import batch_idat

def main():

    # hard coded paths
    iaap = f'./gp2_genotools_data/ilmn_files/iaap-cli/iaap-cli'
    bpm = f'./gp2_genotools_data/ilmn_files/NeuroBooster_20042459_A2.bpm'
    egt = f'./gp2_genotools_data/ilmn_files/recluster_09272022.egt'
    # bpm = "~/gp2-microservices/services/idat_utils/data/ilmn_utils/NeuroBooster_20042459_A2.bpm"
    # egt = "~/gp2-microservices/services/idat_utils/data/ilmn_utils/recluster_09092022.egt"
    # iaap = "/usr/local/bin/iaap-cli-linux-x64-1.1.0-sha.80d7e5b3d9c1fdfc2e99b472a90652fd3848bbc7/iaap-cli/iaap-cli"

    output_dir = "~/Desktop/gp2-microservices/services/idat_utils/data/output"

     # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)

    # input args
    key_path = sys.argv[1]
    study = sys.argv[2]

    print(f"Processing IDATs for study: {study}")

    key = pd.read_csv(key_path, sep='\t', low_memory=False)
    key = key[key['study'] == study]

    # check if previous samples exist
    # TODO: path????
    test = pd.read_csv(f'./gp2_genotools_data/merged_by_cohort_r10/GP2_merge_{study}.fam', sep='\s+', header=None)
    list_included = list(test[1])
    # the samples we still need to process are the ones that were NOT included in the prevous merged file
    key = key[~key['GP2sampleID'].isin(list_included)]
    barcode_list = list(set(list(key['SentrixBarcode_A'])))

    # Create shell script to execute -- convert IDATs to ped
    raw_plink_path = ""
    idat_path = f'./gp2_idats/${{BARCODE}}'

    idat_to_ped_cmd = f'\
    {iaap} gencall \
    {bpm} \
    {egt} \
    {raw_plink_path}/ \
    -f {idat_path} \
    -p \
    -t 8'

    with open(f'./gp2_genotools_data/batch_files/convert_idats_to_ped.sh', 'w') as f:
        f.write('#!/bin/bash\n\n')
        f.write('BARCODE=$1\n\n')
        f.write(f'{idat_to_ped_cmd}\n')
        f.close()


    count = 0 # Update the count with the last job # that you ran:
    codes_per_job = 5 # 8
    # Loop through each chunk of 8 codes and create a single job
    for chunk in batch_idat.chunk_list(barcode_list, codes_per_job):
        # Filter out any None values from the last incomplete chunk
        codes = [code for code in chunk if code is not None]

        # Generate a unique job name
        count += 1
        job_name = f'idattoped{study.lower()}{count}'

        script = ""
        # Add commands for each code in the chunk
        for code in codes:
            script += f"""
            # Make analysis script executable and run for a specific code
            chmod +x ./gp2_genotools_data/batch_files/convert_idats_to_ped.sh
            ./gp2_genotools_data/batch_files/convert_idats_to_ped.sh {code}
            """

        # Create the job with the combined script
        batch_idat.create_script_job_with_buckets(
            project_id = "gp2-release-terra",
            region = 'europe-west4',
            job_name = job_name,
            bucket_name_input= 'gp2_idats',
            bucket_name_output= 'gp2_genotools_data',
            script_text = script
        )
