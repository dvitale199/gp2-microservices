import os
import shutil
import glob
import subprocess
from itertools import zip_longest
import pandas as pd
from google.cloud import storage,batch_v1

# Supress copy warning.
pd.options.mode.chained_assignment = None


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


# 2 buckets mounted: 1 for inputs, 1 for outputs
def create_script_job_with_buckets(project_id: str, region: str, job_name: str, bucket_name_input: str,
                                   bucket_name_output: str, script_text: str) -> batch_v1.Job:
# def create_script_job_with_buckets(project_id: str, region: str, job_name: str, bucket_name_input: str,
#                                    bucket_name_output: str, docker_image: str) -> batch_v1.Job:
    """
    This method shows how to create a sample Batch Job that will run
    a simple command on Cloud Compute instances.

    Args:
        project_id: project ID or project number of the Cloud project you want to use.
        region: name of the region you want to use to run the job. Regions that are
            available for Batch are listed on: https://cloud.google.com/batch/docs/get-started#locations
        job_name: the name of the job that will be created.
            It needs to be unique for each project and region pair.
        bucket_name: name of the bucket to be mounted for your Job.

    Returns:
        A job object representing the job created.
    """
    client = batch_v1.BatchServiceClient()

    # Define what will be done as part of the job.
    runnable = batch_v1.Runnable()

    # This would replace runnable script
    # runnable.container = batch_v1.Runnable.Container()
    # runnable.container.image_uri = docker_image
    # runnable.container.entry_point = "/bin/sh"

    runnable.script = batch_v1.Runnable.Script()
    runnable.script.text = script_text

    task = batch_v1.TaskSpec()
    task.runnables = [runnable]

    # Define the first GCS bucket
    gcs_bucket_input = batch_v1.GCS(remote_path=bucket_name_input)
    gcs_volume_input = batch_v1.Volume(gcs=gcs_bucket_input, mount_path="./gp2_idats")

    # Define the second GCS bucket
    gcs_bucket_output = batch_v1.GCS(remote_path=bucket_name_output)
    gcs_volume_output = batch_v1.Volume(gcs=gcs_bucket_output, mount_path="./gp2_genotools_data")

    # Add both volumes to the task
    task.volumes = [gcs_volume_input, gcs_volume_output]

    # We can specify what resources are requested by each task.
    resources = batch_v1.ComputeResource()
    resources.cpu_milli = 8000  # in milliseconds per cpu-second. This means the task requires 50% of a single CPUs; 2000 = 2 CPUs
    resources.memory_mib = 32768 # 4 GiB of memory (4096 MiB - mebibyte); 32768 MiB = 32 GB; 16384 = 16 GB
    task.compute_resource = resources

    task.max_retry_count = 0 # I'd suggest changing this to 0; otherwise, it tries two times before it fails
    task.max_run_duration = "3600s" # 24 hours = 86400s; 1 hour = 3600s; 8 hours = 28800s

    # Tasks are grouped inside a job using TaskGroups.
    # Currently, it's possible to have only one task group.
    group = batch_v1.TaskGroup()
    group.task_count = 1
    group.task_spec = task

    # Policies are used to define on what kind of virtual machines the tasks will run on.
    # Read more about machine types here: https://cloud.google.com/compute/docs/machine-types
    allocation_policy = batch_v1.AllocationPolicy()
    policy = batch_v1.AllocationPolicy.InstancePolicy()
    policy.machine_type = "e2-standard-8"
    instances = batch_v1.AllocationPolicy.InstancePolicyOrTemplate()
    instances.policy = policy
    allocation_policy.instances = [instances]

    job = batch_v1.Job()
    job.task_groups = [group]
    job.allocation_policy = allocation_policy
    # job.labels = {"env": "testing", "type": "script", "mount": "bucket"}
    job.labels = {"env": "testing", "type": "container", "mount": "bucket"}
    # We use Cloud Logging as it's an out of the box available option
    job.logs_policy = batch_v1.LogsPolicy()
    job.logs_policy.destination = batch_v1.LogsPolicy.Destination.CLOUD_LOGGING

    create_request = batch_v1.CreateJobRequest()
    create_request.job = job
    create_request.job_id = job_name
    # The job's parent is the region in which the job will run
    create_request.parent = f"projects/{project_id}/locations/{region}"

    return client.create_job(create_request)


def chunk_list(iterable, n):
    """Helper function to chunk the list into groups of 8."""
    args = [iter(iterable)] * n
    return zip_longest(*args)


def convert_idat_to_ped(key, study, iaap, bpm, egt, raw_plink_path, idat_path, missing_idat_dir, map_file):
    """Convert IDAT files to PED format."""

    # Set environment variable for .NET Core to run without globalization support
    env = os.environ.copy()
    env["DOTNET_SYSTEM_GLOBALIZATION_INVARIANT"] = "1"

    # Get initial list of PED files before conversion
    initial_ped_files = set(glob.glob(os.path.join(raw_plink_path, "*.ped")))

    idat_to_ped_cmd = f'\
    {iaap} gencall \
    {bpm} \
    {egt} \
    {raw_plink_path}/ \
    -f {idat_path} \
    -p \
    -t 8'

    # create script to convert idat to ped
    with open(f'./gp2_genotools_data/batch_files/convert_idats_to_ped.sh', 'w') as f:
        f.write('#!/bin/bash\n')
        f.write('BARCODE=$1\n\n')
        f.write(f'{idat_to_ped_cmd}\n')
    f.close()

    # convert idats to ped
    count = 0 # Update the count with the last job # that you ran
    codes_per_job = 5 # 8
    # Loop through each chunk of 8 codes and create a single job
    barcode_list = list(set(list(key['SentrixBarcode_A'])))
    for chunk in chunk_list(barcode_list, codes_per_job):
        # Filter out any None values from the last incomplete chunk
        codes = [code for code in chunk if code is not None]

        # Generate a unique job name
        count += 1
        job_name = f'idattoped{study.lower()}{count}'

        script = """
            #!/bin/bash

            # Install and activate PLINK
            chmod +x ./gp2_genotools_data/batch_files/plink2_module.sh
            ./gp2_genotools_data/batch_files/plink2_module.sh
        """

        # Add commands for each code in the chunk
        for code in codes:
            script += f"""
            # Make analysis script executable and run for a specific code
            chmod +x ./gp2_genotools_data/batch_files/convert_idats_to_ped.sh
            ./gp2_genotools_data/batch_files/convert_idats_to_ped.sh {code}
            """

        # Create the job with the combined script
        create_script_job_with_buckets(
            project_id = "gp2-release-terra",
            region = 'europe-west4',
            job_name = job_name,
            bucket_name_input= 'gp2_idats',
            bucket_name_output= 'gp2_genotools_data',
            script_text = script
        )

    # Get the list of PED files after conversion
    all_ped_files = glob.glob(os.path.join(raw_plink_path, "*.ped"))

    # Find new PED files by comparing with the initial set
    new_ped_files = [f for f in all_ped_files if f not in initial_ped_files]

    # copy map file to match name of each ped
    # check for missing ped files
    missing_peds = []
    missing_cnt = 0 # check for missing ped files
    for filename in key.IID:
        ped = f'{raw_plink_path}/{filename}.ped'
        out_map = f'{raw_plink_path}/{filename}.map'
        if os.path.isfile(ped):
            shutil.copyfile(src=map_file, dst=out_map)
        else:
            missing_cnt += 1
            missing_peds.append(filename)

    # create missing_ped file for reference
    with open(f'{missing_idat_dir}/missing_peds_{study}.txt', 'w') as f:
        for m_ped in missing_peds:
            f.write(f'{m_ped}\n')
    f.close()

    # return list of generated ped files
    return new_ped_files


def convert_ped_to_bed(key, study, raw_plink_path, missing_idat_dir):
    """Convert PED files to BED format."""

    # create script to convert ped to bed
    ped_file = f'./gp2_genotools_data/ped_bed/${{FILENAME}}'
    make_bed_cmd = f"./exec/plink_module.sh --file {ped_file} --make-bed --out {ped_file}"
    with open(f'./gp2_genotools_data/batch_files/convert_ped_to_bed.sh', 'w') as f:
        f.write('#!/bin/bash\n\n')
        f.write('FILENAME=$1\n\n')
        f.write(f'{make_bed_cmd}\n')
    f.close()

    filename_list = list(key['IID'])

    # Define the number of codes to process per job (8)(Up this amount x3)
    codes_per_job = 60
    # Loop through each chunk of 8 codes and create a single job
    for chunk in chunk_list(filename_list, codes_per_job):
        # Filter out any None values from the last incomplete chunk
        codes = [code for code in chunk if code is not None]

        # Generate a unique job name
        count += 1
        job_name = f'pedtobed{study.lower()}{count}'

        # Create a combined script for the N codes
        script = """
            #!/bin/bash

            # Install and activate PLINK
            chmod +x ./gp2_genotools_data/batch_files/plink_module.sh
            ./gp2_genotools_data/batch_files/plink_module.sh
        """

        # Add commands for each code in the chunk
        for code in codes:
            script += f"""
            # Make analysis script executable and run for a specific code
            chmod +x ./gp2_genotools_data/batch_files/convert_ped_to_bed.sh
            ./gp2_genotools_data/batch_files/convert_ped_to_bed.sh {code}
            """

        # Create the job with the combined script
        create_script_job_with_buckets(
            project_id="gp2-release-terra",
            region='europe-west4',
            job_name=job_name,
            bucket_name_input='gp2_genotools_data',
            bucket_name_output='gp2_genotools_data',
            script_text=script
        )

    # create file of samples to merge
    missing_beds = []
    df = key[key['study']==study]
    missing_cnt = 0
    with open(f"{raw_plink_path}/merge_bed_{study}.list", 'w') as f:
        for filename in df.IID:
            bed = f'{raw_plink_path}/{filename}'
            if os.path.isfile(f'{bed}.bed'):
                f.write(f'{bed}\n')
            else:
                missing_cnt += 1
                missing_beds.append(filename)
    f.close()

    # create missing_beds file for reference
    with open(f'{missing_idat_dir}/missing_beds_{study}.txt', 'w') as f:
        for m_bed in missing_beds:
            f.write(f'{m_bed}\n')
    f.close()

    # return path to merge_list
    return f"{raw_plink_path}/merge_bed_{study}.list"


def merge_bed_files(study, raw_plink_path, clin_key_dir):
    """merge cohort files"""

    # Create shell script to execute -- merge cohort files

    # Establish paths for inputs/outputs
    out_path = f'./gp2_genotools_data/merged_by_cohort_r10/GP2_merge_{study}'
    out_path = f'./gp2_genotools_data/tests/GP2_merge_{study}_extra'

    plink_merge_cmd = f"/home/levineks/bin/plink1.9/plink --merge-list {raw_plink_path}/merge_bed_{study}.list --update-ids {clin_key_dir}/update_ids_{study}.txt --make-bed --out {out_path}"

    with open(f'./gp2_genotools_data/batch_files/merge_by_cohort.sh', 'w') as f:
        f.write('#!/bin/bash\n\n')
        f.write(f'{plink_merge_cmd}\n')
    f.close()

    # Loop through chromosomes and create separate jobs
    # Generate a unique job name
    job_name = f'mergebycohort{study.lower()}'

    # Create a script for the specific chromosome
    script = f"""
        #!/bin/bash

        # Install and activate PLINK2
        chmod +x ./gp2_genotools_data/batch_files/plink_module.sh
        ./gp2_genotools_data/batch_files/plink_module.sh

        # Make analysis script executable and run for a specific chromosome
        chmod +x ./gp2_genotools_data/batch_files/merge_by_cohort.sh
        ./gp2_genotools_data/batch_files/merge_by_cohort.sh
    """

    # Create the job
    create_script_job_with_buckets(
        project_id="gp2-release-terra",
        region='europe-west4',
        job_name=job_name,
        bucket_name_input='gp2_genotools_data',
        bucket_name_output='gp2_genotools_data',
        script_text=script
    )