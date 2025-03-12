import os
import os
import subprocess
import numpy as np
import pandas as pd
import dask.dataframe as dd
import sys
from concurrent.futures import ProcessPoolExecutor
import psutil
import glob
# Supress copy warning.
pd.options.mode.chained_assignment = None

def shell_do(command, print_cmd=False, log=False, return_log=False, err=False):
    if print_cmd:
        print(f'Executing: {(" ").join(command.split())}', file=sys.stderr)

    res=subprocess.run(command.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    output = res.stdout.decode('utf-8') + res.stderr.decode('utf-8')

    if log:
        print(output)
    if return_log:
        return output
    if err:
        return res.stderr.decode('utf-8')


def get_vcf_names(vcf_path):
    with open(vcf_path, "r") as ifile:
        for line in ifile:
            if line.startswith("#CHROM"):
                vcf_names = [x.strip('\n') for x in line.split('\t')]
                break
    ifile.close()
    return vcf_names
            

def extract_info(info, idx, pattern):
    split_info = info.split(";")
    if idx < len(split_info):
        return split_info[idx].replace(pattern, "")
    return None


def process_vcf_snps_dask(vcf_path, out_path=None):
    
    out_colnames = [
        'chromosome', 
        'position', 
        'snpID', 
        'Sample_ID', 
        'Ref', 
        'Alt',
        'ALLELE_A',
        'ALLELE_B', 
        'BAlleleFreq', 
        'LogRRatio', 
        'R', 
        'Theta', 
        'GenTrain_Score', 
        'GType'
    ]

    names = get_vcf_names(vcf_path)
    vcf = dd.read_csv(vcf_path, comment='#', blocksize='100MB', delim_whitespace=True, header=None, names=names, dtype={'#CHROM':str}, assume_missing=True)
    IIDs = [x for x in names if x not in ['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT']]

    vcf = vcf.rename(columns={'#CHROM':'CHROM'})
    vcf_melt = vcf.melt(id_vars=['CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO'], value_vars=IIDs, value_name='metrics')
    vcf_melt[['GT','GQ','IGC','BAF','LRR','NORMX','NORMY','R','THETA','X','Y']] = vcf_melt['metrics'].str.split(':', expand=True, n=10)

    vcf_melt = vcf_melt.assign(
        GenTrain_Score=vcf_melt['INFO'].map(lambda info: extract_info(info, idx=10, pattern='GenTrain_Score='), meta=('INFO', 'object')),
        ALLELE_A=vcf_melt['INFO'].map(lambda info: extract_info(info, idx=1, pattern='ALLELE_A='), meta=('INFO', 'object')),
        ALLELE_B=vcf_melt['INFO'].map(lambda info: extract_info(info, idx=2, pattern='ALLELE_B='), meta=('INFO', 'object'))
    )

    vcf_melt = vcf_melt.drop(columns=['QUAL','FILTER','INFO','GQ','IGC','NORMX','NORMY','X','Y','metrics'])
    vcf_melt = vcf_melt.rename(columns={'variable':'sampleid'})
    vcf_melt['CHROM'] = vcf_melt['CHROM'].astype(str).str.replace('chr','')
    vcf_final = vcf_melt.loc[:,['CHROM','POS','ID','sampleid','REF','ALT','GT','ALLELE_A','ALLELE_B','BAF','LRR', 'R', 'THETA', 'GenTrain_Score']]

    gtype_map = {'0/0':'AA', '0/1':'AB', '1/1':'BB', './.':'NC'}
    vcf_final['GType'] = vcf_final['GT'].map(gtype_map)
    vcf_final = vcf_final.drop(columns=['GT'])

    vcf_final.columns = out_colnames
    
    if out_path:
        vcf_final.to_parquet(out_path, compression='brotli')
    
    return vcf_final
    

def update_gt(df):
    df.loc[(df['GType'] == 'AA') & (df['ALLELE_A'] == 1), 'GT'] = 'BB'
    df.loc[(df['GType'] == 'AA') & (df['ALLELE_A'] == 0), 'GT'] = 'AA'
    df.loc[(df['GType'] == 'BB') & (df['ALLELE_B'] == 1), 'GT'] = 'BB'
    df.loc[(df['GType'] == 'BB') & (df['ALLELE_B'] == 0), 'GT'] = 'AA'
    df.loc[(df['GType'] == 'AB'), 'GT'] = 'AB'
    df.loc[(df['GType'] == 'NC'), 'GT'] = 'NC'
    df['GT'] = df['GT'].fillna('NC')
    return df


def split_alt(row):
    return pd.Series(row['Alt'].split(','))


def process_partition(partition):
    alt_split = partition.apply(split_alt, axis=1)
    alt_split.columns = [f'Alt{i+1}' for i in range(alt_split.shape[1])]
    
    # Add missing columns with NaN values if necessary
    for i in range(1, 4):
        if f'Alt{i}' not in alt_split.columns:
            alt_split[f'Alt{i}'] = np.nan
    
    return pd.concat([partition, alt_split], axis=1)


def clean_snp_metrics_dask(snp_metrics, out_path):
    '''splits snp metrics files by chromosome and individual'''
    
    snp_metrics = snp_metrics.astype(dtype={
        'chromosome':str,
        'position':int,
        'snpID':str,
        'Sample_ID':str,
        'Ref':str,
        'Alt':str,
        'ALLELE_A':int,
        'ALLELE_B':int,
        'BAlleleFreq':float,
        'LogRRatio':float,
        'R':float,
        'Theta':float,
        'GenTrain_Score':float,
        'GType':str
    })
    snp_metrics = snp_metrics.persist()
    
    # With this updated line to include the new columns in the metadata:
    new_columns = snp_metrics.columns.tolist() + ['Alt1', 'Alt2', 'Alt3']
    meta = snp_metrics._meta.assign(Alt1=pd.Series(dtype=str), Alt2=pd.Series(dtype=str), Alt3=pd.Series(dtype=str))[new_columns]
    # With this updated line to include the new columns in the metadata:
    new_columns = snp_metrics.columns.tolist() + ['Alt1', 'Alt2', 'Alt3']
    meta = snp_metrics._meta.assign(Alt1=pd.Series(dtype=str), Alt2=pd.Series(dtype=str), Alt3=pd.Series(dtype=str))[new_columns]
    snp_metrics = snp_metrics.map_partitions(process_partition, meta=meta)

    snp_metrics['GT'] = snp_metrics['GType']
    snp_metrics = snp_metrics.map_partitions(update_gt)

    snp_metrics['a1'] = snp_metrics['Ref']
    snp_metrics['a1'] = snp_metrics['a1'].mask((snp_metrics['ALLELE_A']==1), snp_metrics['Alt1'])
    snp_metrics['a2'] = snp_metrics['Ref']
    snp_metrics['a2'] = snp_metrics['a2'].mask((snp_metrics['ALLELE_B']==1), snp_metrics['Alt1'])
    snp_metrics['a2'] = snp_metrics['a2'].mask((snp_metrics['ALLELE_B']==2), snp_metrics['Alt2'])

    dd.to_parquet(snp_metrics, out_path, compression='brotli', partition_on=['Sample_ID', 'chromosome'], write_index=False)
    

def setup_directories(idat_path, out_path):
    """Set up directories for processing IDAT files."""
    out_tmp = f'{out_path}/tmp'
    os.makedirs(out_tmp, exist_ok=True)
    barcode = idat_path.split('/')[-1].split('_')[0]
    barcode_out_path = f'{out_path}/{barcode}'
    os.makedirs(barcode_out_path, exist_ok=True)
    return barcode, barcode_out_path, out_tmp




def process_snp_data(barcode_out_path, barcode):
    """Process SNP information from VCF file."""
    vcf_in = f'{barcode_out_path}/{barcode}_sorted_snps.recode.vcf'
    snp_metrics_out = f'{barcode_out_path}/snp_metrics_{barcode}'
    snp_metrics = process_vcf_snps_dask(vcf_in)
    clean_snp_metrics_dask(snp_metrics, snp_metrics_out)

def cleanup_intermediates(barcode_out_path, clean_up):
    """Clean up intermediate files if requested."""
    if clean_up:
        extensions = ['.vcf', '.gtc', '.vcf.gz']
        for file in os.listdir(barcode_out_path):
            for extension in extensions:
                if file.endswith(extension):
                    os.remove(os.path.join(barcode_out_path, file))







def convert_idat_to_gtc(iaap, bpm, egt, barcode_out_path, idat_path):
    """Convert IDAT files to GTC format.
    
    Args:
        iaap: Path to IAAP CLI executable
        bpm: Path to BPM file
        egt: Path to EGT file
        barcode_out_path: Output directory for GTC files
        idat_path: Path to IDAT directory
        
    Returns:
        list or bool: List of generated GTC files if successful, False if failed
    """
    # Set environment variable for .NET Core to run without globalization support
    env = os.environ.copy()
    env["DOTNET_SYSTEM_GLOBALIZATION_INVARIANT"] = "1"
    
    # Get initial list of GTC files before conversion
    initial_gtc_files = set(glob.glob(os.path.join(barcode_out_path, "*.gtc")))
    
    idat_to_gtc_cmd = f"{iaap} gencall {bpm} {egt} {barcode_out_path} -f {idat_path} -g -t 8"
    
    # Use env parameter to pass environment variables
    result = subprocess.run(idat_to_gtc_cmd, 
                          shell=True, 
                          stdout=subprocess.PIPE, 
                          stderr=subprocess.PIPE,
                          env=env)
    
    if result.returncode != 0:
        print(f"Command failed with exit code {result.returncode}")
        print(f"Error: {result.stderr.decode('utf-8')}")
        return False
    
    # Get the list of GTC files after conversion
    all_gtc_files = glob.glob(os.path.join(barcode_out_path, "*.gtc"))
    
    # Find new GTC files by comparing with the initial set
    new_gtc_files = [f for f in all_gtc_files if f not in initial_gtc_files]
    
    print(f"Generated {len(new_gtc_files)} GTC files")
    return new_gtc_files


def gtc_to_vcf(gtc_directory, vcf_out, bpm, bpm_csv, egt, ref_fasta, out_tmp, bcftools_plugins_path, threads=8, memory="4G"):
    """Convert all GTC files in a directory to a single VCF file.
    
    Args:
        gtc_directory: Path to directory containing GTC files
        vcf_out: Path to output VCF file
        bpm: Path to BPM file
        bpm_csv: Path to CSV file
        egt: Path to EGT file
        ref_fasta: Path to reference FASTA file
        out_tmp: Path to temporary directory
        bcftools_plugins_path: Path to bcftools plugins
        threads: Number of threads to use
        memory: Memory to allocate for sorting
    """
    
    ram_tmp = "/dev/shm/vcf_tmp" if os.path.exists("/dev/shm") else out_tmp
    os.makedirs(ram_tmp, exist_ok=True)
    
    cmd = f"""export BCFTOOLS_PLUGINS={bcftools_plugins_path} && \
bcftools +gtc2vcf \
--no-version -Ou \
--bpm {bpm} \
--csv {bpm_csv} \
--egt {egt} \
--gtcs {gtc_directory} \
--fasta-ref {ref_fasta} | \
bcftools norm \
-Ou \
--no-version \
-c w \
-f {ref_fasta} \
--threads {threads} | \
bcftools sort \
-T {ram_tmp}/ \
-m {memory} | \
bcftools view \
--threads {threads} \
-Oz \
-o {vcf_out} && \
bcftools index --threads {threads} {vcf_out}"""
    
    print(f"Processing all GTC files in: {gtc_directory}")
    result = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    
    if result.returncode != 0:
        print(f"Command failed: {result.stderr.decode('utf-8')}")
        print(f"Error message: {result.stderr.decode('utf-8')}")
        return False
    
    return True


def process_idat_files(idat_path, output_directory, bpm, bpm_csv, egt, ref_fasta, iaap, bcftools_plugins_path):
    """Process a single IDAT directory to generate SNP metrics.
    
    This function handles the full pipeline:
    1. Convert IDAT files to GTC files using IAAP tool
    2. Convert all GTC files to a single VCF using bcftools
    3. Process VCF to extract SNP metrics
    4. Clean up intermediate files
    
    Args:
        idat_path: Path to a single IDAT directory
        output_directory: Directory to output results
        bpm: Path to BPM file
        bpm_csv: Path to CSV file
        egt: Path to EGT file
        ref_fasta: Path to reference FASTA file
        iaap: Path to IAAP CLI executable
        bcftools_plugins_path: Path to bcftools plugins
        
    Returns:
        bool: True if processing was successful, False otherwise
    """
    # Verify the IDAT directory exists and contains IDAT files
    if not os.path.isdir(idat_path):
        print(f"IDAT directory not found: {idat_path}")
        return False
        
    if not any(f.endswith('.idat') for f in os.listdir(idat_path)):
        print(f"No IDAT files found in {idat_path}")
        return False
    
    # Get barcode from directory name
    barcode = os.path.basename(idat_path)
    barcode_out_path = os.path.join(output_directory, barcode)
    os.makedirs(barcode_out_path, exist_ok=True)
    out_tmp = os.path.join(output_directory, f"tmp_{barcode}")
    os.makedirs(out_tmp, exist_ok=True)
    
    print(f"Processing IDAT directory: {idat_path}")
    
    # Determine optimal resource allocation
    cpu_count = os.cpu_count()
    total_memory_gb = psutil.virtual_memory().total / (1024**3)
    threads = cpu_count
    memory = f"{int(total_memory_gb // 2)}G"
    
    print(f"Using {threads} threads and {memory} memory")
    
    # Step 1: Convert IDAT files to GTC
    print(f"Converting IDAT to GTC for {barcode}...")
    if not convert_idat_to_gtc(iaap, bpm, egt, barcode_out_path, idat_path):
        print(f"Failed to convert IDAT to GTC for {barcode}")
        return False
    
    # Verify GTC files were created
    gtc_files = glob.glob(os.path.join(barcode_out_path, "*.gtc"))
    if not gtc_files:
        print(f"No GTC files found in {barcode_out_path} after conversion")
        return False
    
    print(f"Successfully created {len(gtc_files)} GTC files")
    
    # Step 2: Process all GTC files in the directory at once
    vcf_out = f"{barcode_out_path}/{barcode}_sorted.recode.vcf.gz"
    print(f"Converting GTC files to VCF for {barcode}...")
    if not gtc_to_vcf(barcode_out_path, vcf_out, bpm, bpm_csv, egt, ref_fasta, out_tmp, bcftools_plugins_path, threads, memory):
        print(f"Failed to convert GTC files to VCF for {barcode}")
        return False
    
    # Verify VCF file was created
    if not os.path.exists(vcf_out):
        print(f"Expected VCF file not found: {vcf_out}")
        return False
    
    # Step 3: Process SNP data
    # print(f"Processing SNP data for {barcode}...")
    # try:
    #     process_snp_data(barcode_out_path, barcode)
    # except Exception as e:
    #     print(f"Error processing SNP data for {barcode}: {e}")
    #     return False
    
    # Step 4: Clean up
    # print(f"Cleaning up for {barcode}...")
    # try:
    #     # Clean up intermediates
    #     cleanup_intermediates(barcode_out_path, True)
        
    #     # Clean temp directory
    #     if os.path.exists(out_tmp):
    #         import shutil
    #         shutil.rmtree(out_tmp)
    # except Exception as e:
    #     print(f"Warning: Error during cleanup: {e}")
        # Continue despite cleanup errors
    
    print(f"Processing complete for {barcode}")
    return True



    
