import os
import subprocess
import numpy as np
from numpy.core.numeric import NaN
import pandas as pd
import glob
import shutil
from sklearn.preprocessing import MinMaxScaler
import statsmodels.api as sm
import dask.dataframe as dd
import configparser
import sys

# Supress copy warning.
pd.options.mode.chained_assignment = None

# Default configuration
DEFAULT_CONFIG = {
    'bcftools_plugins_path': "/data/vitaled2/bin",
    'dask_blocksize': "100MB",
    'num_threads': "8"
}

# Load configuration from file or environment
def load_config(config_file=None):
    config = DEFAULT_CONFIG.copy()
    
    # Try loading from environment variables first
    for key in config:
        env_var = f"SNP_METRICS_{key.upper()}"
        if env_var in os.environ:
            config[key] = os.environ[env_var]
    
    # Then try config file if provided
    if config_file and os.path.exists(config_file):
        parser = configparser.ConfigParser()
        parser.read(config_file)
        if 'SNP_METRICS' in parser:
            for key in config:
                if key in parser['SNP_METRICS']:
                    config[key] = parser['SNP_METRICS'][key]
    
    return config

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


def process_vcf_snps_dask(vcf_path, out_path=None, blocksize='100MB'):
    """Process VCF SNPs using Dask with configurable blocksize"""
    out_colnames = ['chromosome', 'position', 'snpID', 'Sample_ID', 'Ref', 'Alt','ALLELE_A','ALLELE_B', 'BAlleleFreq', 'LogRRatio', 'R', 'Theta', 'GenTrain_Score', 'GType']

    names = get_vcf_names(vcf_path)
    vcf = dd.read_csv(vcf_path, comment='#', blocksize=blocksize, delim_whitespace=True, header=None, names=names, dtype={'#CHROM':str}, assume_missing=True)
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


def clean_snp_metrics_dask(snp_metrics, out_path, partition_size='100MB'):
    """Clean SNP metrics with configurable partition size"""
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
    snp_metrics = snp_metrics.repartition(partition_size=partition_size).persist()
    
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
    

def convert_idat_to_gtc(iaap, bpm, egt, barcode_out_path, idat_path):
    """Convert IDAT files to GTC format"""
    cmd = f'{iaap} gencall {bpm} {egt} {barcode_out_path} -f {idat_path} -g -t 8'
    return shell_do(cmd)

def convert_gtc_to_vcf(bcftools_plugins_path, bpm, bpm_csv, egt, barcode_out_path, ref_fasta, barcode):
    """Convert GTC files to VCF format"""
    cmd = f'export BCFTOOLS_PLUGINS={bcftools_plugins_path}; \
    bcftools +gtc2vcf \
    --no-version -Ob \
    --bpm {bpm} \
    --csv {bpm_csv} \
    --egt {egt} \
    --gtcs {barcode_out_path} \
    --fasta-ref {ref_fasta} | \
    bcftools norm --no-version -Oz -c w -f {ref_fasta} > {barcode_out_path}/{barcode}.vcf.gz'
    return subprocess.call(cmd, shell=True)

def idat_snp_metrics(idat_path, bpm, bpm_csv, egt, ref_fasta, out_path, iaap, clean_up=True, config_file=None):
    """
    Process IDAT files to extract SNP metrics.
    
    Args:
        idat_path (str): Path to IDAT files
        bpm (str): Path to BPM file
        bpm_csv (str): Path to BPM CSV file
        egt (str): Path to EGT file
        ref_fasta (str): Path to reference FASTA
        out_path (str): Output directory
        iaap (str): Path to IAAP executable
        clean_up (bool): Whether to clean up intermediate files
        config_file (str, optional): Path to configuration file
    """
    config = load_config(config_file)
    bcftools_plugins_path = config['bcftools_plugins_path']
    num_threads = config['num_threads']
    
    try:
        out_tmp = f'{out_path}/tmp'
        os.makedirs(out_tmp, exist_ok=True)
        barcode = idat_path.split('/')[-1].split('_')[0]
        barcode_out_path = f'{out_path}/{barcode}'
        os.makedirs(barcode_out_path, exist_ok=True)
        
        # Step 1: Convert IDAT to GTC
        convert_idat_to_gtc(iaap, bpm, egt, barcode_out_path, idat_path)
        
        # Step 2: Convert GTC to VCF
        convert_gtc_to_vcf(bcftools_plugins_path, bpm, bpm_csv, egt, barcode_out_path, ref_fasta, barcode)
        
        # Step 3: Sort VCF
        sort_cmd = f'bcftools sort {barcode_out_path}/{barcode}.vcf.gz -T {out_tmp}/ -Oz -o {barcode_out_path}/{barcode}_sorted.vcf.gz'
        shell_do(sort_cmd)
        
        # Step 4: Extract SNPs
        ext_snps_cmd = f'vcftools --gzvcf {barcode_out_path}/{barcode}_sorted.vcf.gz --remove-indels --recode --recode-INFO-all --out {barcode_out_path}/{barcode}_sorted_snps'
        shell_do(ext_snps_cmd)
        
        # Step 5: Process SNP metrics
        vcf_in = f'{barcode_out_path}/{barcode}_sorted_snps.recode.vcf'
        snp_metrics_out = f'{barcode_out_path}/snp_metrics_{barcode}'
        snp_metrics = process_vcf_snps_dask(vcf_in)
        clean_snp_metrics_dask(snp_metrics, snp_metrics_out)
        
        if clean_up:
            cleanup_intermediate_files(barcode_out_path)
            
    except Exception as e:
        print(f"Error in idat_snp_metrics: {e}")
        raise

def cleanup_intermediate_files(directory):
    """Clean up intermediate files from the processing pipeline"""
    extensions = ['.vcf', '.gtc', '.vcf.gz']
    for file in os.listdir(directory):
        for extension in extensions:
            if file.endswith(extension):
                os.remove(os.path.join(directory, file))