import os
import os
import subprocess
import numpy as np
import pandas as pd
import dask.dataframe as dd
import sys
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

def convert_idat_to_gtc(iaap, bpm, egt, barcode_out_path, idat_path):
    """Convert IDAT files to GTC format."""
    # Set environment variable for .NET Core to run without globalization support
    env = os.environ.copy()
    env["DOTNET_SYSTEM_GLOBALIZATION_INVARIANT"] = "1"
    
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

def convert_gtc_to_vcf(barcode_out_path, barcode, bpm, bpm_csv, egt, ref_fasta, bcftools_plugins_path):
    """Convert GTC files to VCF format."""
    gtc2vcf_cmd = f"""export BCFTOOLS_PLUGINS={bcftools_plugins_path}; \
bcftools +gtc2vcf \
--no-version -Ob \
--bpm {bpm} \
--csv {bpm_csv} \
--egt {egt} \
--gtcs {barcode_out_path} \
--fasta-ref {ref_fasta} | \
bcftools norm --no-version -Oz -c w -f {ref_fasta} > {barcode_out_path}/{barcode}.vcf.gz"""
    subprocess.run(gtc2vcf_cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

def sort_vcf(barcode_out_path, barcode, out_tmp):
    """Sort VCF files."""
    sort_cmd = f"""bcftools \
sort {barcode_out_path}/{barcode}.vcf.gz \
-T {out_tmp}/ \
-Oz -o {barcode_out_path}/{barcode}_sorted.vcf.gz"""
    subprocess.run(sort_cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

def extract_snps(barcode_out_path, barcode):
    """Extract SNPs from sorted VCF."""
    ext_snps_cmd = f"""vcftools --gzvcf \
{barcode_out_path}/{barcode}_sorted.vcf.gz \
--recode \
--recode-INFO-all \
--out {barcode_out_path}/{barcode}_sorted_snps"""
    subprocess.run(ext_snps_cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

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

def idat_snp_metrics(idat_path, bpm, bpm_csv, egt, ref_fasta, out_path, iaap, clean_up=True, bcftools_plugins_path="/data/vitaled2/bin"):
    '''
    Process IDAT files to extract SNP metrics.
    
    current structure of idat storage is such that a directory of each SentrixBarcode_A with all idats for that barcode in it
    for ex.
    1112223334
        --> 1112223334_R01C01_Red.idat
        --> 1112223334_R01C01_Grn.idat
        --> 1112223334_R01C02_Red.idat
        --> 1112223334_R01C02_Grn.idat
        etc.
    '''
    # Setup directories
    barcode, barcode_out_path, out_tmp = setup_directories(idat_path, out_path)
    
    # Convert IDAT to GTC
    # convert_idat_to_gtc(iaap, bpm, egt, barcode_out_path, idat_path)
    
    # Convert GTC to VCF
    convert_gtc_to_vcf(barcode_out_path, barcode, bpm, bpm_csv, egt, ref_fasta, bcftools_plugins_path)
    
    # # Sort VCF
    # sort_vcf(barcode_out_path, barcode, out_tmp)
    
    # # Extract SNPs
    # extract_snps(barcode_out_path, barcode)
    
    # # Process SNP data
    # process_snp_data(barcode_out_path, barcode)
    
    # # Clean up intermediate files
    # cleanup_intermediates(barcode_out_path, clean_up)

    
