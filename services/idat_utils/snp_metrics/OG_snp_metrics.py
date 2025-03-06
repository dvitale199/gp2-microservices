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


# def process_vcf_snps(vcf, out_path):

#     out_colnames = ['chromosome', 'position', 'snpID', 'Sample_ID', 'Ref', 'Alt','ALLELE_A','ALLELE_B', 'BAlleleFreq', 'LogRRatio', 'R', 'Theta', 'GenTrain_Score', 'GType']

#     variant_metrics_out_df = pd.DataFrame(columns=out_colnames)
#     variant_metrics_out_df.to_csv(out_path, header=True, index=False)


#     names = get_vcf_names(vcf)        
#     vcf = pd.read_csv(vcf, comment='#', chunksize=10000, delim_whitespace=True, header=None, names=names, dtype={'#CHROM':str})
#     IIDs = [x for x in names if x not in ['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT']]

#     for chunk in vcf:
#         chunk.rename(columns={'#CHROM':'CHROM'}, inplace=True)
#         chunk_melt = chunk.melt(id_vars=['CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO'], value_vars=IIDs, value_name='metrics')
#         chunk_melt[['GT','GQ','IGC','BAF','LRR','NORMX','NORMY','R','THETA','X','Y']] = chunk_melt.metrics.str.split(':', expand=True)
#         chunk_melt.loc[:,'GenTrain_Score'] = chunk_melt.INFO.str.split(';',expand=True).iloc[:,10].str.replace('GenTrain_Score=','')
#         chunk_melt.loc[:,'ALLELE_A'] = chunk_melt.INFO.str.split(';',expand=True).iloc[:,1].str.replace('ALLELE_A=','')
#         chunk_melt.loc[:,'ALLELE_B'] = chunk_melt.INFO.str.split(';',expand=True).iloc[:,2].str.replace('ALLELE_B=','')
#         chunk_melt.drop(columns=['QUAL','FILTER','INFO','GQ','IGC','NORMX','NORMY','X','Y','metrics'], inplace=True)
#         chunk_melt.rename(columns={'variable':'sampleid'}, inplace=True)
#         chunk_melt.loc[:,'CHROM'] = chunk_melt['CHROM'].astype(str).str.replace('chr','')
#         chunk_final = chunk_melt.loc[:,['CHROM','POS','ID','sampleid','REF','ALT','GT','ALLELE_A','ALLELE_B','BAF','LRR', 'R', 'THETA', 'GenTrain_Score']]
#         gtype_map = {'0/0':'AA', '0/1':'AB', '1/1':'BB', './.':'NC'}
        
#         chunk_final.loc[:,'GType'] = chunk_final['GT'].map(gtype_map)
#         chunk_final.drop(columns=['GT'], inplace=True)
#         chunk_final.columns = ['chromosome', 'position', 'snpID', 'Sample_ID', 'Ref', 'Alt','ALLELE_A','ALLELE_B', 'BAlleleFreq', 'LogRRatio', 'R', 'Theta', 'GenTrain_Score', 'GType']

#         chunk_final.to_csv(out_path, header=False, index=False, mode='a')

        
# def calculate_maf(gtype_df):
    
#     gtypes_map = {
#         'AA': 0,
#         'AB': 1,
#         'BA': 1,
#         'BB': 2,
#         'NC': np.nan
#         }

#     gtypes = gtype_df.pivot(index='snpID', columns='Sample_ID', values='GT').replace(gtypes_map)
    
#     # count only called genotypes
#     N = gtypes.shape[1]-gtypes.isna().sum(axis=1)
#     freq = pd.DataFrame({'freq': gtypes.sum(axis=1)/(2*N)})
#     freq.loc[:,'maf'] = np.where(freq < 0.5, freq, 1-freq)
#     maf_out = freq.drop(columns=['freq']).reset_index()

#     return maf_out


# def clean_snp_metrics(metrics_in, out_path):
#     '''splits snp metrics files by chromosome and individual'''
    
#     snp_metrics = pd.read_csv(metrics_in,
#                      dtype={
#                          'chromosome':str,
#                          'position':int,
#                          'snpID':str,
#                          'Sample_ID':str,
#                          'Ref':str,
#                          'Alt':str,
#                          'ALLELE_A':int,
#                          'ALLELE_B':int,
#                          'BAlleleFreq':float,
#                          'LogRRatio':float,
#                          'R':float,
#                          'Theta':float,
#                          'GenTrain_Score':float,
#                          'GType':str
#                      })
    
    
#     alt_split = snp_metrics.loc[:,'Alt'].str.split(',', expand=True)
#     snp_metrics.loc[:,'Alt1'], snp_metrics.loc[:,'Alt2'] = alt_split.loc[:,0], alt_split.loc[:,1]

#     snp_metrics.loc[(snp_metrics['GType']=='AA') & (snp_metrics['ALLELE_A']==1), 'GT'] = 'BB'
#     snp_metrics.loc[(snp_metrics['GType']=='AA') & (snp_metrics['ALLELE_A']==0), 'GT'] = 'AA'
#     snp_metrics.loc[(snp_metrics['GType']=='BB') & (snp_metrics['ALLELE_B']==1), 'GT'] = 'BB'
#     snp_metrics.loc[(snp_metrics['GType']=='BB') & (snp_metrics['ALLELE_B']==0), 'GT'] = 'AA'

#     snp_metrics.loc[(snp_metrics['GType']=='AB'), 'GT'] = 'AB'
#     snp_metrics.loc[(snp_metrics['GType']=='NC'), 'GT'] = 'NC'
#     snp_metrics.loc[:,'GT'] = snp_metrics.loc[:,'GT'].fillna('NC')

#     # drop snps where gentrain score, theta, and r isna
# #     snp_metrics = snp_metrics.loc[(~snp_metrics['GenTrain_Score'].isna()) & (~snp_metrics['Theta'].isna()) & (~snp_metrics['R'].isna())]

#     snp_metrics.loc[snp_metrics['ALLELE_A']==0, 'a1'] = snp_metrics.loc[snp_metrics['ALLELE_A']==0,'Ref']
#     snp_metrics.loc[snp_metrics['ALLELE_A']==1, 'a1'] = snp_metrics.loc[snp_metrics['ALLELE_A']==1,'Alt1']
#     snp_metrics.loc[snp_metrics['ALLELE_B']==0, 'a2'] = snp_metrics.loc[snp_metrics['ALLELE_B']==0,'Ref']
#     snp_metrics.loc[snp_metrics['ALLELE_B']==1, 'a2'] = snp_metrics.loc[snp_metrics['ALLELE_B']==1,'Alt1']
#     snp_metrics.loc[snp_metrics['ALLELE_B']==2, 'a2'] = snp_metrics.loc[snp_metrics['ALLELE_B']==2,'Alt2']
    
#     # calculate maf for full 
#     maf_df = calculate_maf(snp_metrics)
#     snp_metrics_full = snp_metrics.merge(maf_df, how='left', on='snpID')
        
    
#     # output metrics file per sample, per chrom
#     for iid in snp_metrics_full.Sample_ID.unique():
#         for chrom in sorted(snp_metrics_full.chromosome.unique()):

#             outfile_name = f'{out_path}_{iid}_chr{chrom}.parquet'
#             out_df = snp_metrics_full.loc[(snp_metrics_full.chromosome==chrom) & (snp_metrics_full.Sample_ID==iid)]
#             out_df.to_parquet(outfile_name, compression='brotli')
            

def extract_info(info, idx, pattern):
    split_info = info.split(";")
    if idx < len(split_info):
        return split_info[idx].replace(pattern, "")
    return None


def process_vcf_snps_dask(vcf_path, out_path=None):
    
    out_colnames = ['chromosome', 'position', 'snpID', 'Sample_ID', 'Ref', 'Alt','ALLELE_A','ALLELE_B', 'BAlleleFreq', 'LogRRatio', 'R', 'Theta', 'GenTrain_Score', 'GType']

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
    


def idat_snp_metrics(idat_path, bpm, bpm_csv, egt, ref_fasta, out_path, iaap, clean_up=True, bcftools_plugins_path="/data/vitaled2/bin"):
    '''
    current structure of idat storage is such that a directory of each SentrixBarcode_A with all idats for that barcode in it
    for ex.
    1112223334
        --> 1112223334_R01C01_Red.idat
        --> 1112223334_R01C01_Grn.idat
        --> 1112223334_R01C02_Red.idat
        --> 1112223334_R01C02_Grn.idat
        etc.
        
    '''
    out_tmp = f'{out_path}/tmp'
    os.makedirs(out_tmp, exist_ok=True)
    barcode = idat_path.split('/')[-1].split('_')[0]
    barcode_out_path = f'{out_path}/{barcode}'
    os.makedirs(barcode_out_path, exist_ok=True)
    
    
    idat_to_gtc_cmd = f'\
{iaap} gencall \
{bpm} \
{egt} \
{barcode_out_path} \
-f {idat_path} \
-g \
-t 8'

    # export path to plugins temporarily for biowulf. will figure this out later
    gtc2vcf_cmd = f'\
export BCFTOOLS_PLUGINS={bcftools_plugins_path}; \
bcftools +gtc2vcf \
--no-version -Ob \
--bpm {bpm} \
--csv {bpm_csv} \
--egt {egt} \
--gtcs {barcode_out_path} \
--fasta-ref {ref_fasta} | \
bcftools norm --no-version -Oz -c w -f {ref_fasta} > {barcode_out_path}/{barcode}.vcf.gz'
# use --extra to output .tsv of other info from gtc

    # sort vcf
    sort_cmd = f'\
bcftools \
sort {barcode_out_path}/{barcode}.vcf.gz \
-T {out_tmp}/ \
-Oz -o {barcode_out_path}/{barcode}_sorted.vcf.gz'

    # split indels and snps in vcf
    ext_snps_cmd = f'\
vcftools --gzvcf \
{barcode_out_path}/{barcode}_sorted.vcf.gz \
--remove-indels \
--recode \
--recode-INFO-all \
--out {barcode_out_path}/{barcode}_sorted_snps'

#     split indels and snps in vcf
# can bring this in later if needed. for now, only snps
#     keep_indels_cmd = f'\
# vcftools --gzvcf \
# {barcode_out_path}/{barcode}_sorted.vcf.gz \
# --keep-only-indels \
# --recode \
# --recode-INFO-all \
# --out {barcode_out_path}/{barcode}_sorted_indels'

    cmds = [idat_to_gtc_cmd, gtc2vcf_cmd, sort_cmd, ext_snps_cmd]
    for cmd in cmds:
        if cmd == gtc2vcf_cmd:
            subprocess.call(cmd, shell=True)
        else:
            shell_do(cmd)
            
    # get snp info from each vcf
    vcf_in = f'{barcode_out_path}/{barcode}_sorted_snps.recode.vcf'

    # path to directory that will store parquet partitioned on sampleid
    snp_metrics_out = f'{barcode_out_path}/snp_metrics_{barcode}'
    snp_metrics = process_vcf_snps_dask(vcf_in)
    clean_snp_metrics_dask(snp_metrics, snp_metrics_out)

    if clean_up:
        # clean up large intermediates
        extensions = ['.vcf', '.gtc', '.vcf.gz']

        for file in os.listdir(barcode_out_path):
            for extension in extensions:
                if file.endswith(extension):
                    os.remove(os.path.join(barcode_out_path, file))
            

    
def cnv_qc(geno_path, out_path, maf=0.01, geno=0.02, hwe=5e-6, indep_pairwise=[1000, 10, 0.01], samples_path=None):
        
    if samples_path:
        cmd1 = f'\
plink \
--bfile {geno_path} \
--keep {samples_path} \
--make-bed \
--out {out_path}_tmp1'
        
    else:
        cmd1 = f'\
plink \
--bfile {geno_path} \
--make-bed \
--out {out_path}_tmp1'

    cmd2 = f'\
plink \
--bfile {out_path}_tmp1 \
--maf {maf} \
--geno {geno} \
--hwe {hwe} \
--autosome \
--make-bed \
--out {out_path}_tmp2'

    cmd3 = f'\
plink \
--bfile {out_path}_tmp2 \
--indep-pairwise {indep_pairwise[0]} \
{indep_pairwise[1]} \
{indep_pairwise[2]} \
--autosome \
--out {out_path}_tmp2'
        
    cmd4 = f'\
plink \
--bfile {out_path}_tmp2 \
--extract {out_path}_tmp2.prune.in \
--pca \
--make-bed \
--out {out_path}'

    cmds = [cmd1, cmd2, cmd3, cmd4]

    for cmd in cmds:
        shell_do(f'{cmd}')

    # Define a list of temporary file names and extensions to remove
    tmp_files = [f'{out_path}_tmp1', f'{out_path}_tmp2']
    tmp_exts = ['bed', 'bim', 'fam', 'log','hh', 'prune.in', 'prune.out']

    # Loop over the temporary files and extensions and remove the files if they exist
    for tmp_file in tmp_files:
        for tmp_ext in tmp_exts:
            tmp_path = f'{tmp_file}.{tmp_ext}'
            if os.path.exists(tmp_path):
                os.remove(tmp_path)


def call_cnvs(snp_metrics_file, bim_path, out_path, intervals_file, min_variants=10, kb_window=100, min_gentrain=0.2):

    # Load in the data.
    metrics_df = pd.read_parquet(snp_metrics_file)
    bim = pd.read_csv(bim_path, sep='\s+', header=None, names=['chr','id','pos','bp','a1','a2'], usecols=['id'])
    sample_df = metrics_df.loc[(metrics_df.snpID.isin(bim.id)) & (metrics_df.GenTrain_Score>=min_gentrain)]

    temp_interval_df = pd.read_csv(intervals_file, engine='c')
    temp_interval_df.drop_duplicates(subset = ["NAME"], inplace=True, keep='first')
    intervals_df = temp_interval_df.copy()


    """# Now reduce just to the intervals of interest and summarize each interval."""

    # Break down L2R and BAF per gene.

    results = []

    interval_list = intervals_df['NAME'].unique()
    
    for INTERVAL in interval_list:
      interval_CHR = intervals_df.loc[intervals_df['NAME'] == INTERVAL, 'CHR'].item()
      interval_START_gene = intervals_df.loc[intervals_df['NAME'] == INTERVAL, 'START'].item()
      interval_STOP_gene = intervals_df.loc[intervals_df['NAME'] == INTERVAL, 'STOP'].item()
      interval_START = interval_START_gene - (kb_window*1000)
      interval_STOP = interval_STOP_gene + (kb_window*1000)
      temp_df = sample_df.loc[(sample_df['chromosome'] == interval_CHR) & (sample_df['position'] >= interval_START) & (sample_df['position'] <= interval_STOP)]

      if temp_df.shape[0] < min_variants:

        results.append((INTERVAL, temp_df.shape[0], NaN, NaN, NaN, interval_START, interval_START_gene, interval_STOP_gene, interval_STOP))
        
      else:
        temp_df['BAF_insertion'] = np.where( (temp_df['BAlleleFreq'].between(0.65, 0.85, inclusive='neither')) | (temp_df['BAlleleFreq'].between(0.15, 0.35, inclusive='neither')), 1, 0)
        temp_df['L2R_deletion'] = np.where( temp_df['LogRRatio'] < -0.2, 1, 0)
        temp_df['L2R_duplication'] = np.where( temp_df['LogRRatio'] > 0.2, 1, 0)
        PERCENT_BAF_INSERTION = temp_df['BAF_insertion'].mean()
        PERCENT_L2R_DELETION = temp_df['L2R_deletion'].mean()
        PERCENT_L2R_DUPLICATION = temp_df['L2R_duplication'].mean()
        results.append((INTERVAL, temp_df.shape[0], PERCENT_BAF_INSERTION, PERCENT_L2R_DELETION, PERCENT_L2R_DUPLICATION, interval_START, interval_START_gene, interval_STOP_gene, interval_STOP))

    output = pd.DataFrame(results, columns=('INTERVAL', 'NUM_VARIANTS', 'PERCENT_BAF_INSERTION', 'PERCENT_L2R_DELETION','PERCENT_L2R_DUPLICATION','START_PLUS_WINDOW','START','STOP','STOP_PLUS_WINDOW'))
    output.to_parquet(out_path, compression='brotli')
    
    
def CNV_WAS(cnv_dosage_file, pheno, covar, out_path):
    scaler = MinMaxScaler()
    dosage_df = pd.read_csv(cnv_dosage_file)
    #fix column names
    dosage_df.columns = [x.replace('-','_') for x in dosage_df.columns]
    dosage_df.columns = [x.replace('.','_') for x in dosage_df.columns]
    dosage_df.columns = [x.replace(' ','_') for x in dosage_df.columns]

    pheno_df = pd.read_csv(pheno, sep='\t')
    covar_df = pd.read_csv(covar)


    if covar_df.age_of_onset.isna().all():
        covar_df.drop(columns=['age_of_onset'], inplace=True)
    else:
        covar_df.loc[:,'age_of_onset'] = scaler.fit_transform(covar_df[['age_of_onset']])

    if covar_df.age.isna().all():
        covar_df.drop(columns=['age'], inplace=True)
    else:
        covar_df.loc[:,'age'] = scaler.fit_transform(covar_df[['age']])

    if covar_df.sex_for_qc.isna().all():
        covar_df.drop(columns=['sex_for_qc'], inplace=True)

    covar_df.drop(columns=['FID'], inplace=True)
    covar_df.rename(columns={'GP2sampleID':'sampleid'}, inplace=True)

    data_df = dosage_df.merge(covar_df, on='sampleid', how='left').merge(pheno_df, on='sampleid', how='left').set_index('sampleid')

    rm_pred = [f'PC{i}' for i in range(1,21)] + ['sex_for_qc','age_of_onset','age','pheno']

    pred_list = [x for x in data_df.columns if x not in rm_pred]
    covars_list = [x for x in data_df.columns if x not in pred_list + [f'PC{i}' for i in range(11,21)] + ['pheno']]


    genes_list = dosage_df.columns.drop('sampleid')
    data_df_final = data_df.copy()
    ctrl_stds = data_df.loc[data_df['pheno']==0, genes_list].std()
    ctrl_means = data_df.loc[data_df['pheno']==0, genes_list].mean()

    for gene in genes_list:
        lower = (ctrl_means[gene]-(ctrl_stds[gene]*2))
        upper = (ctrl_means[gene]+(ctrl_stds[gene]*2))         
        data_df_final.loc[:, gene] = data_df_final.loc[:, gene].between(lower, upper, inclusive='both')

    data_df_final.loc[:, genes_list] = data_df_final.loc[:, genes_list].astype(int)

    results = []
    fails = []

    for pred in range(len(pred_list)):
        pred_name = pred_list[pred]
        formula = "pheno ~ " + pred_name + " + " + ' + '.join(covars_list)

        fitted = sm.formula.glm(formula=formula, family=sm.families.Binomial(), data=data_df_final).fit()
        beta_coef  = fitted.params.loc[pred_name]
        beta_se  = fitted.bse.loc[pred_name]
        p_val = fitted.pvalues.loc[pred_name]

        results.append((pred_name, beta_coef, beta_se, p_val))


    output = pd.DataFrame(results, columns=('PREDICTOR', 'BETA_COEF', 'BETA_SE','P_VAL'))
    output.to_csv(out_path, sep='\t', header=True, index=False)