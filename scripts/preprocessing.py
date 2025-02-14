import pandas as pd
import glob

plink2 = '/home/vitaled2/.genotools/misc/executables/plink2'
project_dir = '/home/vitaled2/gp2_carriers'
data_dir = f'{project_dir}/data'
report_path = f'{data_dir}/NBA_Report.csv'
geno_dir = f'{data_dir}/raw_genotypes'
labels = ['AAC', 'AFR', 'AJ', 'AMR', 'CAH', 'CAS', 'EAS', 'EUR', 'FIN','MDE', 'SAS']
key = pd.read_csv(f'{data_dir}/nba_app_key.csv')

report = pd.read_csv(report_path)
report.columns = [x.strip() for x in report.columns]
report.drop(columns=['GP2ID','Cohort Code', 'Local ID',	'PI','PD EUR PRS (GWAS vXX)','NBA derived Ancestry'], inplace=True)

snp_df = report.loc[~report['Unnamed: 0'].isna()]
snp_df.set_index('Unnamed: 0', inplace=True)
snp_df = snp_df.transpose()
snp_df.reset_index(inplace=True)
snp_df.rename(columns={'index':'snp'}, inplace=True)
snp_df.drop(columns=['Alternative variant name'], inplace=True)
snp_df.columns.name = None
snp_df.columns = ['snp','rsid','hg19','hg38','ancestry']
snp_df = snp_df.loc[~snp_df.hg38.isna()]

pvar_total = pd.DataFrame()
for label in labels:
    label_geno_path = f'{geno_dir}/{label}/{label}_release9_vwb'
    pvar = pd.read_csv(f'{label_geno_path}.pvar', sep='\t', dtype={'#CHROM':str})
    pvar.columns = ['chrom', 'pos', 'id', 'a1', 'a2']
    pvar.loc[:,'snpid1'] = pvar['chrom'] + ':' + pvar['pos'].astype(str) + ':' + pvar['a1'] + ':' + pvar['a2']
    pvar.loc[:,'snpid2'] = pvar['chrom'] + ':' + pvar['pos'].astype(str) + ':' + pvar['a2'] + ':' + pvar['a1']
    pvar_sub1 = pvar.loc[pvar['snpid1'].isin(snp_df['hg38'])].copy()
    pvar_sub1.loc[:,'snpid'] = pvar_sub1['snpid1']
    pvar_sub2 = pvar.loc[pvar['snpid2'].isin(snp_df['hg38'])].copy()
    pvar_sub2.loc[:,'snpid'] = pvar_sub2['snpid2']
    pvar_total = pd.concat([pvar_total, pvar_sub1, pvar_sub2], ignore_index=True)
pvar_total.drop(columns=['snpid1','snpid2'], inplace=True)

snp_df_out = pvar_total.merge(snp_df, how='left', left_on='snpid', right_on='hg38')
snp_df_out.drop(columns=['chrom', 'pos', 'a1', 'a2'], inplace=True)
snp_df_out.rename(columns={'snp':'variant'}, inplace=True)

snp_df_out['snp_name'] = snp_df_out['variant'].str.extract(f'\((.*?)\)')

snp_df_out.loc[:,'locus'] = snp_df_out.variant.str.split(' ', expand=True)[:][0]
snp_df_out.loc[:,'snp_name_full'] = snp_df_out['locus'] + '_' + snp_df_out['snp_name'] + '_' + snp_df_out['id']
snp_df_out.drop(columns=['snpid'], inplace=True)
snp_df_out.drop_duplicates(inplace=True)
snp_df_out.to_csv(f'{data_dir}/snps_out.csv', index=False)





######### TESTING #########
# from services.carriers_api.app.carriers import extract_carriers, combine_carrier_files

# output_dir = f'{data_dir}/api_test/outputs'

# results_by_label = {}
# for label in labels:
#     results = extract_carriers(
#         geno_path=f'{geno_dir}/{label}/{label}_release9_vwb',
#         snplist_path=f'{data_dir}/snps_out.csv',
#         out_path=f'{output_dir}/{label}',
#         return_dfs=True
#     )
#     results_by_label[label] = results

# combined_results = combine_carrier_files(
#     results_by_label=results_by_label,
#     key_file=f'{data_dir}/nba_app_key.csv',
#     output_dir=output_dir
# )