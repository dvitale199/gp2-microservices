import pandas as pd
import subprocess
import os

def genotype_to_string(g, snp_allele):
    if pd.isna(g):
        return g
    g = int(g)
    if g == 2:
        return "WT/WT"
    elif g == 1:
        return f"WT/{snp_allele}"
    elif g == 0:
        return f"{snp_allele}/{snp_allele}"
    else:
        return ":/:"
    

def extract_carriers(geno_path: str, snplist_path: str, out_path: str, return_dfs: bool = False) -> dict:
    """
    Extract carrier information for given SNPs from a PLINK2 dataset.
    
    Args:
        geno_path: Path to PLINK2 files prefix (without .pgen/.pvar/.psam extension)
        snplist_path: Path to file containing SNP information
        out_path: Output path prefix for generated files
        return_dfs: If True, return DataFrames along with file paths
    
    Returns:
        dict: Paths to generated carrier files and optionally the DataFrames themselves
    """
    snp_df = pd.read_csv(snplist_path)
    
    temp_snps_path = f"{geno_path}_temp_snps.txt"
    snp_df['id'].to_csv(temp_snps_path, header=False, index=False)
    
    plink_out = f"{geno_path}_snps"
    extract_cmd = f"plink2 --pfile {geno_path} --extract {temp_snps_path} --export Av --freq --out {plink_out}"
    subprocess.run(extract_cmd, shell=True, check=True)
    
    traw = pd.read_csv(f"{plink_out}.traw", sep='\t')
    traw_merged = snp_df.merge(traw, how='left', left_on='id', right_on='SNP')
    
    # eventually accept file with only id, chr, pos, a1, a2 and merge results later
    colnames = [
        'id', 'rsid', 'hg19', 'hg38', 'ancestry',
        'CHR', 'SNP', '(C)M', 'POS', 'COUNTED', 'ALT',
        'variant', 'snp_name', 'locus', 'snp_name_full'
    ]
    var_cols = [x for x in colnames if x not in ['snp_name_full']]
    sample_cols = list(traw_merged.drop(columns=colnames).columns)
    
    # Process final traw data
    traw_final = traw_merged.loc[:, colnames + sample_cols]
    
    # Create string format output
    traw_out = traw_final.copy()
    traw_out[sample_cols] = traw_out.apply(
        lambda row: [genotype_to_string(row[col], row['snp_name']) for col in sample_cols],
        axis=1,
        result_type='expand'
    )
    
    # Process and save frequency info
    freq = pd.read_csv(f"{plink_out}.afreq", sep='\t')
    freq.rename(columns={'ID':'SNP'}, inplace=True)
    var_info_df = traw_final.loc[:,var_cols]
    var_info_df = var_info_df.merge(freq, how='left', on='SNP')
    var_info_df.to_csv(f"{out_path}_var_info.csv", index=False)
    
    # Process and save string format
    carriers_string = traw_out.drop(columns=var_cols).set_index('snp_name_full').T.reset_index()
    carriers_string.columns.name = None
    carriers_string = carriers_string.fillna('./.')
    carriers_string = carriers_string.astype(str)
    carriers_string.rename(columns={'index':'IID'}, inplace=True)
    carriers_string.to_csv(f"{out_path}_carriers_string.csv", index=False)
    
    # Process and save integer format
    carriers_int = traw_final.drop(columns=var_cols).set_index('snp_name_full').T.reset_index()
    carriers_int.columns.name = None
    carriers_int.rename(columns={'index':'IID'}, inplace=True)
    carriers_int.to_csv(f"{out_path}_carriers_int.csv", index=False)
    
    os.remove(temp_snps_path)
    
    result = {
        'var_info': f"{out_path}_var_info.csv",
        'carriers_string': f"{out_path}_carriers_string.csv",
        'carriers_int': f"{out_path}_carriers_int.csv"
    }
    
    if return_dfs:
        result.update({
            'var_info_df': var_info_df,
            'carriers_string_df': carriers_string,
            'carriers_int_df': carriers_int
        })
    
    return result


def combine_carrier_files(results_by_label: dict, key_file: str, out_path: str, temp_dir: str) -> dict:
    """
    Combine carrier files from multiple ancestry labels into consolidated output files.
    
    Args:
        results_by_label: Dictionary mapping ancestry labels to their extract_carriers results
        key_file: Path to key file containing study information
        out_path: Full GCS path including desired prefix
        temp_dir: Temporary directory to store files before GCS upload
    
    Returns:
        dict: Paths to local temporary files
    """
    carriers_string_full = pd.DataFrame()
    carriers_int_full = pd.DataFrame()
    
    # Get base variant info from first label and drop frequency columns
    var_info_base = next(iter(results_by_label.values()))['var_info_df' if 'var_info_df' in next(iter(results_by_label.values())) else 'var_info']
    if isinstance(var_info_base, str):
        var_info_base = pd.read_csv(var_info_base)
    freq_cols = ['ALT_FREQS', 'OBS_CT']
    var_info_base = var_info_base.drop(columns=freq_cols, errors='ignore')
    
    # Read key file
    key = pd.read_csv(key_file)
    
    # Process each ancestry label's results
    for label, results in results_by_label.items():
        # Get DataFrames (either directly or from files)
        if 'var_info_df' in results:
            label_var_info = results['var_info_df']
            carriers_string = results['carriers_string_df']
            carriers_int = results['carriers_int_df']
        else:
            label_var_info = pd.read_csv(results['var_info'])
            carriers_string = pd.read_csv(results['carriers_string'])
            carriers_int = pd.read_csv(results['carriers_int'])
        
        # Add frequency data for this population
        var_info_base[f'ALT_FREQS_{label}'] = label_var_info['ALT_FREQS']
        var_info_base[f'OBS_CT_{label}'] = label_var_info['OBS_CT']
        
        # Process string format carriers
        carriers_string['IID'] = carriers_string['IID'].str.replace('0_', '')
        carriers_string.loc[:,'ancestry'] = label
        carriers_string_full = pd.concat([carriers_string_full, carriers_string], ignore_index=True)
        
        # Process integer format carriers
        carriers_int['IID'] = carriers_int['IID'].str.replace('0_', '')
        carriers_int.loc[:,'ancestry'] = label
        carriers_int_full = pd.concat([carriers_int_full, carriers_int], ignore_index=True)
    
    # Get variant columns (excluding metadata columns)
    variant_columns = [x for x in carriers_string_full.columns if x not in ['IID','ancestry']]
    
    # Process string format output
    carriers_string_full_out = carriers_string_full[['IID', 'ancestry'] + variant_columns]
    carriers_string_full_out[variant_columns] = carriers_string_full_out[variant_columns].fillna('./.')
    carriers_string_full_out_merge = carriers_string_full_out.merge(key[['IID','study']], how='left', on='IID')
    carriers_string_final = carriers_string_full_out_merge[['IID', 'study', 'ancestry'] + variant_columns]
    
    # Process integer format output
    carriers_int_full_out = carriers_int_full[['IID', 'ancestry'] + variant_columns]
    carriers_int_full_out_merge = carriers_int_full_out.merge(key[['IID','study']], how='left', on='IID')
    carriers_int_final = carriers_int_full_out_merge[['IID', 'study', 'ancestry'] + variant_columns]
    
    # Use only the base name for local files
    base_name = os.path.basename(out_path)
    
    # Create local paths without gs:// prefix
    carriers_string_path = os.path.join(temp_dir, f'{base_name}_string.csv')
    carriers_int_path = os.path.join(temp_dir, f'{base_name}_int.csv')
    var_info_path = os.path.join(temp_dir, f'{base_name}_info.csv')
    
    # Save to local temporary files
    carriers_string_final.to_csv(carriers_string_path, index=False)
    carriers_int_final.to_csv(carriers_int_path, index=False)
    var_info_base.to_csv(var_info_path, index=False)
    
    return {
        'carriers_string': carriers_string_path,
        'carriers_int': carriers_int_path,
        'var_info': var_info_path
    }


def verify_genotype(raw_value, string_value, snp_name):
    if pd.isna(raw_value):
        return string_value == './.'
    raw_value = int(raw_value)
    if raw_value == 2:
        return string_value == "WT/WT"
    elif raw_value == 1:
        return string_value.startswith("WT/")
    elif raw_value == 0:
        return string_value.count(snp_name.split('_')[1]) == 2
    else:
        return string_value == ":/:"

def validate_carrier_data(traw_dir, combined_file, snp_info_file, samples_to_check=None):
    """
    Validate that combined carrier data matches original traw files.
    
    Args:
        traw_dir: Directory containing traw files
        combined_file: Path to combined carriers file
        snp_info_file: Path to SNP info file
        samples_to_check: List of sample IDs to check (None for all)
    """
    carriers = pd.read_csv(combined_file)
    snp_df = pd.read_csv(snp_info_file)
    
    if samples_to_check is None:
        samples_to_check = carriers['IID'].unique()[:5]  # Check first 5 samples by default
    
    mismatches = []
    for sample_id in samples_to_check:
        # Get ancestry for this sample
        ancestry = carriers.loc[carriers['IID'] == sample_id, 'ancestry'].iloc[0]
        
        # Read corresponding traw file
        traw_path = f'{traw_dir}/{ancestry}/{ancestry}_release9_vwb_snps.traw'
        traw = pd.read_csv(traw_path, sep='\t')
        traw_merged = snp_df.merge(traw, how='left', left_on='id', right_on='SNP')
        
        # Check each SNP
        for snp in snp_df['snp_name_full']:
            raw_value = traw_merged.loc[
                traw_merged['snp_name_full'] == snp, 
                f'0_{sample_id}'
            ].iloc[0]
            
            combined_value = carriers.loc[
                carriers['IID'] == sample_id,
                snp
            ].iloc[0]
            
            if not verify_genotype(raw_value, combined_value, snp):
                mismatches.append({
                    'sample': sample_id,
                    'ancestry': ancestry,
                    'snp': snp,
                    'raw': raw_value,
                    'combined': combined_value
                })
    
    if mismatches:
        print("Found mismatches:")
        print(pd.DataFrame(mismatches))
    else:
        print("All checked genotypes match!")


# Example usage:
# if __name__ == "__main__":
#     # Example configuration
#     data_dir = '/path/to/data'
#     geno_dir = f'{data_dir}/raw_genotypes'
#     output_dir = f'{data_dir}/outputs'
#     key_file = f'{data_dir}/key.csv'
#     snplist_path = f'{data_dir}/snps_out.csv'
#     labels = ['AAC', 'AFR', 'AJ', 'AMR', 'CAH', 'CAS', 'EAS', 'EUR', 'FIN', 'MDE', 'SAS']
    
#     # Extract carriers for each ancestry label
#     results_by_label = {}
#     for label in labels:
#         results = extract_carriers(
#             geno_path=f'{geno_dir}/{label}/{label}_release9_vwb',
#             snplist_path=snplist_path,
#             out_path=f'{output_dir}/{label}',
#             return_dfs=True  # Get DataFrames directly
#         )
#         results_by_label[label] = results
    
#     # Combine results
#     combined_results = combine_carrier_files(
#         results_by_label=results_by_label,
#         key_file=key_file,
#         output_dir=output_dir
#     )