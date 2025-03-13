import os
import subprocess
import numpy as np
import pandas as pd
import sys
import psutil
import glob
import gzip
import time
from concurrent.futures import ProcessPoolExecutor
# Suppress copy warning.
pd.options.mode.chained_assignment = None

# def shell_do(command, print_cmd=False, log=False, return_log=False, err=False):
    # if print_cmd:
    #     print(f'Executing: {(" ").join(command.split())}', file=sys.stderr)

    # res=subprocess.run(command.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    # output = res.stdout.decode('utf-8') + res.stderr.decode('utf-8')

    # if log:
    #     print(output)
    # if return_log:
    #     return output
    # if err:
    #     return res.stderr.decode('utf-8')


def cleanup_intermediates(barcode_out_path, clean_up):
    """Clean up intermediate files if requested."""
    if clean_up:
        extensions = ['.vcf', '.gtc']
        for file in os.listdir(barcode_out_path):
            for extension in extensions:
                if file.endswith(extension):
                    os.remove(os.path.join(barcode_out_path, file))


def convert_idat_to_gtc(iaap, bpm, egt, barcode_out_path, idat_path):
    """Convert IDAT files to GTC format."""
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


def gtc_to_vcf(gtc_directory, vcf_directory, bpm, bpm_csv, egt, ref_fasta, out_tmp, bcftools_plugins_path, threads=8, memory="4G"):
    """Convert all GTC files in a directory to multiple VCF files (one per sample)."""
    
    ram_tmp = "/dev/shm/vcf_tmp" if os.path.exists("/dev/shm") else out_tmp
    os.makedirs(ram_tmp, exist_ok=True)
    os.makedirs(vcf_directory, exist_ok=True)
    
    # First, generate a combined VCF in a temporary location
    temp_vcf = os.path.join(out_tmp, "temp_combined.vcf.gz")
    
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
-o {temp_vcf}"""
    
    print(f"Processing all GTC files in: {gtc_directory}")
    result = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    
    if result.returncode != 0:
        print(f"Command failed: {result.stderr.decode('utf-8')}")
        print(f"Error message: {result.stderr.decode('utf-8')}")
        return False
    
    # Index the combined VCF
    index_cmd = f"bcftools index --tbi {temp_vcf}"
    subprocess.run(index_cmd, shell=True)
    
    # Get sample names from VCF
    sample_ids = [x for x in get_vcf_names(temp_vcf) 
                  if x not in ['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT']]
    
    print(f"Splitting combined VCF into {len(sample_ids)} sample files")
    
    # Split the combined VCF by sample
    vcf_files = []
    for sample_id in sample_ids:
        sample_vcf = os.path.join(vcf_directory, f"{sample_id}.vcf.gz")
        split_cmd = f"bcftools view -Oz -s {sample_id} -o {sample_vcf} {temp_vcf}"
        split_result = subprocess.run(split_cmd, shell=True)
        
        if split_result.returncode != 0:
            print(f"Failed to extract sample {sample_id}")
            continue
            
        vcf_files.append(sample_vcf)
    
    # Clean up temporary combined VCF
    os.remove(temp_vcf)
    if os.path.exists(temp_vcf + ".tbi"):
        os.remove(temp_vcf + ".tbi")
    
    print(f"Created {len(vcf_files)} sample VCF files in {vcf_directory}")
    return vcf_files

def get_vcf_names(vcf_path):
    """Get column names from VCF file."""
    opener = gzip.open if vcf_path.endswith('.gz') else open
    mode = 'rt' if vcf_path.endswith('.gz') else 'r'
    
    with opener(vcf_path, mode) as ifile:
            for line in ifile:
                if line.startswith("#CHROM"):
                    vcf_names = [x.strip('\n') for x in line.split('\t')]
                    break
    return vcf_names

def extract_info(info, idx, pattern):
    """Extract information from VCF INFO field."""
    split_info = info.split(";")
    if idx < len(split_info):
        return split_info[idx].replace(pattern, "")
    return None



########################################################
# Process VCF files using chunked reading to minimize memory usage.
########################################################    
# def process_vcf_snps_chunked(vcf_path, out_path=None, chunk_size=50000):
#     """Process VCF files using chunked reading to minimize memory usage."""
    
#     print(f"Processing VCF file: {vcf_path} using chunked reading")
#     start_time = time.time()
    
#     # Get column names
#     names = get_vcf_names(vcf_path)
    
#     # Identify sample columns (not metadata columns)
#     metadata_cols = ['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT']
#     sample_ids = [x for x in names if x not in metadata_cols]
    
#     # Create output directory if needed
#     if out_path:
#         os.makedirs(os.path.dirname(out_path), exist_ok=True)
    
#     # Open file for reading
#     opener = gzip.open if vcf_path.endswith('.gz') else open
#     mode = 'rt' if vcf_path.endswith('.gz') else 'r'
    
#     # Count total lines for progress reporting
#     total_lines = 0
#     with opener(vcf_path, mode) as f:
#         for line in f:
#             if not line.startswith('#'):
#                 total_lines += 1
    
#     print(f"Total data lines to process: {total_lines}")
    
#     # Process file in chunks
#     chunk_count = 0
    
#     with opener(vcf_path, mode) as f:
#         # Skip header lines
#         for line in f:
#             if line.startswith('#') and not line.startswith('#CHROM'):
#                 continue
#             elif line.startswith('#CHROM'):
#                 # Found column headers
#                 break
        
#         # Process data lines in chunks
#         current_chunk = []
#         processed_lines = 0
        
#         for line in f:
#             line_data = line.strip().split('\t')
#             current_chunk.append(line_data)
#             processed_lines += 1
            
#             # Process when chunk reaches desired size or at end of file
#             if len(current_chunk) >= chunk_size:
#                 # Process this chunk
#                 chunk_df = pd.DataFrame(current_chunk, columns=names)
#                 processed_chunk = process_vcf_chunk(chunk_df, sample_ids, metadata_cols)
                
#                 # Write chunk to parquet
#                 if out_path:
#                     chunk_file = f"{out_path}_chunk_{chunk_count}.parquet"
#                     processed_chunk.to_parquet(chunk_file)
                
#                 # Report progress
#                 progress = (processed_lines / total_lines) * 100
#                 elapsed = time.time() - start_time
#                 print(f"Progress: {progress:.2f}% ({processed_lines}/{total_lines} lines) - {elapsed:.2f} seconds elapsed")
                
#                 # Clear memory
#                 del chunk_df, processed_chunk, current_chunk
#                 current_chunk = []
#                 chunk_count += 1
        
#         # Process the last chunk if any
#         if current_chunk:
#             chunk_df = pd.DataFrame(current_chunk, columns=names)
#             processed_chunk = process_vcf_chunk(chunk_df, sample_ids, metadata_cols)
            
#             if out_path:
#                 chunk_file = f"{out_path}_chunk_{chunk_count}.parquet"
#                 processed_chunk.to_parquet(chunk_file)
    
#     print(f"Processed {chunk_count + 1} chunks in {time.time() - start_time:.2f} seconds")
    
#     # Merge all chunk files
#     if out_path:
#         merge_parquet_chunks(f"{out_path}_chunk_*.parquet", f"{out_path}_merged")
    
#     return True

# def process_vcf_chunk(chunk_df, sample_ids, metadata_cols):
#     """Process a single chunk of VCF data."""
    
#     # Fix metadata columns to match actual dataframe columns
#     fixed_metadata_cols = []
#     for col in metadata_cols:
#         if col == '#CHROM' and '#CHROM' not in chunk_df.columns:
#             if 'CHROM' in chunk_df.columns:
#                 fixed_metadata_cols.append('CHROM')
#             else:
#                 # Look for a column that might be CHROM (like chr1)
#                 chrom_cols = [c for c in chunk_df.columns if c.startswith('chr')]
#                 if chrom_cols:
#                     fixed_metadata_cols.append(chrom_cols[0])
#         else:
#             fixed_metadata_cols.append(col)
    
#     # Melt dataframe to convert from wide to long format
#     vcf_melt = pd.melt(
#         chunk_df, 
#         id_vars=fixed_metadata_cols, 
#         value_vars=sample_ids,
#         var_name='sampleid', 
#         value_name='metrics'
#     )
    
#     # Get the chromosome column name
#     chrom_col = 'CHROM'
#     if '#CHROM' in vcf_melt.columns:
#         chrom_col = '#CHROM'
#     elif any(col.startswith('chr') for col in vcf_melt.columns):
#         potential_chrom_cols = [col for col in vcf_melt.columns if col.startswith('chr')]
#         if potential_chrom_cols:
#             chrom_col = potential_chrom_cols[0]
    
#     # If the column isn't CHROM or #CHROM, rename it
#     if chrom_col != 'CHROM' and chrom_col != '#CHROM':
#         vcf_melt = vcf_melt.rename(columns={chrom_col: 'CHROM'})
#     elif chrom_col == '#CHROM':
#         vcf_melt = vcf_melt.rename(columns={'#CHROM': 'CHROM'})
    
#     # Split metrics column
#     metric_cols = ['GT','GQ','IGC','BAF','LRR','NORMX','NORMY','R','THETA','X','Y']
#     vcf_melt[metric_cols] = vcf_melt['metrics'].str.split(':', expand=True, n=10)
    
#     # Extract information from INFO column
#     vcf_melt['GenTrain_Score'] = vcf_melt['INFO'].apply(
#         lambda info: extract_info(info, idx=10, pattern='GenTrain_Score=')
#     )
#     vcf_melt['ALLELE_A'] = vcf_melt['INFO'].apply(
#         lambda info: extract_info(info, idx=1, pattern='ALLELE_A=')
#     )
#     vcf_melt['ALLELE_B'] = vcf_melt['INFO'].apply(
#         lambda info: extract_info(info, idx=2, pattern='ALLELE_B=')
#     )
    
#     # Drop unused columns
#     vcf_melt = vcf_melt.drop(columns=[
#         'QUAL', 'FILTER', 'INFO', 'GQ', 'IGC', 
#         'NORMX', 'NORMY', 'X', 'Y', 'metrics', 'FORMAT'
#     ])
    
#     # Clean up chromosome column
#     vcf_melt['CHROM'] = vcf_melt['CHROM'].astype(str).str.replace('chr','')
    
#     # Map GT values
#     gtype_map = {'0/0':'AA', '0/1':'AB', '1/1':'BB', './.':'NC'}
#     vcf_melt['GType'] = vcf_melt['GT'].map(gtype_map)
    
#     # Process GT values
#     vcf_melt['GT'] = vcf_melt['GType']
    
#     # Apply update_gt function
#     vcf_melt.loc[(vcf_melt['GType'] == 'AA') & (vcf_melt['ALLELE_A'].astype(str) == '1'), 'GT'] = 'BB'
#     vcf_melt.loc[(vcf_melt['GType'] == 'AA') & (vcf_melt['ALLELE_A'].astype(str) == '0'), 'GT'] = 'AA'
#     vcf_melt.loc[(vcf_melt['GType'] == 'BB') & (vcf_melt['ALLELE_B'].astype(str) == '1'), 'GT'] = 'BB'
#     vcf_melt.loc[(vcf_melt['GType'] == 'BB') & (vcf_melt['ALLELE_B'].astype(str) == '0'), 'GT'] = 'AA'
#     vcf_melt.loc[(vcf_melt['GType'] == 'AB'), 'GT'] = 'AB'
#     vcf_melt.loc[(vcf_melt['GType'] == 'NC'), 'GT'] = 'NC'
#     vcf_melt['GT'] = vcf_melt['GT'].fillna('NC')
    
#     # Process alternate alleles
#     vcf_melt['a1'] = vcf_melt['REF']
#     vcf_melt.loc[vcf_melt['ALLELE_A'].astype(str) == '1', 'a1'] = vcf_melt.loc[vcf_melt['ALLELE_A'].astype(str) == '1', 'ALT']
    
#     vcf_melt['a2'] = vcf_melt['REF']
#     vcf_melt.loc[vcf_melt['ALLELE_B'].astype(str) == '1', 'a2'] = vcf_melt.loc[vcf_melt['ALLELE_B'].astype(str) == '1', 'ALT']
    
#     # Rename columns to match expected output format
#     final_df = vcf_melt.rename(columns={
#         'CHROM': 'chromosome',
#         'POS': 'position',
#         'ID': 'snpID',
#         'sampleid': 'Sample_ID',
#         'REF': 'Ref',
#         'ALT': 'Alt'
#     })
    
#     # Convert types
#     final_df = final_df.astype({
#         'chromosome': str,
#         'position': int,
#         'snpID': str,
#         'Sample_ID': str,
#         'Ref': str,
#         'Alt': str,
#         'ALLELE_A': int,
#         'ALLELE_B': int,
#         'BAF': float,
#         'LRR': float,
#         'R': float,
#         'THETA': float,
#         'GenTrain_Score': float,
#         'GType': str
#     })
    
#     # Select and order columns
#     out_colnames = [
#         'chromosome', 'position', 'snpID', 'Sample_ID', 'Ref', 'Alt',
#         'ALLELE_A', 'ALLELE_B', 'BAF', 'LRR', 'R', 'THETA', 
#         'GenTrain_Score', 'GType', 'GT', 'a1', 'a2'
#     ]
    
#     return final_df[out_colnames]

########################################################
# end of chunked processing
########################################################


def merge_parquet_chunks(chunk_pattern, output_directory):
    """Merge multiple parquet chunks into a single parquet file."""
    import glob
    import pandas as pd
    import os
    
    chunk_files = sorted(glob.glob(chunk_pattern))
    if not chunk_files:
        print(f"No chunk files found matching pattern: {chunk_pattern}")
        return False
    
    print(f"Merging {len(chunk_files)} chunk files...")
    
    # Create a directory for the merged dataset
    os.makedirs(output_directory, exist_ok=True)
    
    # Process each chunk file individually
    for i, chunk_file in enumerate(chunk_files):
        print(f"Processing chunk {i+1}/{len(chunk_files)}: {chunk_file}")
        
        try:
            # Read chunk
            chunk_df = pd.read_parquet(chunk_file)
            
            # Always partition by chromosome
            chunk_df.to_parquet(
                output_directory,
                partition_cols=['chromosome'],
                compression='brotli'
            )
            
            # Clean up memory
            del chunk_df
            
            # Remove processed chunk file to save space
            os.remove(chunk_file)
            
        except Exception as e:
            print(f"Error processing {chunk_file}: {e}")
    
    print(f"All chunks merged into: {output_directory}")
    return True


def process_idat_files(idat_path, output_directory, bpm, bpm_csv, egt, ref_fasta, iaap, bcftools_plugins_path):
    """Process a single IDAT directory to generate SNP metrics."""
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
    
    # Optimize resource allocation based on system specs
    cpu_count = os.cpu_count()
    total_memory_gb = psutil.virtual_memory().total / (1024**3)
    
    # Allocate resources for different processing stages
    iaap_threads = min(8, cpu_count)
    bcftools_threads = max(1, cpu_count - 1)
    
    print(f"Resource allocation:")
    print(f"- IAAP threads: {iaap_threads}")
    print(f"- bcftools threads: {bcftools_threads}")
    
    # Step 1: Convert IDAT files to GTC
    # print(f"Converting IDAT to GTC for {barcode}...")
    # gtc_files = convert_idat_to_gtc(iaap, bpm, egt, barcode_out_path, idat_path)
    # if not gtc_files:
    #     print(f"Failed to convert IDAT to GTC for {barcode}")
    #     return False
    
    # print(f"Successfully created {len(gtc_files)} GTC files")
    
    # Step 2: Process all GTC files and create per-sample VCF files
    vcf_directory = os.path.join(barcode_out_path, "sample_vcfs")
    # print(f"Converting GTC files to per-sample VCF files...")
    # vcf_files = gtc_to_vcf(barcode_out_path, vcf_directory, bpm, bpm_csv, egt, ref_fasta, out_tmp, 
    #                  bcftools_plugins_path, threads=bcftools_threads, 
    #                  memory=f"{int(total_memory_gb * 0.7)}G")  # Use 70% of available memory
    
    # if not vcf_files:
    #     print(f"Failed to convert GTC files to VCF for {barcode}")
    #     return False
    
    # Step 3: Process each sample VCF file in parallel
    vcf_files = sorted(glob.glob(os.path.join(vcf_directory, '*.vcf.gz')))
    print(f"Processing {len(vcf_files)} sample VCF files for {barcode}...")
    
    # Set chunk size based on available memory (can be higher for single sample files)
    if total_memory_gb < 16:
        chunk_size = 50000  # Doubled from before since we're only processing one sample at a time
    elif total_memory_gb < 32:
        chunk_size = 150000
    else:
        chunk_size = 200000
    
    success = True
    with ProcessPoolExecutor(max_workers=min(cpu_count, len(vcf_files))) as executor:
        futures = []
        
        for vcf_file in vcf_files:
            sample_id = os.path.basename(vcf_file).replace('.vcf.gz', '')
            # Change output path to be per-sample directory
            sample_output_dir = os.path.join(output_directory, sample_id)
            os.makedirs(sample_output_dir, exist_ok=True)
            
            # Use temporary directory for chunks
            temp_output_path = os.path.join(out_tmp, f"temp_{sample_id}")
            
            # Submit processing task to executor
            future = executor.submit(
                process_single_sample_vcf, 
                vcf_file,
                temp_output_path,
                chunk_size,
                sample_output_dir  # Pass final output directory
            )
            futures.append((sample_id, future))
        
        # Wait for all tasks to complete and check for errors
        for sample_id, future in futures:
            try:
                if not future.result():
                    print(f"Failed to process SNP data for sample {sample_id}")
                    success = False
            except Exception as e:
                print(f"Error processing sample {sample_id}: {e}")
                import traceback
                traceback.print_exc()
                success = False
    
    if not success:
        print(f"Some samples failed to process for {barcode}")
        return False
    
    print(f"Processing complete for {barcode}")
    return True


def process_single_sample_vcf(vcf_file, out_path, chunk_size=100000, final_output_dir=None):
    """Process a single-sample VCF file."""
    print(f"Processing sample VCF file: {vcf_file}")
    start_time = time.time()
    
    # Get sample ID from filename
    sample_id = os.path.basename(vcf_file).replace('.vcf.gz', '')
    
    # Create output directory if needed
    if out_path:
        os.makedirs(os.path.dirname(out_path), exist_ok=True)
    
    # Open file for reading
    opener = gzip.open if vcf_file.endswith('.gz') else open
    mode = 'rt' if vcf_file.endswith('.gz') else 'r'
    
    # Get column names
    names = get_vcf_names(vcf_file)
    
    # Ensure we have the right sample column
    metadata_cols = ['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT']
    sample_ids = [x for x in names if x not in metadata_cols]
    
    if len(sample_ids) != 1:
        print(f"Expected single sample VCF, but found {len(sample_ids)} samples: {sample_ids}")
        return False
    
    # Count total lines for progress reporting
    total_lines = 0
    with opener(vcf_file, mode) as f:
        for line in f:
            if not line.startswith('#'):
                total_lines += 1
    
    print(f"Total data lines to process: {total_lines}")
    
    # Process file in chunks
    chunk_count = 0
    
    with opener(vcf_file, mode) as f:
        # Skip header lines
        for line in f:
            if line.startswith('#') and not line.startswith('#CHROM'):
                continue
            elif line.startswith('#CHROM'):
                # Found column headers
                break
        
        # Process data lines in chunks
        current_chunk = []
        processed_lines = 0
        
        for line in f:
            line_data = line.strip().split('\t')
            current_chunk.append(line_data)
            processed_lines += 1
            
            # Process when chunk reaches desired size or at end of file
            if len(current_chunk) >= chunk_size:
                # Process this chunk
                chunk_df = pd.DataFrame(current_chunk, columns=names)
                
                # For a single sample VCF, we don't need to melt - just extract the data directly
                processed_chunk = process_single_sample_chunk(chunk_df, sample_ids[0])
                
                # Write chunk to parquet
                if out_path:
                    chunk_file = f"{out_path}_chunk_{chunk_count}.parquet"
                    processed_chunk.to_parquet(chunk_file)
                
                # Report progress
                progress = (processed_lines / total_lines) * 100
                elapsed = time.time() - start_time
                print(f"Sample {sample_id} - Progress: {progress:.2f}% ({processed_lines}/{total_lines} lines) - {elapsed:.2f} seconds elapsed")
                
                # Clear memory
                del chunk_df, processed_chunk, current_chunk
                current_chunk = []
                chunk_count += 1
        
        # Process the last chunk if any
        if current_chunk:
            chunk_df = pd.DataFrame(current_chunk, columns=names)
            processed_chunk = process_single_sample_chunk(chunk_df, sample_ids[0])
            
            if out_path:
                chunk_file = f"{out_path}_chunk_{chunk_count}.parquet"
                processed_chunk.to_parquet(chunk_file)
    
    print(f"Sample {sample_id} - Processed {chunk_count + 1} chunks in {time.time() - start_time:.2f} seconds")
    
    # Merge all chunk files
    if out_path:
        # If final_output_dir is specified, use that instead of the default merged location
        output_dir = final_output_dir if final_output_dir else f"{out_path}_merged"
        merge_parquet_chunks(f"{out_path}_chunk_*.parquet", output_dir)
    
    return True

def process_single_sample_chunk(chunk_df, sample_id):
    """Process a single chunk of a single-sample VCF data."""
    
    # Fix metadata columns to match actual dataframe columns
    chrom_col = '#CHROM'
    if '#CHROM' not in chunk_df.columns and 'CHROM' in chunk_df.columns:
        chrom_col = 'CHROM'
    elif '#CHROM' not in chunk_df.columns:
        # Look for a column that might be CHROM (like chr1)
        chrom_cols = [c for c in chunk_df.columns if c.startswith('chr')]
        if chrom_cols:
            chrom_col = chrom_cols[0]
    
    # Extract data for this sample - no need to melt since we only have one sample
    sample_data = chunk_df[[chrom_col, 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', sample_id]].copy()
    
    # Rename CHROM column if needed
    if chrom_col != 'CHROM':
        sample_data = sample_data.rename(columns={chrom_col: 'CHROM'})
    
    # Split metrics column
    metric_cols = ['GT','GQ','IGC','BAF','LRR','NORMX','NORMY','R','THETA','X','Y']
    sample_data[metric_cols] = sample_data[sample_id].str.split(':', expand=True, n=10)
    
    # Extract information from INFO column
    sample_data['GenTrain_Score'] = sample_data['INFO'].apply(
        lambda info: extract_info(info, idx=10, pattern='GenTrain_Score=')
    )
    sample_data['ALLELE_A'] = sample_data['INFO'].apply(
        lambda info: extract_info(info, idx=1, pattern='ALLELE_A=')
    )
    sample_data['ALLELE_B'] = sample_data['INFO'].apply(
        lambda info: extract_info(info, idx=2, pattern='ALLELE_B=')
    )
    
    # Drop unused columns
    sample_data = sample_data.drop(columns=[
        'QUAL', 'FILTER', 'INFO', 'GQ', 'IGC', 
        'NORMX', 'NORMY', 'X', 'Y', sample_id, 'FORMAT'
    ])
    
    # Clean up chromosome column
    sample_data['CHROM'] = sample_data['CHROM'].astype(str).str.replace('chr','')
    
    # Map GT values
    gtype_map = {'0/0':'AA', '0/1':'AB', '1/1':'BB', './.':'NC'}
    sample_data['GType'] = sample_data['GT'].map(gtype_map)
    
    # Store current GT values and then drop the original GT column
    orig_gt = sample_data['GT'].copy()  # Save original values
    sample_data = sample_data.drop(columns=['GT'])  # Drop original GT column
    
    # Process GT values (now starting fresh)
    sample_data['GT'] = sample_data['GType']  # Start with GType values
    
    # Apply the transformations directly to the GT column
    sample_data.loc[(sample_data['GType'] == 'AA') & (sample_data['ALLELE_A'].astype(str) == '1'), 'GT'] = 'BB'
    sample_data.loc[(sample_data['GType'] == 'AA') & (sample_data['ALLELE_A'].astype(str) == '0'), 'GT'] = 'AA'
    sample_data.loc[(sample_data['GType'] == 'BB') & (sample_data['ALLELE_B'].astype(str) == '1'), 'GT'] = 'BB'
    sample_data.loc[(sample_data['GType'] == 'BB') & (sample_data['ALLELE_B'].astype(str) == '0'), 'GT'] = 'AA'
    sample_data.loc[(sample_data['GType'] == 'AB'), 'GT'] = 'AB'
    sample_data.loc[(sample_data['GType'] == 'NC'), 'GT'] = 'NC'
    sample_data['GT'] = sample_data['GT'].fillna('NC')
    
    # Process alternate alleles
    sample_data['a1'] = sample_data['REF']
    sample_data.loc[sample_data['ALLELE_A'].astype(str) == '1', 'a1'] = sample_data.loc[sample_data['ALLELE_A'].astype(str) == '1', 'ALT']
    
    sample_data['a2'] = sample_data['REF']
    sample_data.loc[sample_data['ALLELE_B'].astype(str) == '1', 'a2'] = sample_data.loc[sample_data['ALLELE_B'].astype(str) == '1', 'ALT']
    
    # Add sample ID column
    sample_data['Sample_ID'] = sample_id
    
    # Rename columns to match expected output format
    final_df = sample_data.rename(columns={
        'CHROM': 'chromosome',
        'POS': 'position',
        'ID': 'snpID',
        'REF': 'Ref',
        'ALT': 'Alt'
    })
    
    # Convert types
    final_df = final_df.astype({
        'chromosome': str,
        'position': int,
        'snpID': str,
        'Sample_ID': str,
        'Ref': str,
        'Alt': str,
        'ALLELE_A': int,
        'ALLELE_B': int,
        'BAF': float,
        'LRR': float,
        'R': float,
        'THETA': float,
        'GenTrain_Score': float,
        'GType': str,
        'GT': str
    })
    
    # Select and order columns
    out_colnames = [
        'chromosome', 'position', 'snpID', 'Sample_ID', 'Ref', 'Alt',
        'ALLELE_A', 'ALLELE_B', 'BAF', 'LRR', 'R', 'THETA', 
        'GenTrain_Score', 'GType', 'GT', 'a1', 'a2'
    ]
    
    return final_df[out_colnames]

def extract_vcf_columns(vcf_file, output_path=None, num_rows=10, columns="all"):
    """
    Extract rows and specific columns from a VCF file, split INFO and FORMAT fields.
    
    Args:
        vcf_file: Path to VCF file
        output_path: Optional path to save the extracted data to parquet format
        num_rows: Number of rows to extract (default: 10, use None for all rows)
        columns: Columns to extract. Options:
            - "all": Extract all columns (default)
            - "metadata": Extract only metadata columns
            - "sample": Extract only sample-specific columns
            - List of specific column names to extract
        
    Returns:
        DataFrame containing the extracted data with selected columns
    """
    # Define metadata columns
    metadata_cols_list = [
        'CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER',
        'ASSAY_TYPE', 'devR_AB', 'FRAC_T', 'FRAC_G', 'meanTHETA_BB', 
        'meanR_AB', 'devTHETA_AB', 'GC', 'N_AA', 'Orig_Score', 
        'FRAC_C', 'GenTrain_Score', 'devR_BB', 'NORM_ID', 
        'devR_AA', 'Intensity_Threshold', 'meanR_AA', 'devTHETA_AA', 
        'ALLELE_A', 'N_AB', 'meanR_BB', 'meanTHETA_AA', 
        'meanTHETA_AB', 'devTHETA_BB', 'N_BB', 'ALLELE_B', 
        'FRAC_A', 'BEADSET_ID', 'Cluster_Sep'
    ]
    
    # Define sample-specific columns
    sample_specific_list = ['Sample_ID', 'GT', 'GQ', 'IGC', 'BAF', 'LRR', 'NORMX', 'NORMY', 'R', 'THETA', 'X', 'Y']
    
    # Print extraction info
    if num_rows is None:
        print(f"Extracting ALL rows from VCF: {vcf_file}")
    else:
        print(f"Extracting first {num_rows} rows from VCF: {vcf_file}")
    
    start_time = time.time()
    
    # Open file for reading
    opener = gzip.open if vcf_file.endswith('.gz') else open
    mode = 'rt' if vcf_file.endswith('.gz') else 'r'
    
    # Get column names from VCF
    vcf_names = get_vcf_names(vcf_file)
    vcf_metadata_cols = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT']
    vcf_sample_cols = [x for x in vcf_names if x not in vcf_metadata_cols]
    
    # Read rows
    current_chunk = []
    row_count = 0
    
    with opener(vcf_file, mode) as f:
        # Skip header lines
        for line in f:
            if line.startswith('#') and not line.startswith('#CHROM'):
                continue
            elif line.startswith('#CHROM'):
                # Found column headers
                break
        
        # Read data lines
        for line in f:
            # Break if we've reached the requested number of rows
            if num_rows is not None and row_count >= num_rows:
                break
                
            line_data = line.strip().split('\t')
            current_chunk.append(line_data)
            row_count += 1
            
            # For large files, print progress periodically
            if num_rows is None and row_count % 100000 == 0:
                print(f"Read {row_count} rows so far...")
    
    # Create DataFrame
    if not current_chunk:
        print("No data rows found in VCF file")
        return pd.DataFrame()
        
    # Create initial DataFrame with all columns
    result_df = pd.DataFrame(current_chunk, columns=vcf_names)
    print(f"Read {len(result_df)} total rows")
    
    # Fix CHROM column name if needed
    chrom_col = '#CHROM'
    if '#CHROM' not in result_df.columns and 'CHROM' in result_df.columns:
        chrom_col = 'CHROM'
    elif '#CHROM' not in result_df.columns:
        chrom_cols = [c for c in result_df.columns if c.startswith('chr')]
        if chrom_cols:
            chrom_col = chrom_cols[0]
    
    # Rename chromosome column for consistency
    if chrom_col != 'CHROM':
        result_df = result_df.rename(columns={chrom_col: 'CHROM'})
    
    # Clean up chromosome column
    result_df['CHROM'] = result_df['CHROM'].astype(str).str.replace('chr', '')
    
    # Convert numeric columns to appropriate types
    if 'POS' in result_df.columns:
        result_df['POS'] = result_df['POS'].astype(int)
    
    # Keep track of columns to drop
    columns_to_drop = []
    
    # Extract sample ID from the filename
    sample_id = os.path.basename(vcf_file).replace('.vcf.gz', '')
    
    # Step 1: Parse INFO column into separate columns
    if 'INFO' in result_df.columns:
        print("Parsing INFO column into separate fields...")
        # Get all INFO keys first
        info_keys = set()
        for info_str in result_df['INFO']:
            for field in info_str.split(';'):
                if '=' in field:
                    key = field.split('=')[0]
                    info_keys.add(key)
        
        # Create a column for each INFO key - direct column name without prefix
        for key in info_keys:
            result_df[key] = result_df['INFO'].apply(
                lambda info_str: next(
                    (item.split('=')[1] for item in info_str.split(';') 
                     if item.startswith(f"{key}=")), 
                    None
                )
            )
        
        # Mark INFO column for deletion
        columns_to_drop.append('INFO')
    
    # Step 2: Parse FORMAT column and sample genotype data
    if 'FORMAT' in result_df.columns and vcf_sample_cols:
        print("Parsing FORMAT column and sample data...")
        # Get all FORMAT fields from the first row (usually consistent)
        format_fields = result_df['FORMAT'].iloc[0].split(':')
        
        for sample in vcf_sample_cols:
            # Check if the sample column exists before processing
            if sample in result_df.columns:
                # Split the sample data by ':'
                sample_data = result_df[sample].str.split(':', expand=True)
                
                # Create new columns for each format field - without sample prefix
                for i, field in enumerate(format_fields):
                    if i < sample_data.shape[1]:  # Only process if there's data
                        # Use just the field name - if there are multiple samples, 
                        # the last one will overwrite previous ones
                        result_df[field] = sample_data[i]
                
                # Mark this sample column for deletion
                columns_to_drop.append(sample)
        
        # Mark FORMAT column for deletion
        columns_to_drop.append('FORMAT')
    
    # Drop all the columns we've processed
    if columns_to_drop:
        result_df = result_df.drop(columns=columns_to_drop)
    
    # Add Sample_ID column if it doesn't exist yet
    if 'Sample_ID' not in result_df.columns:
        result_df['Sample_ID'] = sample_id
    
    # Add IID column (same as Sample_ID for now)
    result_df['IID'] = result_df['Sample_ID']
    
    # Ensure 'ID' is in the columns for all filtering options
    if columns == "all":
        # Keep all columns
        filtered_df = result_df
    elif columns == "metadata":
        # Filter to only metadata columns, ensure ID is included
        available_metadata = [col for col in metadata_cols_list if col in result_df.columns]
        if 'ID' not in available_metadata and 'ID' in result_df.columns:
            available_metadata.append('ID')
        filtered_df = result_df[available_metadata]
    elif columns == "sample":
        # Filter to only sample-specific columns, ensure ID and IID are included
        available_sample_cols = [col for col in sample_specific_list if col in result_df.columns]
        required_cols = ['ID', 'IID']
        for col in required_cols:
            if col not in available_sample_cols and col in result_df.columns:
                available_sample_cols.append(col)
        filtered_df = result_df[available_sample_cols]
    elif isinstance(columns, list):
        # Filter to user-specified columns, ensure ID is included
        available_columns = [col for col in columns if col in result_df.columns]
        
        # Add ID if not already in the list
        if 'ID' not in available_columns and 'ID' in result_df.columns:
            available_columns.append('ID')
        
        # If any sample-specific column is requested, include IID
        has_sample_column = any(col in sample_specific_list for col in available_columns)
        if has_sample_column and 'IID' not in available_columns and 'IID' in result_df.columns:
            available_columns.append('IID')
            
        filtered_df = result_df[available_columns]
    else:
        print(f"Warning: Unrecognized columns option '{columns}'. Returning all columns.")
        filtered_df = result_df
    
    print(f"Extracted and processed {len(filtered_df)} rows with {len(filtered_df.columns)} columns in {time.time() - start_time:.2f} seconds")
    
    # Save result if a path is provided
    if output_path:
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        filtered_df.to_csv(output_path, index=False)
        print(f"Saved VCF data to: {output_path}")
    
    return filtered_df






    
