import os
import argparse
from snp_metrics import process_idat_files

def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='Process IDAT files to generate SNP metrics.')
    
    # Required arguments
    parser.add_argument('--idat_path', required=True, help='Path to directory containing IDAT files')
    parser.add_argument('--output_directory', required=True, help='Directory to output processed files')
    parser.add_argument('--bpm', required=True, help='Path to BPM file')
    parser.add_argument('--bpm_csv', required=True, help='Path to BPM CSV file')
    parser.add_argument('--egt', required=True, help='Path to EGT file')
    parser.add_argument('--ref_fasta', required=True, help='Path to reference FASTA file')
    parser.add_argument('--iaap', required=True, help='Path to IAAP CLI executable')
    parser.add_argument('--bcftools_plugins_path', required=True, help='Path to bcftools plugins directory')
    
    # Optional arguments
    parser.add_argument('--cleanup', action='store_true', help='Delete intermediate files after processing')
    
    return parser.parse_args()

def main():
    """Run the IDAT processing pipeline."""
    # Parse command line arguments
    args = parse_arguments()
    
    # Create output directory if it doesn't exist
    os.makedirs(args.output_directory, exist_ok=True)
    
    print(f"Processing IDAT directory: {args.idat_path}")
    print(f"Output directory: {args.output_directory}")
    
    success = process_idat_files(
        idat_path=args.idat_path,
        output_directory=args.output_directory,
        bpm=args.bpm,
        bpm_csv=args.bpm_csv,
        egt=args.egt,
        ref_fasta=args.ref_fasta,
        iaap=args.iaap,
        bcftools_plugins_path=args.bcftools_plugins_path,
        cleanup_intermediate_files=args.cleanup
    )
    
    if success:
        barcode = os.path.basename(args.idat_path)
        print(f"Successfully processed IDAT files for {barcode}")
        print(f"Results are available in: {os.path.join(args.output_directory, barcode)}")
        
        if not args.cleanup:
            tmp_dir = os.path.join(args.output_directory, f"tmp_{barcode}")
            print(f"Intermediate files are stored in: {tmp_dir}")
    else:
        print("Failed to process IDAT files")
        return 1
    
    return 0

if __name__ == "__main__":
    exit(main())