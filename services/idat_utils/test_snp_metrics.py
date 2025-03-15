import os
from snp_metrics.snp_metrics import process_idat_files

def main():
    # Hardcoded paths
    idat_path = "/home/vitaled2/gp2-microservices/services/idat_utils/data/207847320055"
    bpm = "/home/vitaled2/gp2-microservices/services/idat_utils/data/ilmn_utils/NeuroBooster_20042459_A2.bpm"
    bpm_csv = "/home/vitaled2/gp2-microservices/services/idat_utils/data/ilmn_utils/NeuroBooster_20042459_A2.csv"
    egt = "/home/vitaled2/gp2-microservices/services/idat_utils/data/ilmn_utils/recluster_09092022.egt"
    ref_fasta = "/home/vitaled2/gp2-microservices/services/idat_utils/data/ref/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
    output_directory = "/home/vitaled2/gp2-microservices/services/idat_utils/data/output"
    bcftools_plugins_path = "/home/vitaled2/gp2-microservices/services/idat_utils/bin"
    iaap = "/home/vitaled2/gp2-microservices/services/idat_utils/exec/iaap-cli-linux-x64-1.1.0-sha.80d7e5b3d9c1fdfc2e99b472a90652fd3848bbc7/iaap-cli/iaap-cli"
    
    # Create output directory if it doesn't exist
    os.makedirs(output_directory, exist_ok=True)
    
    print(f"Processing IDAT directory: {idat_path}")
    
    # Process using the updated process_idat_files function
    success = process_idat_files(
        idat_path=idat_path,
        output_directory=output_directory,
        bpm=bpm,
        bpm_csv=bpm_csv,
        egt=egt,
        ref_fasta=ref_fasta,
        iaap=iaap,
        bcftools_plugins_path=bcftools_plugins_path,
        cleanup_intermediate_files=False  # Keep intermediate files for testing
    )
    
    if success:
        print("IDAT processing completed successfully")
        # Print locations of intermediate files for reference
        barcode = os.path.basename(idat_path)
        tmp_dir = os.path.join(output_directory, f"tmp_{barcode}")
        print(f"Intermediate GTC files are stored in: {tmp_dir}")
        print(f"Intermediate VCF files are stored in: {tmp_dir}")
    else:
        print("IDAT processing failed")

if __name__ == "__main__":
    main() 