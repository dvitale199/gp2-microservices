import os
from snp_metrics.snp_metrics import idat_snp_metrics

def main():
    # Hardcoded paths
    idat_path = "/home/vitaled2/gp2-microservices/services/idat_utils/data/207847320055"
    bpm = "/home/vitaled2/gp2-microservices/services/idat_utils/data/ilmn_utils/NeuroBooster_20042459_A2.bpm"
    bpm_csv = "/home/vitaled2/gp2-microservices/services/idat_utils/data/ilmn_utils/NeuroBooster_20042459_A2.csv"
    egt = "/home/vitaled2/gp2-microservices/services/idat_utils/data/ilmn_utils/recluster_09092022.egt"
    ref_fasta = "/home/vitaled2/gp2-microservices/services/idat_utils/data/ref/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
    out_path = "/home/vitaled2/gp2-microservices/services/idat_utils/data/output"
    iaap = "/home/vitaled2/gp2-microservices/services/idat_utils/exec/iaap-cli-linux-x64-1.1.0-sha.80d7e5b3d9c1fdfc2e99b472a90652fd3848bbc7/iaap-cli/iaap-cli"
    
    # Processing options
    # num_threads = 8
    # clean_up = True
    # verbose = True
    
    # Set environment variables
    # os.environ['SNP_METRICS_NUM_THREADS'] = str(num_threads)
    # if verbose:
    #     os.environ['SNP_METRICS_VERBOSE'] = 'True'
    
    # print(f"Processing barcode with parallel samples: {idat_path}")
    # print(f"Using {num_threads} threads")
    
    # Process the barcode
    idat_snp_metrics(
        idat_path=idat_path,
        bpm=bpm,
        bpm_csv=bpm_csv,
        egt=egt,
        ref_fasta=ref_fasta,
        out_path=out_path,
        iaap=iaap,
        clean_up=True
    )
    
    print("Barcode processing complete")

if __name__ == "__main__":
    main() 