# test run on single sample
python main.py process  \
--idat_path "/home/vitaled2/gp2-microservices/services/idat_utils/data/207847320055"   \
--output_directory "/home/vitaled2/gp2-microservices/services/idat_utils/data/output/snp_metrics"   \
--bpm "/home/vitaled2/gp2-microservices/services/idat_utils/data/ilmn_utils/NeuroBooster_20042459_A2.bpm"   \
--bpm_csv "/home/vitaled2/gp2-microservices/services/idat_utils/data/ilmn_utils/NeuroBooster_20042459_A2.csv"   \
--egt "/home/vitaled2/gp2-microservices/services/idat_utils/data/ilmn_utils/recluster_09092022.egt"   \
--ref_fasta "/home/vitaled2/gp2-microservices/services/idat_utils/data/ref/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"   \
--iaap "/home/vitaled2/gp2-microservices/services/idat_utils/exec/iaap-cli-linux-x64-1.1.0-sha.80d7e5b3d9c1fdfc2e99b472a90652fd3848bbc7/iaap-cli/iaap-cli"   --bcftools_plugins_path "/home/vitaled2/gp2-microservices/services/idat_utils/bin"

# extract metadata from vcf
python main.py extract  \
--vcf_path "/home/vitaled2/gp2-microservices/services/idat_utils/data/output/snp_metrics/tmp_207847320055/207847320055_R01C01.vcf.gz"   \
--output_path "/home/vitaled2/gp2-microservices/services/idat_utils/data/output/NBA_metadata"   \
--columns "metadata"   \
--output_format "parquet"   \
--partition