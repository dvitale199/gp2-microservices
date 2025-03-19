sudo apt install wget unzip git g++ zlib1g-dev bwa unzip samtools vcftools msitools cabextract mono-devel libgdiplus icu-devtools bcftools

cd services/idat_utils/data/ref/

wget -O- ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz | \
  gzip -d > GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
samtools faidx GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
bwa index GCA_000001405.15_GRCh38_no_alt_analysis_set.fna