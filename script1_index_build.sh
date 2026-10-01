#!/bin/bash
#SBATCH --account=XXXXX
#SBATCH --partition=default_free
#SBATCH --mem=50G
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=40
#SBATCH --job-name=dna_index_prep
#SBATCH --output=slurm_log_%j.out

# Script to prepare DNA panel sequencing mutation analysis environment including downloading references and build BWA indexes
# Ruth Cranston 2026

[ $# -ne 0 ] && { echo -en \
"\nRuth Cranston 2026\n\n
*** Script to build DNA panel sequencing mutation analysis environment including downloading references and build BWA indexes.
Prepared both GRCh38 and GRCh37
Runs in current directory *** \n\n" ; exit 1; }

# Set variables
BASE_DIR="$PWD"
REFERENCE_DIR=${BASE_DIR}/"References"
REF38=${REFERENCE_DIR}/GRCh38
REF37=${REFERENCE_DIR}/GRCh37
BWA_INDEX_DIR=${BASE_DIR}/"BWA_indexes"

# Load modules
echo -en "Loading modules...\n"
module --force purge
module load BWA/0.7.18-GCCcore-13.3.0
module load BEDTools
module load GATK/4.6.0.0-GCCcore-13.2.0-Java-17
module load VEP/113.3-GCC-13.3.0
module load Python/3.11.5-GCCcore-13.2.0

set -euo pipefail

echo -en "Environment set up.\n"

# make directories
mkdir -p ${REF38}
mkdir -p ${REF37}
mkdir -p ${BWA_INDEX_DIR}
mkdir -p logs

# Load Python >=3.10 for gcloud CLI compatibility
export CLOUDSDK_PYTHON=$(which python3)

# install google cloud sdk locally for download of files from google repo
GCLOUD_DIR="$(pwd)/google-cloud-sdk"

if [ ! -d "${GCLOUD_DIR}" ]; then
    echo "Installing Google Cloud SDK..."
    curl -s https://sdk.cloud.google.com | bash -s -- --disable-prompts --install-dir="$(pwd)"
fi

# Add to PATH for this script's session
export PATH="${GCLOUD_DIR}/bin:${PATH}"

# Confirm it's available
which gsutil
gsutil --version

# Configure anonymous access for public buckets
gcloud config set auth/disable_credentials True
gcloud config unset project 2>/dev/null || true

######################
### GRCh38 section ###
######################

# checking if reference files are present - if not then download these
echo "Detecting if GRCh38 references are present"
if [[ -f ${REF38}/Homo_sapiens_assembly38.fasta ]];
then

    echo -en " * GRCh38 references already exist in ${REF38}, no need to re-download\n\n"

else

    echo -en " * Downloading GRCh38 references to ${REF38} now...\n"

    # downloads: gatk specific genome build, indexes and dict files, gencodev44 gtf files, BQSR files, Mutect2 resources
    # gatk genome files
    gsutil cp gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta ${REF38}/
    gsutil cp gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai ${REF38}/
    gsutil cp gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dict ${REF38}/

    # BQSR files
    curl -L "https://storage.googleapis.com/storage/v1/b/gcp-public-data--broad-references/o/hg38%2Fv0%2FHomo_sapiens_assembly38.dbsnp138.vcf?alt=media" -o ${REF38}/Homo_sapiens_assembly38.dbsnp138.vcf
    curl -L "https://storage.googleapis.com/storage/v1/b/gcp-public-data--broad-references/o/hg38%2Fv0%2FHomo_sapiens_assembly38.dbsnp138.vcf.idx?alt=media" -o ${REF38}/Homo_sapiens_assembly38.dbsnp138.vcf.idx
    gsutil cp gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.known_indels.vcf.gz ${REF38}/
    gsutil cp gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.known_indels.vcf.gz.tbi ${REF38}/
    gsutil cp gs://gcp-public-data--broad-references/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz ${REF38}/
    gsutil cp gs://gcp-public-data--broad-references/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi ${REF38}/

    # Mutect2 files
    # NB: this PoN is WGS-derived. For panel data, consider building a matched PoN from your own normals for best filtering.
    gsutil cp gs://gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz "${REF38}/"
    gsutil cp gs://gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz.tbi "${REF38}/"
    gsutil cp gs://gatk-best-practices/somatic-hg38/1000g_pon.hg38.vcf.gz "${REF38}/"
    gsutil cp gs://gatk-best-practices/somatic-hg38/1000g_pon.hg38.vcf.gz.tbi "${REF38}/"
    gsutil cp gs://gatk-best-practices/somatic-hg38/small_exac_common_3.hg38.vcf.gz "${REF38}/"
    gsutil cp gs://gatk-best-practices/somatic-hg38/small_exac_common_3.hg38.vcf.gz.tbi "${REF38}/"

    # gencode gtf
    wget "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.primary_assembly.annotation.gtf.gz" -O "${REF38}/gencode.v44.primary_assembly.annotation.gtf.gz"
    gunzip "${REF38}/gencode.v44.primary_assembly.annotation.gtf.gz"

    # make chromosome sizes file
    cut -f1,2 ${REF38}/Homo_sapiens_assembly38.fasta.fai > ${REF38}/GRCh38.genome
    echo "GRCh38 reference download complete."

    # Download files for vep cache
    mkdir -p ${REF38}/vep_cache
    wget https://ftp.ensembl.org/pub/release-113/variation/indexed_vep_cache/homo_sapiens_vep_113_GRCh38.tar.gz \
         -O ${REF38}/vep_cache/homo_sapiens_vep_113_GRCh38.tar.gz
    tar xzf ${REF38}/vep_cache/homo_sapiens_vep_113_GRCh38.tar.gz -C ${REF38}/vep_cache

    # BWA reference genome preparation:
    echo "Detecting if BWA indexes are present"
    if [[ -d ${BWA_INDEX_DIR}/BWA_GRCh38 ]]
    then
        echo -en " * ${BWA_INDEX_DIR}/BWA_GRCh38 exists, no need to re-create\n\n"
    else
        echo -en " * ${BWA_INDEX_DIR}/BWA_GRCh38 not found, preparing BWA index...\n"
        mkdir ${BWA_INDEX_DIR}/BWA_GRCh38
        # -p sets the index prefix; index files are written here while the fasta stays in REF38 (no duplication)
        bwa index -p ${BWA_INDEX_DIR}/BWA_GRCh38/Homo_sapiens_assembly38 ${REF38}/Homo_sapiens_assembly38.fasta
    fi

fi

######################
### GRCh37 section ###
######################

echo "Detecting if GRCh37 references are present"
if [[ -f ${REF37}/Homo_sapiens_assembly19.fasta ]];
then

    echo -en " * GRCh37 references already exist in ${REF37}, no need to re-download\n\n"

else

    echo -en " * Downloading GRCh37 references to ${REF37} now...\n"
    # downloads: gatk specific genome build, indexes and dict files, Ensembl release 87 gtf, BQSR files, Mutect2 resources

    # gatk genome files
    gsutil cp gs://gcp-public-data--broad-references/hg19/v0/Homo_sapiens_assembly19.fasta ${REF37}/
    gsutil cp gs://gcp-public-data--broad-references/hg19/v0/Homo_sapiens_assembly19.fasta.fai ${REF37}/
    gsutil cp gs://gcp-public-data--broad-references/hg19/v0/Homo_sapiens_assembly19.dict ${REF37}/

    # BQSR files
    gsutil cp gs://gcp-public-data--broad-references/hg19/v0/dbsnp_138.b37.vcf.gz ${REF37}/
    gsutil cp gs://gcp-public-data--broad-references/hg19/v0/dbsnp_138.b37.vcf.gz.tbi ${REF37}/
    gsutil cp gs://gatk-legacy-bundles/b37/Mills_and_1000G_gold_standard.indels.b37.vcf ${REF37}/
    gsutil cp gs://gatk-legacy-bundles/b37/Mills_and_1000G_gold_standard.indels.b37.vcf.idx ${REF37}/
    gsutil cp gs://gatk-legacy-bundles/b37/1000G_phase1.indels.b37.vcf ${REF37}/
    gsutil cp gs://gatk-legacy-bundles/b37/1000G_phase1.indels.b37.vcf.idx ${REF37}/

    # Mutect2 files
    # NB: this PoN is WGS-derived. For panel data, consider building a matched PoN from your own normals for best filtering.
    gsutil cp gs://gatk-best-practices/somatic-b37/af-only-gnomad.raw.sites.vcf ${REF37}/
    gsutil cp gs://gatk-best-practices/somatic-b37/af-only-gnomad.raw.sites.vcf.idx ${REF37}/
    gsutil cp gs://gatk-best-practices/somatic-b37/Mutect2-WGS-panel-b37.vcf ${REF37}/
    gsutil cp gs://gatk-best-practices/somatic-b37/Mutect2-WGS-panel-b37.vcf.idx ${REF37}/
    gsutil cp gs://gatk-best-practices/somatic-b37/small_exac_common_3.vcf ${REF37}/
    gsutil cp gs://gatk-best-practices/somatic-b37/small_exac_common_3.vcf.idx ${REF37}/

    # Download GRCh37 release 87 GTF
    # Release 87 = last native Ensembl GRCh37 release. No chr prefix
    # Using Ensembl GTF rather than GENCODE because GENCODE uses chr prefix
    wget "https://ftp.ensembl.org/pub/grch37/release-87/gtf/homo_sapiens/Homo_sapiens.GRCh37.87.gtf.gz" \
         -O "${REF37}/Homo_sapiens.GRCh37.87.gtf.gz"
    gunzip "${REF37}/Homo_sapiens.GRCh37.87.gtf.gz"

    # make chromosome sizes file
    cut -f1,2 ${REF37}/Homo_sapiens_assembly19.fasta.fai > ${REF37}/GRCh37.genome
    echo "GRCh37 reference download complete."

    # Download files for vep cache
    mkdir -p ${REF37}/vep_cache
    wget https://ftp.ensembl.org/pub/release-113/variation/indexed_vep_cache/homo_sapiens_vep_113_GRCh37.tar.gz \
         -O ${REF37}/vep_cache/homo_sapiens_vep_113_GRCh37.tar.gz
    tar xzf ${REF37}/vep_cache/homo_sapiens_vep_113_GRCh37.tar.gz \
         -C ${REF37}/vep_cache

    # BWA index build
    echo "Detecting if BWA indexes are present"
    if [[ -d ${BWA_INDEX_DIR}/BWA_GRCh37 ]]
    then
        echo -en " * ${BWA_INDEX_DIR}/BWA_GRCh37 exists, no need to re-create\n\n"
    else
        echo -en " * ${BWA_INDEX_DIR}/BWA_GRCh37 not found, preparing BWA index...\n"
        mkdir ${BWA_INDEX_DIR}/BWA_GRCh37
        bwa index -p ${BWA_INDEX_DIR}/BWA_GRCh37/Homo_sapiens_assembly19 ${REF37}/Homo_sapiens_assembly19.fasta
    fi
fi
echo -en "*** All done ***"
