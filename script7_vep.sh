#!/bin/bash
#SBATCH --account=XXXXX
#SBATCH --partition=default_free
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=8
#SBATCH --job-name=vep_annotate
#SBATCH --output=logs/vep_%A_%a.out
#SBATCH --array=1-100%100

# Script to run an array of VEP annotation jobs on merged, PASS-filtered DNA panel variants
# Ruth Cranston 2026

[ $# -ne 3 ] && { echo -en \
"\nRuth Cranston 2026\n\n
*** Script to run VEP annotation jobs on a list of sample ids from the original sample sheet \
[sample name] [fastq1] [fastq2] (tab delimited sheet).
Runs in current directory. Input dir is location of merged phased-variant VCFs (script6 output). \
Output directory is created.
<sample sheet> <input dir (relative)> <output dir (relative)>
example run: sbatch ./script7_vep.sh sample_sheet.txt output_merged_phased/ output_vep/ *** \n\n" ; exit 1; }

# --array=1-5%10 means run array job IDs 1-5 with a maximum of 10 running at once

# Set variables
BASE_DIR="$PWD"
ASSEMBLY="GRCh37"
SAMPLE_SHEET=$1
INPUT_DIR=${BASE_DIR}/$2
OUTPUT_DIR=${BASE_DIR}/$3
REFERENCE_DIR=${BASE_DIR}/References/${ASSEMBLY}

# Load modules

echo -en " * Loading modules...\n"
module --force purge
module load VEP/113.3-GCC-13.3.0
module load HTSlib/1.21-GCC-13.3.0

set -euo pipefail

echo -en " * Environment set up.\n"

# make output dir
mkdir -p ${OUTPUT_DIR}
mkdir -p logs

# Get the correct row for this array task
LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${SAMPLE_SHEET})
SAMPLE_ID=$(echo $LINE | awk '{print $1}')

if [[ -z "${SAMPLE_ID}" ]]; then
    echo "ERROR: SAMPLE_ID is empty for array task ${SLURM_ARRAY_TASK_ID} - check row ${SLURM_ARRAY_TASK_ID} of ${SAMPLE_SHEET} exists and is correctly formatted" >&2
    exit 1
fi

echo "Processing sample: ${SAMPLE_ID}"
echo "Task ID: ${SLURM_ARRAY_TASK_ID}"

# Set reference fasta
if [ "${ASSEMBLY}" == "GRCh38" ]; then
    REF_FASTA=${REFERENCE_DIR}/Homo_sapiens_assembly38.fasta
else
    REF_FASTA=${REFERENCE_DIR}/Homo_sapiens_assembly19.fasta
fi
INPUT_VCF=${INPUT_DIR}${SAMPLE_ID}_tumor_filtered_PASS_merged.vcf.gz
OUTPUT_VCF=${OUTPUT_DIR}${SAMPLE_ID}_tumor_annotated.vcf.gz
if [[ ! -f "${INPUT_VCF}" ]]; then
    echo "ERROR: expected input VCF not found at ${INPUT_VCF} (check script6_merge_phased_variants.sh has run for this sample)" >&2
    exit 1
fi

# VEP annotation of merged, PASS-filtered Mutect2 variants (VEP reads .vcf.gz natively, no need to decompress)
vep \
    --input_file ${INPUT_VCF} \
    --output_file ${OUTPUT_VCF} \
    --compress_output bgzip \
    --format vcf --vcf \
    --cache --offline \
    --dir_cache ${REFERENCE_DIR}/vep_cache \
    --assembly ${ASSEMBLY} \
    --everything \
    --fork ${SLURM_CPUS_PER_TASK} \
    --fasta ${REF_FASTA}

echo -ne "*** VEP annotation complete! ***\n"

tabix -p vcf ${OUTPUT_VCF}
echo -ne "*** Indexing complete! ***\n"

bgzip -d -k ${OUTPUT_VCF}
echo -ne "*** Uncompressed copy written! ***\n"

echo -ne "*** All done! ***\n"
