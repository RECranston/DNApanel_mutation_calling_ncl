#!/bin/bash
#SBATCH --account=XXXXX
#SBATCH --partition=default_free
#SBATCH --mem=60G
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=8
#SBATCH --job-name=align
#SBATCH --output=logs/align_%A_%a.out
#SBATCH --array=1-100%100

# Script to run an array of BWA-MEM alignment jobs on a paired list of trimmed fastq files (DNA panel sequencing)
# Ruth Cranston 2026

[ $# -ne 3 ] && { echo -en \
"\nRuth Cranston 2026\n\n
*** Script to run BWA-MEM alignment for a paired list of trimmed fastq files [sample name] [fastq1] [fastq2] (\
tab delimited sheet).
Runs in current directory. Output directory is created.
<sample sheet> <input dir (relative)> <output dir (relative)>
example run: sbatch ./script3_align_bwa.sh trimmed_sample_sheet.txt test_trimmed_fastq/ test_aligned_array/ *** \
\n\n" ; exit 1; }

# example run
# sbatch --array=1-5 script3_align_bwa.sh test_trimmed_sample_sheet.txt test_trimmed_fastq/ test_aligned_array/
# --array=1-5%10 means run array job IDs 1-5 with a maximum of 10 running at once

# Set variables
BASE_DIR="$PWD"
ASSEMBLY="GRCh37"
BWA_INDEX_DIR=${BASE_DIR}/BWA_indexes/BWA_${ASSEMBLY}
if [[ "${ASSEMBLY}" == "GRCh38" ]]; then
    INDEX_PREFIX=${BWA_INDEX_DIR}/Homo_sapiens_assembly38
elif [[ "${ASSEMBLY}" == "GRCh37" ]]; then
    INDEX_PREFIX=${BWA_INDEX_DIR}/Homo_sapiens_assembly19
else
    echo "ASSEMBLY must be GRCh38 or GRCh37" >&2
    exit 1
fi
SAMPLE_SHEET=$1
INPUT_DIR=${BASE_DIR}/$2
OUTPUT_DIR=${BASE_DIR}/$3

# Set this to PAIRED_END or SINGLE_ENDED
FORMAT=PAIRED_END

# Load modules
echo -en " * Loading modules...\n"

set -euo pipefail

module --force purge
module load BWA/0.7.18-GCCcore-13.3.0
module load SAMtools
echo -en " * Environment set up.\n"

mkdir -p ${OUTPUT_DIR}
mkdir -p logs

# Get the correct row for this array task
LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${SAMPLE_SHEET})

SAMPLE_ID=$(echo $LINE | awk '{print $1}')
FILE1=$(echo $LINE | awk '{print $2}')
if [ "${FORMAT}" == "PAIRED_END" ]; then
    FILE2=$(echo $LINE | awk '{print $3}')
fi

echo "Processing sample: ${SAMPLE_ID}"
echo "Task ID: ${SLURM_ARRAY_TASK_ID}"
echo "Format: ${FORMAT}"

if [ "${FORMAT}" == "PAIRED_END" ]; then
    READS_IN="${INPUT_DIR}${FILE1} ${INPUT_DIR}${FILE2}"
else
    READS_IN="${INPUT_DIR}${FILE1}"
fi

# Read group line - required by GATK downstream. LB/PU are placeholders as in the original script;
# update if your samples are multiplexed across libraries/lanes.
RG="@RG\tID:${SAMPLE_ID}\tSM:${SAMPLE_ID}\tPL:ILLUMINA\tLB:lib1\tPU:unit1"

# Run BWA-MEM alignment, piped straight into coordinate sort
# -M marks shorter split hits as secondary, for Picard/GATK compatibility downstream
bwa mem \
    -M \
    -t ${SLURM_CPUS_PER_TASK} \
    -R "${RG}" \
    ${INDEX_PREFIX} \
    ${READS_IN} \
    | samtools sort -@ ${SLURM_CPUS_PER_TASK} -m 1G -o ${OUTPUT_DIR}${SAMPLE_ID}_sorted.bam -
echo -ne "*** BWA-MEM alignment done! ***\n"

# Index sorted bam
samtools index ${OUTPUT_DIR}${SAMPLE_ID}_sorted.bam
echo -ne "*** Bam indexing done! ***\n"
