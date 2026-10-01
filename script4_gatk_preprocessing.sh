#!/bin/bash
#SBATCH --account=XXXXX
#SBATCH --partition=default_free
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=8
#SBATCH --job-name=gatk_preprocess
#SBATCH --output=logs/gatk_preprocess_%A_%a.out
#SBATCH --array=1-100%100

# Script to run an array of GATK preprocessing jobs (MarkDuplicates, BQSR) on aligned DNA panel BAMs
# Ruth Cranston 2026

[ $# -ne 3 ] && { echo -en \
"\nRuth Cranston 2026\n\n
*** Script to run GATK preprocessing (MarkDuplicates + BQSR) on a list of sample ids from the original sample \
sheet [sample name] [fastq1] [fastq2] (tab delimited sheet).
Runs in current directory. Output directory is created.
<sample sheet> <input dir (relative, BWA output)> <output dir (relative)>
example run: sbatch ./script4_gatk_preprocessing.sh sample_sheet.txt test_aligned_array/ output_preprocessing/ *** \n\n" ; exit 1; }

# Set variables - must be GRCh37 for bed file to work
BASE_DIR="$PWD"
ASSEMBLY="GRCh37"
REFERENCE_DIR=${BASE_DIR}/References/${ASSEMBLY}
PANEL_BED="${BASE_DIR}/bed_files/NPHD2019A_Covered_paddel_fixed.sorted.bed"
INTERVAL_PADDING=100
TMPDIR=${BASE_DIR}/tmp
SAMPLE_SHEET=$1
INPUT_DIR=${BASE_DIR}/$2
OUTPUT_DIR=${BASE_DIR}/$3

# Set to False for amplicon/PCR-based panels, where shared read start positions are expected
# and should NOT be marked as duplicates. Leave True for hybrid-capture panels.
MARK_DUPLICATES=True

# Load modules
echo -en " * Loading modules...\n"
module --force purge
module load GATK/4.6.0.0-GCCcore-13.2.0-Java-17
module load SAMtools

set -euo pipefail

echo -en " * Environment set up.\n"

# make output/tmp dirs
mkdir -p ${OUTPUT_DIR}
mkdir -p logs
mkdir -p ${TMPDIR}

if [[ ! -f "${PANEL_BED}" ]]; then
    echo "Panel BED not found at ${PANEL_BED}" >&2
    exit 1
fi

# Get the correct row for this array task
LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${SAMPLE_SHEET})
SAMPLE_ID=$(echo $LINE | awk '{print $1}')

echo "Processing sample: ${SAMPLE_ID}"
echo "Task ID: ${SLURM_ARRAY_TASK_ID}"

# Set genome reference files
if [ "${ASSEMBLY}" == "GRCh38" ]; then

    REF_FASTA=${REFERENCE_DIR}/Homo_sapiens_assembly38.fasta
    DBSNP=${REFERENCE_DIR}/Homo_sapiens_assembly38.dbsnp138.vcf
    KNOWN_INDELS_1=${REFERENCE_DIR}/Homo_sapiens_assembly38.known_indels.vcf.gz
    KNOWN_INDELS_2=${REFERENCE_DIR}/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz

else

    REF_FASTA=${REFERENCE_DIR}/Homo_sapiens_assembly19.fasta
    DBSNP=${REFERENCE_DIR}/dbsnp_138.b37.vcf.gz
    KNOWN_INDELS_1=${REFERENCE_DIR}/Mills_and_1000G_gold_standard.indels.b37.vcf
    KNOWN_INDELS_2=${REFERENCE_DIR}/1000G_phase1.indels.b37.vcf

fi

SORTED_BAM=${INPUT_DIR}${SAMPLE_ID}_sorted.bam
MARKED_DUP_BAM=${OUTPUT_DIR}${SAMPLE_ID}_marked_dup.bam

if [ "${MARK_DUPLICATES}" == "True" ]; then

    # Mark duplicates
    gatk MarkDuplicates \
         --java-options "-Xmx28g -Djava.io.tmpdir=${TMPDIR}" \
         -I ${SORTED_BAM} \
         -O ${MARKED_DUP_BAM} \
         -M ${OUTPUT_DIR}${SAMPLE_ID}_marked_dup_metrics.txt \
         --CREATE_INDEX true \
         --TMP_DIR ${TMPDIR}
    echo -ne "*** Mark duplicates done! ***\n"

else

    echo -ne " * MARK_DUPLICATES is false (amplicon panel) - carrying sorted BAM through unmarked\n"
    cp ${SORTED_BAM} ${MARKED_DUP_BAM}
    samtools index ${MARKED_DUP_BAM}

fi

# catch here if MarkDuplicates (or the copy step above) failed to produce usable output
if [ ! -s ${MARKED_DUP_BAM} ]; then
    echo "ERROR: no usable BAM at ${MARKED_DUP_BAM} for ${SAMPLE_ID}" >&2
    exit 1
fi

# BQSR (dbSNP + Mills + second known-indels resource, restricted to panel intervals)
gatk BaseRecalibrator \
     --java-options "-Xmx28g -Djava.io.tmpdir=${TMPDIR}" \
     -R ${REF_FASTA} \
     -I ${MARKED_DUP_BAM} \
     -L ${PANEL_BED} \
     --interval-padding ${INTERVAL_PADDING} \
     --known-sites ${DBSNP} \
     --known-sites ${KNOWN_INDELS_1} \
     --known-sites ${KNOWN_INDELS_2} \
     -O ${OUTPUT_DIR}${SAMPLE_ID}_recal.table \
     --tmp-dir ${TMPDIR}

echo -ne "*** BaseRecalibrator done! ***\n"

gatk ApplyBQSR \
     --java-options "-Xmx28g -Djava.io.tmpdir=${TMPDIR}" \
     -R ${REF_FASTA} \
     -I ${MARKED_DUP_BAM} \
     -L ${PANEL_BED} \
     --interval-padding ${INTERVAL_PADDING} \
     --bqsr-recal-file ${OUTPUT_DIR}${SAMPLE_ID}_recal.table \
     -O ${OUTPUT_DIR}${SAMPLE_ID}_recal.bam \
     --tmp-dir ${TMPDIR}

echo -ne "*** BQSR done! ***\n"

echo -ne "*** All done! ***\n"
