#!/bin/bash
#SBATCH --account=XXXXX
#SBATCH --partition=default_free
#SBATCH --mem=32G
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=8
#SBATCH --job-name=gatk_var_calling
#SBATCH --output=logs/gatk_var_calling_%A_%a.out
#SBATCH --array=1-100%100

# Script to run an array of tumour-only Mutect2 DNA panel variant calling jobs
# Ruth Cranston 2026

[ $# -ne 3 ] && { echo -en \
"\nRuth Cranston 2026\n\n
*** Script to run tumour-only GATK Mutect2 mutation calling on a list of sample ids from the original sample \
sheet [sample name] [fastq1] [fastq2] (tab delimited sheet).
Runs in current directory. Input dir is location of GATK preprocessed (recalibrated) BAMs. Output directory is created.
<sample sheet> <input dir (relative)> <output dir (relative)>
example run: sbatch ./script5_gatk_variant_calling.sh sample_sheet.txt output_preprocessing/ output_mutation_calling/ *** \n\n" ; exit 1; }

# --array=1-5%10 means run array job IDs 1-5 with a maximum of 10 running at once

# Set variables
BASE_DIR="$PWD"
ASSEMBLY="GRCh37"
REFERENCE_DIR=${BASE_DIR}/References/${ASSEMBLY}

# Should match the PANEL_BED used in script4 - keep every stage pointed at the same file
PANEL_BED="${BASE_DIR}/bed_files/NPHD2019A_Covered_paddel_fixed.sorted.bed"
INTERVAL_PADDING=100
SAMPLE_SHEET=$1
INPUT_DIR=${BASE_DIR}/$2
OUTPUT_DIR=${BASE_DIR}/$3

# Load modules
echo -en " * Loading modules...\n"
module --force purge
module load GATK/4.6.0.0-GCCcore-13.2.0-Java-17
module load SAMtools/1.21-GCC-13.3.0
module load HTSlib/1.21-GCC-13.3.0

set -euo pipefail

echo -en " * Environment set up.\n"

# make output dir
mkdir -p ${OUTPUT_DIR}
mkdir -p logs

if [[ ! -f "${PANEL_BED}" ]]; then
    echo "Panel BED not found at ${PANEL_BED}" >&2
    exit 1
fi

# Get the correct row for this array task
LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${SAMPLE_SHEET})
SAMPLE_ID=$(echo $LINE | awk '{print $1}')

if [[ -z "${SAMPLE_ID}" ]]; then
    echo "ERROR: SAMPLE_ID is empty for array task ${SLURM_ARRAY_TASK_ID} - check row ${SLURM_ARRAY_TASK_ID} of ${SAMPLE_SHEET} exists and is correctly formatted" >&2
    exit 1
fi

echo "Processing sample: ${SAMPLE_ID}"
echo "Task ID: ${SLURM_ARRAY_TASK_ID}"

# Set reference files
if [ "${ASSEMBLY}" == "GRCh38" ]; then
    REF_FASTA=${REFERENCE_DIR}/Homo_sapiens_assembly38.fasta
    GNOMAD=${REFERENCE_DIR}/af-only-gnomad.hg38.vcf.gz
    PON=${REFERENCE_DIR}/1000g_pon.hg38.vcf.gz
    EXAC=${REFERENCE_DIR}/small_exac_common_3.hg38.vcf.gz
else
    REF_FASTA=${REFERENCE_DIR}/Homo_sapiens_assembly19.fasta
    GNOMAD=${REFERENCE_DIR}/af-only-gnomad.raw.sites.vcf
    PON=${REFERENCE_DIR}/Mutect2-WGS-panel-b37.vcf
    EXAC=${REFERENCE_DIR}/small_exac_common_3.vcf
fi

RECAL_BAM=${INPUT_DIR}${SAMPLE_ID}_recal.bam

# Mutect2 (tumor-only mode), restricted to panel intervals.

gatk Mutect2 \
     --java-options "-Xmx28g" \
     -R ${REF_FASTA} \
     -I ${RECAL_BAM} \
     --tumor-sample ${SAMPLE_ID} \
     -L ${PANEL_BED} \
     --interval-padding ${INTERVAL_PADDING} \
     --germline-resource ${GNOMAD} \
     --panel-of-normals ${PON} \
     --native-pair-hmm-threads ${SLURM_CPUS_PER_TASK} \
     --f1r2-tar-gz ${OUTPUT_DIR}${SAMPLE_ID}_f1r2.tar.gz \
     -O ${OUTPUT_DIR}${SAMPLE_ID}_tumor_raw.vcf.gz

echo -ne "*** Mutect2 finished! ***\n"

# Learn orientation bias model. This targets a real DNA artefact (FFPE oxidative damage in  particular) (--obs-priors option explicitly
# excluded)

gatk LearnReadOrientationModel \
     --java-options "-Xmx28g" \
     -I ${OUTPUT_DIR}${SAMPLE_ID}_f1r2.tar.gz \
     -O ${OUTPUT_DIR}${SAMPLE_ID}_artifact_prior.tar.gz
echo -ne "*** LearnReadOrientationModel finished! ***\n"

# Get contamination estimate, restricted to sites that are both common-population sites AND  within the panel (small_exac_common_3
# spans the whole exome, most of which we have no coverage for)

gatk GetPileupSummaries \
    --java-options "-Xmx28g" \
    -I ${RECAL_BAM} \
    -V ${EXAC} \
    -L ${EXAC} \
    -L ${PANEL_BED} \
    --interval-padding ${INTERVAL_PADDING} \
    --interval-set-rule INTERSECTION \
    -O ${OUTPUT_DIR}${SAMPLE_ID}_pileup_summaries.table

gatk CalculateContamination \
    --java-options "-Xmx28g" \
    -I ${OUTPUT_DIR}${SAMPLE_ID}_pileup_summaries.table \
    -O ${OUTPUT_DIR}${SAMPLE_ID}_contamination.table

echo -ne "*** CalculateContamination finished! ***\n"

# Apply filters, including orientation bias priors. Not hard-excluded here

gatk FilterMutectCalls \
     --java-options "-Xmx28g" \
     -R ${REF_FASTA} \
     -V ${OUTPUT_DIR}${SAMPLE_ID}_tumor_raw.vcf.gz \
     --contamination-table ${OUTPUT_DIR}${SAMPLE_ID}_contamination.table \
     --ob-priors ${OUTPUT_DIR}${SAMPLE_ID}_artifact_prior.tar.gz \
     --min-allele-fraction 0.05 \
     -O ${OUTPUT_DIR}${SAMPLE_ID}_tumor_filtered.vcf.gz

echo -ne "*** Mutect call filtering finished! ***\n"

# Extract variants: keep PASS, keep variants whose ONLY filter reason is orientation bias (soft-flagged) and keep clustered_events -
# haplotype calls above a TLOD (tumour log odds) confidence threshold (>100 - v confident). (how likely is it a tumour variant vs noise variant)

ORIENTATION_FILTER_TAG="orientation"
CLUSTERED_HAPLOTYPE_TAG="clustered_events;haplotype"
CLUSTERED_HAPLOTYPE_TLOD_THRESHOLD=100

EXTRACT_VARIANTS='
    function get_tlod(info,    n, i, arr, kv) {
        n = split(info, arr, ";")
        for (i = 1; i <= n; i++) {
            if (arr[i] ~ /^TLOD=/) {
                split(arr[i], kv, "=")
                return kv[2] + 0
            }
        }
        return -1
    }
    /^#/ { print; next }
    $7 == "PASS" { print; next }
    $7 == orient_tag { print; next }
    $7 == cluster_tag {
        if (get_tlod($8) > tlod_min) print
        next
    }
'

if command -v bgzip >/dev/null 2>&1 && command -v tabix >/dev/null 2>&1; then

    zcat ${OUTPUT_DIR}${SAMPLE_ID}_tumor_filtered.vcf.gz | \
        awk -F'\t' -v OFS='\t' \
            -v orient_tag="${ORIENTATION_FILTER_TAG}" \
            -v cluster_tag="${CLUSTERED_HAPLOTYPE_TAG}" \
            -v tlod_min="${CLUSTERED_HAPLOTYPE_TLOD_THRESHOLD}" \
            "${EXTRACT_VARIANTS}" | \
        bgzip > ${OUTPUT_DIR}${SAMPLE_ID}_tumor_filtered_PASS.vcf.gz


    tabix -p vcf ${OUTPUT_DIR}${SAMPLE_ID}_tumor_filtered_PASS.vcf.gz
else

    echo "WARNING: bgzip/tabix not found in PATH - writing uncompressed VCF instead (fix the module load once you know where bgzip lives on this cluster)" >&2
    zcat ${OUTPUT_DIR}${SAMPLE_ID}_tumor_filtered.vcf.gz | \
        awk -F'\t' -v OFS='\t' \
            -v orient_tag="${ORIENTATION_FILTER_TAG}" \
            -v cluster_tag="${CLUSTERED_HAPLOTYPE_TAG}" \
            -v tlod_min="${CLUSTERED_HAPLOTYPE_TLOD_THRESHOLD}" \
            "${EXTRACT_VARIANTS}" \
        > ${OUTPUT_DIR}${SAMPLE_ID}_tumor_filtered_PASS.vcf

fi

echo -ne "*** Extract PASS + rescued variants finished! ***\n"

echo -ne "*** All done! ***\n"
