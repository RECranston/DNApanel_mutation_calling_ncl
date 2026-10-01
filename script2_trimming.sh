#!/bin/bash
#SBATCH --account=XXXXX
#SBATCH --partition=default_free
#SBATCH --mem=50G
#SBATCH --time=23:00:00
#SBATCH --cpus-per-task=8
#SBATCH --job-name=trimming
#SBATCH --output=logs/trimming_%A_%a.out
#SBATCH --array=1-100%100

# Script to run Trim_Galore on a paired list of fastq files
# Ruth Cranston 2026

# example run
# sbatch script2_trimming.sh sample_sheet.txt input_dir/ output_dir/
# --array=1-100%20 means run array job IDs 1-100 with a maximum of 20 running at once

[ $# -ne 3 ] && { echo -en \
"\nRuth Cranston 2026\n\n
*** Script to run trim_galore for a list of fastq files - either paired end [sample name] [fastq1] [fastq2] or single ended [sample name] [fastq1] (tab delimited sheet). 
Define single/paired end status by setting the FORMAT parameter within the script. Runs as an array job in current directory. Output directory is created. 
<sample sheet> <input dir (relative)> <output dir (relative)>*** \n\n" ; exit 1; }

# Set variables
BASE_DIR="$PWD"
THREADS=8
SAMPLE_SHEET=$1
INPUT_DIR=${BASE_DIR}/$2
OUTPUT_DIR=${BASE_DIR}/$3

# Set if paired end or single ended (PAIRED_END/SINGLE_ENDED)
FORMAT=PAIRED_END

# Load modules
echo -en " * Loading modules...\n"
module --force purge
module load FastQC
module load MultiQC
module load pigz
module load cutadapt
module load Trim_Galore
echo -en " * Environment set up.\n"

# Remove output dir file if exists
echo -en " * Creating output directory...\n"
mkdir -p ${OUTPUT_DIR}
mkdir -p ${OUTPUT_DIR}fastqc_output

# Get the correct row for this array task
LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${SAMPLE_SHEET})

# Get files from sheet
SAMPLE_ID=$(echo $LINE | awk '{print $1}')
FILE1=$(echo $LINE | awk '{print $2}')
if [ "${FORMAT}" == "PAIRED_END" ]; then
    FILE2=$(echo $LINE | awk '{print $3}')
fi

echo "Processing sample: ${SAMPLE_ID}"
echo "Task ID: ${SLURM_ARRAY_TASK_ID}"
echo "Format: ${FORMAT}"

cd ${OUTPUT_DIR}

if [ "${FORMAT}" == "PAIRED_END" ]; then

    trim_galore --cores ${THREADS} \
		--fastqc \
		--paired \
		--dont_gzip \
		--output_dir ${OUTPUT_DIR} \
		${INPUT_DIR}${FILE1} ${INPUT_DIR}${FILE2}

    # Derive base names
    FILE1_BASE=$(basename ${FILE1} .fastq.gz)
    FILE2_BASE=$(basename ${FILE2} .fastq.gz)

    # Compress outputs
    echo "Compressing ${SAMPLE_ID} trimmed fastq files"
    pigz ${OUTPUT_DIR}${FILE1_BASE}_val_1.fq
    pigz ${OUTPUT_DIR}${FILE2_BASE}_val_2.fq

    # Rename & clean up
    echo "Renaming and cleaning up ${SAMPLE_ID} output files"
    mv ${OUTPUT_DIR}${FILE1_BASE}_val_1.fq.gz ${OUTPUT_DIR}${FILE1_BASE}_val_1.fastq.gz
    mv ${OUTPUT_DIR}${FILE2_BASE}_val_2.fq.gz ${OUTPUT_DIR}${FILE2_BASE}_val_2.fastq.gz

    mv ${OUTPUT_DIR}${FILE1_BASE}.fastq.gz_trimming_report.txt fastqc_output/

    mv ${OUTPUT_DIR}${FILE1_BASE}_val_1_fastqc.html fastqc_output/
    mv ${OUTPUT_DIR}${FILE1_BASE}_val_1_fastqc.zip fastqc_output/

    mv ${OUTPUT_DIR}${FILE2_BASE}.fastq.gz_trimming_report.txt fastqc_output/
    mv ${OUTPUT_DIR}${FILE2_BASE}_val_2_fastqc.html fastqc_output/
    mv ${OUTPUT_DIR}${FILE2_BASE}_val_2_fastqc.zip fastqc_output/

else

    trim_galore --cores ${THREADS} \
		--fastqc \
		--dont_gzip \
		--output_dir ${OUTPUT_DIR} \
		${INPUT_DIR}${FILE1}
 
    # Derive base name
    FILE1_BASE=$(basename ${FILE1} .fastq.gz)
 
    # Compress outputs
    echo "Compressing ${SAMPLE_ID} trimmed fastq file"
    pigz ${OUTPUT_DIR}${FILE1_BASE}_trimmed.fq
 
    # Rename & clean up
    echo "Renaming and cleaning up ${SAMPLE_ID} output files"
    mv ${OUTPUT_DIR}${FILE1_BASE}_trimmed.fq.gz ${OUTPUT_DIR}${FILE1_BASE}_trimmed.fastq.gz
    mv ${OUTPUT_DIR}${FILE1_BASE}.fastq.gz_trimming_report.txt fastqc_output/
    mv ${OUTPUT_DIR}${FILE1_BASE}_trimmed_fastqc.html fastqc_output/
    mv ${OUTPUT_DIR}${FILE1_BASE}_trimmed_fastqc.zip fastqc_output/
 
fi

echo "${SAMPLE_ID} processing complete"

