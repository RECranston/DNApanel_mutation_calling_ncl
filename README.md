# DNApanel_mutation_calling_ncl

Here are seven sequential bash scripts to run an end-to-end tumour-only analysis of targeted DNA panel sequencing mutation detection using [BWA-MEM](https://github.com/lh3/bwa), [gatk](https://github.com/broadinstitute/gatk) and [vep](https://github.com/ensembl/ensembl-vep) on the Newcastle University server (Comet) with slurm scheduler.

The pipeline is set up for the hybrid-capture panel defined by `NPHD2019A_Covered_paddel_fixed.sorted.bed` (downloaded by `script1_index_build.sh`), aligned to GRCh37. See [Important notes before running](#important-notes-before-running) for how to adapt it to a different panel.

The scripts include:
* `script1_index_build.sh`:
    * Downloads reference files for GRCh37 and GRCh38. Including files for use with gatk, BQSR, mutect2 and builds indexes/prepares reference files where needed
    * Creates vep caches
    * Builds BWA indexes
    * Downloads the panel bed file to `bed_files/`
* `script2_trimming.sh`:
    * Trims sequencing adapters from fastq.gz files using Trim Galore
    * Hard-code `FORMAT` parameter in the script to `FORMAT=SINGLE_ENDED` for single ended reads, or `FORMAT=PAIRED_END` for paired end reads
    * Runs FastQC
* `script3_align_bwa.sh`:
    * Performs BWA-MEM alignment on a sample list of trimmed fastq.gz files, then coordinate-sorts and indexes the bam.
    * Adds a read group to each bam (required by gatk). `LB` and `PU` are set to placeholder values (`lib1`, `unit1`) - edit these if your samples are split across multiple libraries or lanes.
    * Hard-code `FORMAT` parameter in the script to `FORMAT=SINGLE_ENDED` for single ended reads, or `FORMAT=PAIRED_END` for paired end reads
    * Hard-code `ASSEMBLY` parameter in the script to `ASSEMBLY=GRCh37` for alignment to GRCh37 genome build, or `ASSEMBLY=GRCh38` for alignment to GRCh38 genome build (must be `ASSEMBLY=GRCh37` for use with the supplied panel bed file)
* `script4_gatk_preprocessing.sh`:
    * Performs gatk preprocessing on named files (names derived from sample sheet)
    * Preprocessing includes: Mark duplicates (optional, see below), Base Quality Score Recalibration (BQSR) model building and application. BQSR is restricted to the panel bed intervals (+100 bp padding).
    * Hard-code `MARK_DUPLICATES` parameter in the script to `MARK_DUPLICATES=True` for hybrid-capture panels, or `MARK_DUPLICATES=False` for amplicon/PCR-based panels (where reads legitimately share start positions and must not be marked as duplicates)
    * Hard-code `ASSEMBLY` parameter in the script to `ASSEMBLY=GRCh37` for setting reference to GRCh37 genome build, or `ASSEMBLY=GRCh38` for GRCh38 genome build (must be `ASSEMBLY=GRCh37` for use with the supplied panel bed file)
* `script5_gatk_variant_calling.sh`:
    * Performs mutation detection on named files (names derived from sample sheet) using Mutect2 in tumour-only mode, with comparison to the reference genome, gnomAD germline resource and a panel of normals. Calling is restricted to the panel bed intervals (+100 bp padding).
    * Considers orientation bias (e.g. FFPE oxidative damage), sample contamination (estimated only at common population sites that fall within the panel) and filters variants accordingly, with a minimum allele fraction of 0.05.
    * Extracts filtered variants (see [Variant filtering](#variant-filtering-in-script5) - note this is not PASS-only)
    * Hard-code `ASSEMBLY` parameter in the script to `ASSEMBLY=GRCh37` for setting reference to GRCh37 genome build, or `ASSEMBLY=GRCh38` for GRCh38 genome build (must be `ASSEMBLY=GRCh37` for use with the supplied panel bed file)
* `script6_merge_phased_variants.sh`:
    * Merges phased Mutect2 records (sharing a phase set ID, `PID`) that sit directly next to each other into a single compound variant before annotation. vep annotates each record independently, so several small linked edits that together make one in-frame event can otherwise each be annotated as a spurious frameshift.
    * Merged records are tagged in the INFO field with `PHASED_MERGE_COUNT` (number of records merged).
    * Hard-code `MAX_PHASED_SPAN` parameter in the script to set the maximum bp span of records that will be merged (currently `10`)
    * Hard-code `ASSEMBLY` parameter in the script to `ASSEMBLY=GRCh37` for setting reference to GRCh37 genome build, or `ASSEMBLY=GRCh38` for GRCh38 genome build (must be `ASSEMBLY=GRCh37` for use with the supplied panel bed file)
* `script7_vep.sh`:
    * Annotation of merged, filtered variants by vep (`--everything`).
    * Outputs a bgzipped, tabix-indexed vcf and an uncompressed copy.
    * Hard-code `ASSEMBLY` parameter in the script to `ASSEMBLY=GRCh37` for setting reference to GRCh37 genome build, or `ASSEMBLY=GRCh38` for GRCh38 genome build (must be `ASSEMBLY=GRCh37` for use with the supplied panel bed file)
* `vcf_conversion.R`:
    * R script which converts the uncompressed vep vcf file output of `script7_vep.sh` into functional tables for downstream analysis.
    * Splits out the vep annotation (`CSQ`) into separate columns, with one row per transcript annotation per variant (so a variant overlapping several transcripts or genes appears on several rows).
    * Adds per-sample read evidence for each variant: total depth (`sample_DP`), reference and alternate read counts (`sample_AD_ref`, `sample_AD_alt`) and Mutect2's variant allele fraction (`sample_VAF`). These are prefixed `sample_` to distinguish them from vep's population allele frequency columns (e.g. `AF`, gnomAD frequencies).
    * Path to input directory (vep vcf file output) and path to output directory (to save output tables) to be defined in the script (`read.dir` and `output.dir`). Only files ending `.vcf` are read.
    * Requires R packages `vcfR` and `stringr` (`parallel` is included with base R). Files are processed in parallel using 20 cores (`mc.cores = 20`) - reduce this as needed.
    * Requires editing to include appropriate file paths for analysis.

### Important notes before running

* Genome build: the supplied panel bed file is in GRCh37 (b37) coordinates with no `chr` prefix, so `ASSEMBLY=GRCh37` must be used in scripts 3-7. GRCh38 references are still downloaded by `script1_index_build.sh` and the scripts will run on GRCh38, but only with a GRCh38 bed file.
* Using a different panel: place your own bed file in `bed_files/` and update the `PANEL_BED` path in both `script4_gatk_preprocessing.sh` and `script5_gatk_variant_calling.sh` (these must point to the same file). The bed file must be sorted and use the same contig names as the reference genome. Scripts 4 and 5 will exit with an error if the bed file is not found.
* Hybrid-capture vs amplicon panels: check `MARK_DUPLICATES` in `script4_gatk_preprocessing.sh` is set correctly for your panel type before running (see above).
* Panel of normals: the panel of normals downloaded by `script1_index_build.sh` is the gatk WGS-derived panel, not one built from normals sequenced on this panel. In tumour-only mode this limits how well panel-specific sequencing artefacts are removed. For best filtering, consider building a matched panel of normals from your own normal samples.
* Tumour-only calling: with no matched normal, some rare germline variants will remain in the output despite gnomAD and panel of normals filtering. Interpret results with this in mind.

### Setup

* Create a new directory and move into it.
* Git clone this repository.
```
git clone https://github.com/RECranston/DNApanel_mutation_calling_ncl.git
```
* Change into the cloned directory `cd DNApanel_mutation_calling_ncl`. Make all shell and R scripts executable
```
chmod 777 *.sh
chmod 777 *.R
```
* Please edit the script header of all scripts to assign the correct account name to the sbatch run.
* Run the setup script. References, indexes, the panel bed file and required files will be downloaded and built as required for all analysis stages.
```
sbatch ./script1_index_build.sh
```
* Move all fastq.gz files for analysis (or symlink using `ln -s`) to a sub-directory within the current directory.
* Create a tab-separated sample sheet of fastq.gz files and sample identifiers and save as a `.txt` file. If paired end, the sample sheet should be three columns including sample identifier followed by paired files, if single ended, the sample sheet should be two columns including sample identifier followed by the fastq.gz file. Examples are shown below:

Paired end:
```
sample_name  sample1_L001_R1.fastq.gz  sample1_L001_R2.fastq.gz
```
Single end:
```
sample_name  sample1.fastq.gz
```

### Run the pipeline
Parameters required for each script can be checked by running `./script_name.sh` in the terminal:
* `<tab delimited sample sheet>` is the name of the sample sheet (`.txt` file) of fastq.gz file locations and associated sample names created during the setup stage.
* `<input dir (relative)>` is the location of the directory containing input files relative to the current directory e.g. `trimmed_fastq/` (note the trailing "/").
* Similarly `<output dir (relative)>` is the location of the directory where the output data is to be stored, relative to the current directory e.g. `vep_output/` (note the trailing "/").
* Please edit the script header to assign the correct account name to the sbatch run and the correct number of jobs in the array, and jobs to be simultaneously performed.
  E.g. this example runs samples 1-10 from the sample sheet, running two samples at a time.
```
#SBATCH --array=1-10%2
```
This can usually be defined by the number of rows in the sample sheet e.g. `cat sample_sheet.txt | wc -l`
* Run scripts in sequential order. The input directory of each script is the output directory of the previous script.
* After `script2_trimming.sh` create a new tab delimited sample sheet including trimmed fastq.gz files and sample identifiers and save as a `.txt` file. This is the sample sheet used by `script3_align_bwa.sh`.
* If paired-end sequencing, the sample sheet of trimmed fastq.gz files should be three columns including sample identifier followed by paired files, if single ended, the sample sheet should be two columns including sample identifier followed by the trimmed fastq.gz file. Examples are shown below:

Paired end:
```
sample_name  sample1_L001_val_1.fastq.gz  sample1_L001_val_2.fastq.gz
```
Single end:
```
sample_name  sample1_trimmed.fastq.gz
```
* Scripts 4-7 only use the first column (sample identifier) of the sample sheet, so either the original or the trimmed sample sheet can be used, as long as the row order is the same.

Example run order:
```
sbatch ./script2_trimming.sh sample_sheet.txt fastq/ trimmed_fastq/
sbatch ./script3_align_bwa.sh trimmed_sample_sheet.txt trimmed_fastq/ aligned/
sbatch ./script4_gatk_preprocessing.sh trimmed_sample_sheet.txt aligned/ output_preprocessing/
sbatch ./script5_gatk_variant_calling.sh trimmed_sample_sheet.txt output_preprocessing/ output_mutation_calling/
sbatch ./script6_merge_phased_variants.sh trimmed_sample_sheet.txt output_mutation_calling/ output_merged_phased/
sbatch ./script7_vep.sh trimmed_sample_sheet.txt output_merged_phased/ output_vep/
```
### Variant filtering in script5

After `FilterMutectCalls`, variants are extracted to `<sample>_tumor_filtered_PASS.vcf.gz`. Despite the file name, this file contains three groups of variants:
1. Variants with `FILTER` = `PASS`.
2. Variants whose only filter reason is `orientation` (possible orientation-bias artefact). These are kept as soft-flagged rather than removed, as real variants can be flagged here at panel depth. The `orientation` label stays in the `FILTER` column so they can be identified and excluded downstream if needed.
3. Variants whose filter reasons are only `clustered_events;haplotype` and that have a tumour log odds score (`TLOD`) above 100. These are typically real compound variants (e.g. several nearby linked changes) that Mutect2 flags because they cluster together. A TLOD above 100 means the evidence for a real variant over sequencing noise is very strong.

All other filtered variants are removed. The thresholds can be changed via `ORIENTATION_FILTER_TAG`, `CLUSTERED_HAPLOTYPE_TAG` and `CLUSTERED_HAPLOTYPE_TLOD_THRESHOLD` in the script.

### Output files (per sample)  

These are the expected output files per sample for each analysis script:  
* `script3_align_bwa.sh`  
    * `<sample>_sorted.bam` (+ `.bai`)  
* `script4_gatk_preprocessing.sh`  
    * `<sample>_recal.bam` (recalibrated, panel regions only)
    * `<sample>_marked_dup_metrics.txt`  
* `script5_gatk_variant_calling.sh`  
    * `<sample>_tumor_raw.vcf.gz`  
    * `<sample>_tumor_filtered.vcf.gz` (all calls with filter labels)  
    * `<sample>_tumor_filtered_PASS.vcf.gz` 
    * `<sample>_contamination.table`
* `script6_merge_phased_variants.sh`  
    * `<sample>_tumor_filtered_PASS_merged.vcf.gz` (+ `.tbi`)  
* `script7_vep.sh`  
    * `<sample>_tumor_annotated.vcf.gz` (+ `.tbi`)  
    * Uncompressed `<sample>_tumor_annotated.vcf`  
* `vcf_conversion.R`  
    * `<sample>_final.txt` (tab-delimited table) 
* Resulting filtered, merged, vep-annotated mutation data is saved to the defined output directory of `script7_vep.sh` as `.vcf` files.  
* Utilise the `vcf_conversion.R` script (outlined above) to process vep vcf output files into readable tables for downstream analysis. Set `read.dir` to the `script7_vep.sh` output directory and set `output.dir` to an output directory to save files to, then run e.g.  
```
Rscript vcf_conversion.R
```
* The `FILTER` column is carried through to the final tables, so soft-flagged `orientation` and rescued `clustered_events;haplotype` variants (see [Variant filtering](#variant-filtering-in-script5)) can be identified and excluded at this stage if required.
* Logs for each array job are written to `logs/`. Merged phased groups are reported in the `script6_merge_phased_variants.sh` log.

