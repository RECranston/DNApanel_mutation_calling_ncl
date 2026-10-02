# BiocManager::install("vcfR", lib = "~/R/x86_64-pc-linux-gnu-library/4.5/")
# BiocManager::install("stringr", lib = "~/R/x86_64-pc-linux-gnu-library/4.5/")
# BiocManager::install("parallel", lib = "~/R/x86_64-pc-linux-gnu-library/4.5/")

library(vcfR)
library(stringr)
library(parallel)

# read in a demo vcf file and get into a workable format
read.dir <- "~/path/to/input/"
file.path.list <- paste0(read.dir, dir(read.dir))
output.dir <- "~/path/to/output/"

# get vcf files only
file.path.list <- file.path.list[grep("vcf$", file.path.list)]

mclapply(file.path.list, mc.cores = 20, function(x){
  
  # to test:
  # x = file.path.list[1]
  
  # Set file path
  temp.file.path <- x
  
  # Read in vcf
  temp.vcf <- read.vcfR(temp.file.path)
  
  # Get file id
  temp.id <- gsub(read.dir, "", x)
  temp.id <- gsub("_tumor_annotated.vcf", "", temp.id)
  
  # Fixed fields: CHROM, POS, ID, REF, ALT, QUAL, FILTER
  fix_df <- as.data.frame(getFIX(temp.vcf), stringsAsFactors = FALSE)
  fix_df$POS <- as.numeric(fix_df$POS)
  fix_df$QUAL <- as.numeric(fix_df$QUAL)
  
  # Genotype matrix (one column per sample; tumor-only here, so one column)
  gt <- extract.gt(temp.vcf, element = "GT")                     # genotypes e.g. 0/1
  dp <- extract.gt(temp.vcf, element = "DP", as.numeric = TRUE)  # total depth
  ad <- extract.gt(temp.vcf, element = "AD", as.numeric = FALSE) # "ref,alt" — can't be numeric, comma-separated
  gq <- extract.gt(temp.vcf, element = "GQ", as.numeric = TRUE)  # genotype quality
  af <- extract.gt(temp.vcf, element = "AF", as.numeric = TRUE)  # Mutect2's own computed allele fraction
  
  # Split AD into ref/alt depth counts
  ad_split <- str_split_fixed(ad[,1], ",", 2)
  
  # Attach to fix_df (row order matches getFIX() 1:1, so this is a direct positional bind)
  # Namespaced as "sample_*" to avoid colliding with the population-level "AF" column
  # that comes from VEP's CSQ (1000G/gnomAD frequencies) later in the script
  fix_df$sample_DP     <- dp[,1]
  fix_df$sample_AD_ref <- as.numeric(ad_split[,1])
  fix_df$sample_AD_alt <- as.numeric(ad_split[,2])
  fix_df$sample_VAF    <- af[,1]
  
  # Full genotype block as tidy long-format data frame
  gt_tidy <- vcfR2tidy(temp.vcf, single_frame = FALSE)
  # gt_tidy$fix   -> fixed fields + INFO already split into columns
  # gt_tidy$gt    -> long-format genotype data (one row per sample-variant)
  # gt_tidy$meta  -> header metadata (field descriptions/types)
  
  # VEP CSQ annotation: extract field names from header
  csq_format <- temp.vcf@meta[grep("ID=CSQ", temp.vcf@meta)]
  csq_fields <- strsplit(sub('.*Format: ([^"]+)".*', "\\1", csq_format), "\\|")[[1]]
  
  csq_strings <- extract.info(temp.vcf, "CSQ")
  
  # Split each variant's CSQ into one string per transcript (comma-separated)
  csq_list <- str_split(csq_strings, ",")
  
  # Track which original variant row each transcript annotation belongs to
  n_transcripts <- lengths(csq_list)
  variant_idx <- rep(seq_along(csq_list), n_transcripts)
  
  # Flatten to one transcript-annotation string per element, then split on |
  csq_flat <- unlist(csq_list)
  csq_split <- str_split_fixed(csq_flat, "\\|", n = length(csq_fields))
  csq_df <- as.data.frame(csq_split, stringsAsFactors = FALSE)
  colnames(csq_df) <- csq_fields
  
  # Expand fix_df to match (one row per transcript, repeating variant info —
  # including the sample_DP/AD/VAF columns just added, since they're now part of fix_df)
  fix_df_expanded <- fix_df[variant_idx, ]
  
  full_df <- cbind(fix_df_expanded, csq_df)
  
  # Write out table
  write.table(full_df, file = paste0(output.dir, temp.id, "_final.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
  
})


