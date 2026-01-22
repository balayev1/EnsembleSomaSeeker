#!/usr/bin/env R

library(dplyr)

args <- commandArgs(trailingOnly = TRUE)
# args[1:3] are files; args[4:13] are thresholds

# Snakemake objects
tsv_path <- args[2]
vcf_in   <- args[1]
vcf_out  <- args[3]

# Load TSV
tsv <- read.table(tsv_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Calculate Metrics
tsv$mq <- rowMeans(tsv[, c("nBAM_REF_MQ", "nBAM_ALT_MQ", "tBAM_REF_MQ", "tBAM_ALT_MQ")], na.rm = TRUE)
tsv$bq <- rowMeans(tsv[, c("nBAM_REF_BQ", "nBAM_ALT_BQ", "tBAM_REF_BQ", "tBAM_ALT_BQ")], na.rm = TRUE)
tsv$n_dp <- (tsv$N_ALT_FOR + tsv$N_ALT_REV + tsv$N_REF_FOR + tsv$N_REF_REV)
tsv$t_dp <- (tsv$T_ALT_FOR + tsv$T_ALT_REV + tsv$T_REF_FOR + tsv$T_REF_REV)
tsv$t_afalt <- (tsv$T_ALT_FOR + tsv$T_ALT_REV) / tsv$t_dp
tsv$n_afalt <- (tsv$N_ALT_FOR + tsv$N_ALT_REV) / tsv$n_dp

# Apply Hard Filters + Caller Consensus (NumCallers >= 2)
tsv$NumCallers <- rowSums(tsv[, c("if_MuTect", "if_Strelka", "if_VarScan2", "if_LoFreq", "MuSE_Tier")], na.rm = TRUE)
filtered_tsv <- tsv %>% filter(
    n_dp > args[4] & t_dp > args[5] & t_afalt > args[6] & n_afalt < args[7] & 
    tBAM_ALT_Concordant >= args[8] & nBAM_ALT_Concordant < args[9] &
    mq >= args[10] & bq >= args[11] & tBAM_ALT_MQ >= args[12] &
    NumCallers >= args[13])

# Create a list of positions to keep (CHROM_POS_REF_ALT)
keep_sites <- paste(filtered_tsv$CHROM, filtered_tsv$POS, filtered_tsv$REF, filtered_tsv$ALT, sep = "_")

# 4. Filter the VCF (Basic header handling)
lines <- readLines(vcf_in)
header <- lines[grepl("^#", lines)]
body <- lines[!grepl("^#", lines)]

# Split body and filter based on the keep_sites list
body_df <- as.data.frame(do.call(rbind, strsplit(body, "\t")), stringsAsFactors = FALSE)
if(nrow(body_df) > 0){
    body_df$site_id <- paste(body_df$V1, body_df$V2, body_df$V4, body_df$V5, sep = "_")
    final_body <- body[body_df$site_id %in% keep_sites]
} else {
    final_body <- character(0)
}

writeLines(c(header, final_body), vcf_out)