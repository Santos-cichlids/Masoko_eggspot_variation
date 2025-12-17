# Load libraries
library(rehh, lib.loc = "~/rds/rds-durbin-group-8b3VcZwY7rY/projects/cichlid/suhaib/mbuna_vcf/vcf_attempt3/vcf-filtering_245_males/GWAS/R_lib/")
library(tibble, lib.loc = "~/rds/rds-durbin-group-8b3VcZwY7rY/projects/cichlid/suhaib/mbuna_vcf/vcf_attempt3/vcf-filtering_245_males/GWAS/R_lib/")
library(R.utils, lib.loc= "~/rds/rds-durbin-group-8b3VcZwY7rY/projects/cichlid/suhaib/mbuna_vcf/vcf_attempt3/vcf-filtering_245_males/GWAS/R_lib/")
library(data.table, lib.loc= "~/rds/rds-durbin-group-8b3VcZwY7rY/projects/cichlid/suhaib/mbuna_vcf/vcf_attempt3/vcf-filtering_245_males/GWAS/R_lib/")

# Get chromosome number from command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("Usage: Rscript calc_iHH_BENTH_LITT.R <chromosome_number>")
}
chromosome <- args[1]

# Create file paths using the chromosome variable
benth_file <- paste0("Masoko_chr", chromosome, "_benthic.vcf.gz")
litt_file <- paste0("Masoko_chr", chromosome, "_littoral.vcf.gz")

# Check if input files exist
if (!file.exists(benth_file)) {
  stop(paste("File not found:", benth_file))
}
if (!file.exists(litt_file)) {
  stop(paste("File not found:", litt_file))
}

Benth_hh <- data2haplohh(hap_file = benth_file, polarize_vcf = FALSE)
Litt_hh  <- data2haplohh(hap_file = litt_file, polarize_vcf = FALSE)

Benth_hh_f <- subset(Benth_hh, min_maf = 0.05)
Litt_hh_f  <- subset(Litt_hh, min_maf = 0.05)

Benth_scan <- scan_hh(Benth_hh_f, polarized = FALSE, discard_integration_at_border = FALSE)
Litt_scan  <- scan_hh(Litt_hh_f, polarized = FALSE, discard_integration_at_border = FALSE)

rm(Benth_hh, Litt_hh, Benth_hh_f, Litt_hh_f)
gc()

Benth_ihh <- ihh2ihs(Benth_scan, standardize = FALSE)
Litt_ihh <- ihh2ihs(Litt_scan, standardize = FALSE)

Benth_ihh_df <- as_tibble(Benth_ihh$ihs)
Litt_ihh_df <- as_tibble(Litt_ihh$ihs)

BENTHvsLITT <- ies2xpehh(Benth_scan, Litt_scan, popname1 = "Benth", popname2 = "Litt", include_freq = TRUE)
BENTHvsLITT_df <- as_tibble(BENTHvsLITT)

rm(Benth_scan, Litt_scan)
gc()

colnames(Benth_ihh_df) <- tolower(colnames(Benth_ihh_df))
colnames(Litt_ihh_df) <- tolower(colnames(Litt_ihh_df))
colnames(BENTHvsLITT_df) <- tolower(colnames(BENTHvsLITT_df))

# Find common column names for merging
benth_cols <- colnames(Benth_ihh_df)
litt_cols <- colnames(Litt_ihh_df)
xpehh_cols <- colnames(BENTHvsLITT_df)

# Use the first two columns (typically chr and position) for merging
merge_cols <- intersect(benth_cols[1:2], xpehh_cols[1:2])

combined_results <- BENTHvsLITT_df

combined_results <- merge(combined_results, 
                         Benth_ihh_df, 
                         by = merge_cols, 
                         all.x = TRUE, 
                         suffixes = c("", "_benth"))

combined_results <- merge(combined_results, 
                         Litt_ihh_df, 
                         by = merge_cols, 
                         all.x = TRUE, 
                         suffixes = c("", "_litt"))

combined_results <- combined_results[order(combined_results$chr, combined_results$position), ]

output_file <- paste0("./HighSEQ_BENTHvsLITT_chr", chromosome, "_IHH_xpEHH_combined.tsv")

write.table(combined_results, file = output_file, sep = "\t", quote = FALSE, row.names = FALSE)