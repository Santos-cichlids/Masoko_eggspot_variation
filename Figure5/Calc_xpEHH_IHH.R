# Load libraries
library(rehh, lib.loc = "~/rds/rds-durbin-group-8b3VcZwY7rY/projects/cichlid/suhaib/mbuna_vcf/vcf_attempt3/vcf-filtering_245_males/GWAS/R_lib/")
library(tibble, lib.loc = "~/rds/rds-durbin-group-8b3VcZwY7rY/projects/cichlid/suhaib/mbuna_vcf/vcf_attempt3/vcf-filtering_245_males/GWAS/R_lib/")
library(R.utils, lib.loc= "~/rds/rds-durbin-group-8b3VcZwY7rY/projects/cichlid/suhaib/mbuna_vcf/vcf_attempt3/vcf-filtering_245_males/GWAS/R_lib/")
library(data.table, lib.loc= "~/rds/rds-durbin-group-8b3VcZwY7rY/projects/cichlid/suhaib/mbuna_vcf/vcf_attempt3/vcf-filtering_245_males/GWAS/R_lib/")

# Get chromosome number from command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("Usage: Rscript calc_EHH.R <chromosome_number>")
(base) [sfy23@login-q-3 IHH]$ cat calc_IHH.R
# Load libraries
library(rehh, lib.loc = "~/rds/rds-durbin-group-8b3VcZwY7rY/projects/cichlid/suhaib/mbuna_vcf/vcf_attempt3/vcf-filtering_245_males/GWAS/R_lib/")
library(tibble, lib.loc = "~/rds/rds-durbin-group-8b3VcZwY7rY/projects/cichlid/suhaib/mbuna_vcf/vcf_attempt3/vcf-filtering_245_males/GWAS/R_lib/")
library(R.utils, lib.loc= "~/rds/rds-durbin-group-8b3VcZwY7rY/projects/cichlid/suhaib/mbuna_vcf/vcf_attempt3/vcf-filtering_245_males/GWAS/R_lib/")
library(data.table, lib.loc= "~/rds/rds-durbin-group-8b3VcZwY7rY/projects/cichlid/suhaib/mbuna_vcf/vcf_attempt3/vcf-filtering_245_males/GWAS/R_lib/")

# Get chromosome number from command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("Usage: Rscript calc_EHH.R <chromosome_number>")
}
chromosome <- args[1]

# Create file paths using the chromosome variable
benth_file <- paste0("../Masoko_chr", chromosome, "_benthic.vcf.gz")
litt_file <- paste0("../Masoko_chr", chromosome, "_littoral.vcf.gz")
whole_pop_file <- paste0("Masoko_chr", chromosome, "_ALL.vcf.gz")


# ========== BENTHIC POPULATION ==========

Benth_hh <- data2haplohh(hap_file = benth_file, 
                         polarize_vcf = FALSE,
                         vcf_reader = "data.table")

Benth_hh_f <- subset(Benth_hh, min_maf = 0.05)

# Use scan_hh to calculate integrated EHH (iHH) for all markers
Benth_scan <- scan_hh(Benth_hh_f, polarized = FALSE, 
                      discard_integration_at_border = FALSE)

# Extract the EHH data from scan results
Benth_ehh_df <- as_tibble(Benth_scan)

benth_ehh_file <- paste0("./EHH_Benthic_chr", chromosome, ".tsv")
write.table(Benth_ehh_df, file = benth_ehh_file, sep = "\t", 
            quote = FALSE, row.names = FALSE)


# Clean up
rm(Benth_hh, Benth_hh_f, Benth_scan, Benth_ehh_df)
gc()

# ========== LITTORAL POPULATION ==========
cat("Processing Littoral population...\n")

Litt_hh <- data2haplohh(hap_file = litt_file, 
                        polarize_vcf = FALSE,
                        vcf_reader = "data.table")

Litt_hh_f <- subset(Litt_hh, min_maf = 0.05)

# Use scan_hh to calculate integrated EHH (iHH) for all markers
Litt_scan <- scan_hh(Litt_hh_f, polarized = FALSE, 
                     discard_integration_at_border = FALSE)

# Extract the EHH data from scan results
Litt_ehh_df <- as_tibble(Litt_scan)

litt_ehh_file <- paste0("./EHH_Littoral_chr", chromosome, ".tsv")
write.table(Litt_ehh_df, file = litt_ehh_file, sep = "\t", 
            quote = FALSE, row.names = FALSE)


# Clean up
rm(Litt_hh, Litt_hh_f, Litt_scan, Litt_ehh_df)
gc()

# ========== WHOLE POPULATION ==========

Whole_hh <- data2haplohh(hap_file = whole_pop_file, 
                         polarize_vcf = FALSE,
                         vcf_reader = "data.table")

Whole_hh_f <- subset(Whole_hh, min_maf = 0.05)

# Use scan_hh to calculate integrated EHH (iHH) for all markers
Whole_scan <- scan_hh(Whole_hh_f, polarized = FALSE, 
                      discard_integration_at_border = FALSE)

# Extract the EHH data from scan results
Whole_ehh_df <- as_tibble(Whole_scan)

whole_ehh_file <- paste0("./EHH_Whole_chr", chromosome, ".tsv")
write.table(Whole_ehh_df, file = whole_ehh_file, sep = "\t", 
            quote = FALSE, row.names = FALSE)

cat(paste0("  ✓ Whole population EHH written to: ", whole_ehh_file, "\n"))
cat(paste0("    Total markers: ", nrow(Whole_ehh_df), "\n\n"))

# Clean up
rm(Whole_hh, Whole_hh_f, Whole_scan, Whole_ehh_df)
gc()
