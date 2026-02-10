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

# Create file paths using the chromosome variable - NOW WITH _AA SUFFIX
benth_file <- paste0("Masoko_chr", chromosome, "_benthic_AA.vcf.gz")
litt_file <- paste0("Masoko_chr", chromosome, "_littoral_AA.vcf.gz")


# Create output directory if it doesn't exist
output_dir <- "Downsampled_xpEHH_IHS_POLARIZED"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Set parameters for downsampling
n_bootstrap <- 5  # Number of bootstrap iterations
n_target <- 98    # Target sample size (matching benthic population)

# Store results from each bootstrap iteration
all_results <- list()

for(boot_iter in 1:n_bootstrap) {
  
  # Set seed for reproducibility
  set.seed(boot_iter * 42)
  
  # Load data fresh for each iteration with polarization

  # Try loading with polarization, catch any errors
  tryCatch({
    Benth_hh <- data2haplohh(hap_file = benth_file, polarize_vcf = TRUE)
    Litt_hh  <- data2haplohh(hap_file = litt_file, polarize_vcf = TRUE)
    cat("Successfully loaded polarized data\n")
  }, error = function(e) {
    cat("ERROR in data2haplohh with polarize_vcf = TRUE:\n")
    cat(conditionMessage(e), "\n")
    cat("\nMake sure you've run add_AA_tag.sh to add ancestral allele info!\n")
    stop("Polarization failed. Check VCF format.")
  })
  
  # Get sample names
  benth_samples <- hap.names(Benth_hh)
  litt_samples <- hap.names(Litt_hh)
  
  # Downsample littorals to match benthics
  litt_subset <- sample(litt_samples, n_target, replace = FALSE)
  Litt_hh <- subset(Litt_hh, select.hap = litt_subset)
  

  # Filter by MAF
  Benth_hh_f <- subset(Benth_hh, min_maf = 0.05)
  Litt_hh_f  <- subset(Litt_hh, min_maf = 0.05)
  
  # Scan haplotypes with polarized = TRUE
  tryCatch({
    Benth_scan <- scan_hh(Benth_hh_f, polarized = TRUE, discard_integration_at_border = FALSE)
    Litt_scan  <- scan_hh(Litt_hh_f, polarized = TRUE, discard_integration_at_border = FALSE)
    cat("Successfully scanned polarized haplotypes\n")
  }, error = function(e) {
    cat("ERROR in scan_hh with polarized = TRUE:\n")
    cat(conditionMessage(e), "\n")
    stop("Haplotype scanning failed.")
  })
  
  # Clean up memory
  rm(Benth_hh, Litt_hh, Benth_hh_f, Litt_hh_f)
  gc()
  
  # Calculate iHS
  Benth_ihh <- ihh2ihs(Benth_scan, standardize = FALSE)
  Litt_ihh <- ihh2ihs(Litt_scan, standardize = FALSE)
  
  Benth_ihh_df <- as_tibble(Benth_ihh$ihs)
  Litt_ihh_df <- as_tibble(Litt_ihh$ihs)
  
  # Calculate XP-EHH
  BENTHvsLITT <- ies2xpehh(Benth_scan, Litt_scan, popname1 = "Benth", popname2 = "Litt", include_freq = TRUE)
  BENTHvsLITT_df <- as_tibble(BENTHvsLITT)
  
  # Clean up memory
  rm(Benth_scan, Litt_scan, Benth_ihh, Litt_ihh, BENTHvsLITT)
  gc()
  
  # Standardize column names
  colnames(Benth_ihh_df) <- tolower(colnames(Benth_ihh_df))
  colnames(Litt_ihh_df) <- tolower(colnames(Litt_ihh_df))
  colnames(BENTHvsLITT_df) <- tolower(colnames(BENTHvsLITT_df))
  
  # Find common column names for merging
  benth_cols <- colnames(Benth_ihh_df)
  litt_cols <- colnames(Litt_ihh_df)
  xpehh_cols <- colnames(BENTHvsLITT_df)
  
  # Use the first two columns (typically chr and position) for merging
  merge_cols <- intersect(benth_cols[1:2], xpehh_cols[1:2])
  
  # Merge results
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
  
  # Store this iteration's results
  all_results[[boot_iter]] <- combined_results
  
  # Clean up
  rm(Benth_ihh_df, Litt_ihh_df, BENTHvsLITT_df, combined_results)
  gc()
}


# Merge all bootstrap iterations by chr and position
# Start with the first iteration
combined_final <- all_results[[1]]
colnames(combined_final)[-(1:2)] <- paste0(colnames(combined_final)[-(1:2)], "_iter1")

# Merge subsequent iterations
for(i in 2:n_bootstrap) {
  temp_df <- all_results[[i]]
  colnames(temp_df)[-(1:2)] <- paste0(colnames(temp_df)[-(1:2)], "_iter", i)
  
  # Get merge column names (first two columns: chr and position)
  merge_cols <- colnames(combined_final)[1:2]
  
  combined_final <- merge(combined_final, temp_df, 
                         by = merge_cols, 
                         all = TRUE)
}

# Now calculate means and standard deviations for each statistic
# Identify the base statistic names (without _iter suffix)
all_cols <- colnames(combined_final)
merge_cols <- all_cols[1:2]  # chr and position
stat_cols_with_iter <- all_cols[-(1:2)]

# Extract unique statistic names
stat_names <- unique(gsub("_iter[0-9]+$", "", stat_cols_with_iter))

# Create final output dataframe starting with chr and position
final_output <- combined_final[, merge_cols, drop = FALSE]

# For each statistic, calculate mean and SD across iterations
for(stat in stat_names) {
  # Find all columns for this statistic across iterations
  iter_cols <- grep(paste0("^", stat, "_iter[0-9]+$"), all_cols, value = TRUE)
  
  if(length(iter_cols) > 0) {
    # Extract the values for these columns
    stat_values <- combined_final[, iter_cols, drop = FALSE]
    
    # Calculate mean and SD
    final_output[[stat]] <- rowMeans(stat_values, na.rm = TRUE)
    final_output[[paste0(stat, "_sd")]] <- apply(stat_values, 1, sd, na.rm = TRUE)
  }
}

# Sort by chromosome and position
final_output <- final_output[order(final_output[[merge_cols[1]]], final_output[[merge_cols[2]]]), ]

# Write averaged results to the output directory
output_file <- file.path(output_dir, paste0("HighSEQ_BENTHvsLITT_chr", chromosome, "_IHH_xpEHH_combined_downsampled_POLARIZED.tsv"))
write.table(final_output, file = output_file, sep = "\t", quote = FALSE, row.names = FALSE)
