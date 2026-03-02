# Load required library
library (qqman, lib.loc="~/rds/rds-durbin-group-8b3VcZwY7rY/projects/cichlid/suhaib/mbuna_vcf/vcf_attempt3/vcf-filtering_245_males/GWAS/R_lib/")
input_path <- "./output/Masoko_GWAS_with_covariates_BlackL_SL_covar.assoc.txt"
# Read input file
gemma_data <- read.table(input_path, header = TRUE)  # Replace with your file name


# Ensure numeric columns are properly formatted
gemma_data$chr <- as.numeric(gemma_data$chr)
gemma_data$ps <- as.numeric(gemma_data$ps)
#gemma_data$p_wald <- as.numeric(gemma_data$p_wald)
#gemma_data$p_lrt <- as.numeric(gemma_data$p_lrt)
gemma_data$p_score <- as.numeric(gemma_data$p_score)

# Function to create Manhattan plot
plot_manhattan <- function(data, p_col, output_name) {
  # Prepare data
  plot_data <- data.frame(
    CHR = data$chr,
    BP = data$ps,
    P = data[[p_col]],
    SNP = data$rs  # Add SNP column if present
  )
  
  # Remove NA values
  plot_data <- na.omit(plot_data)
  
  # Plot and save
  png(output_name, width = 1200, height = 600, res = 150)
  manhattan(plot_data, 
            main = paste("Manhattan Plot -", p_col), 
            col = c("black", "grey"), cex=0.7, 
            ylim = c(0, max(-log10(plot_data$P), na.rm = TRUE) + 1))
  dev.off()
}
# Plot for p_wald
plot_manhattan(gemma_data, "p_wald", "manhattan_Theis_BlckL_SL_p_wald.pdf")

# Plot for p_lrt
plot_manhattan(gemma_data, "p_lrt", "manhattan_Theis_BlckL_SL_p_lrt.pdf")

# Plot for p_score
plot_manhattan(gemma_data, "p_score", "manhattan_Theis_BlckL_SL_p_score.pdf")




# Load required package
library (dplyr, lib.loc="~/rds/rds-durbin-group-8b3VcZwY7rY/projects/cichlid/suhaib/mbuna_vcf/vcf_attempt3/vcf-filtering_245_males/GWAS/R_lib/")

output_path <- "./top_BL_SL_candidate_loci.txt"

# Process data: sort by p_score and keep the first 100 rows
top_100 <- gemma_data %>%
  arrange(p_score) %>%   # Sort in ascending order by p_score
  slice_head(n = 100)     # Keep the first 100 rows

# Save the result to a file
write.table(top_100, file = output_path, sep = "\t", row.names = FALSE, quote = FALSE)