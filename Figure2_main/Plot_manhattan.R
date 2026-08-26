# Load required library
library(qqman)
input_path <- "NAME_assoc.txt"

# Read input file
gemma_data <- read.table(input_path, header = TRUE)

# Ensure numeric columns are properly formatted
gemma_data$chr <- as.numeric(gemma_data$chr)
gemma_data$ps <- as.numeric(gemma_data$ps)
gemma_data$p_score <- as.numeric(gemma_data$p_score)

# ---------- Bonferroni thresholds ----------
total_snps <- 1691704 #[from GWAS file]
bonf <- -log10(0.05 / total_snps)
suggestive <- -log10(1e-5)

# Function to create Manhattan plot
plot_manhattan <- function(data, p_col, output_name) {
  # Prepare data
  plot_data <- data.frame(
    CHR = data$chr,
    BP = data$ps,
    P = data[[p_col]],
    SNP = data$rs
  )
  
  # Remove NA values
  plot_data <- na.omit(plot_data)
  
  # Plot and save
  png(output_name, width = 1200, height = 600, res = 150)
  manhattan(plot_data, 
            main = paste("Manhattan Plot -", p_col), 
            col = c("black", "grey"), 
            cex = 0.7, 
            ylim = c(0, max(-log10(plot_data$P), na.rm = TRUE) + 1),
            genomewideline = -log10(bonferroni_threshold),  # study-specific Bonferroni
            suggestiveline = -log10(suggestive_threshold))  # suggestive threshold
  dev.off()
}

# Plot for p_wald
plot_manhattan(gemma_data, "p_wald", "manhattan_Theis_BlckL_SL_p_wald.png")

# Plot for p_lrt
plot_manhattan(gemma_data, "p_lrt", "manhattan_Theis_BlckL_SL_p_lrt.png")

# Plot for p_score
plot_manhattan(gemma_data, "p_score", "manhattan_Theis_BlckL_SL_p_score.png")