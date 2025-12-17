library(tidyverse)
library(dplyr)
library(ggpubr)
library(readr)

# Input files and parameters
input_file <- "only_SL_as_covar/GT_top_100snps_THEIS_SL.txt"
phenotype_file <- "Masoko_filtered_egg_spot_theis.csv"
sample_names <- "samples_for_theis.txt"
Phenotype <- "Eggspot number (Theis)"

# Load base data once
data <- read.table(input_file, header = FALSE)
samples <- read.table(sample_names, header = FALSE)
colnames(data)[5:ncol(data)] <- samples$V1

pheno <- read.csv(phenotype_file, header = TRUE)

# Extract all unique SNPs from the input file
data <- as.data.frame(data)  # ensure it's a data frame

# Select V1 and V2
df_selected <- data[, c("V1", "V2")]

# Get distinct rows
df_distinct <- unique(df_selected)

# Rename columns
colnames(df_distinct) <- c("chr", "pos")

# Order by chr and pos
snp_list <- df_distinct[order(df_distinct$chr, df_distinct$pos), ]


cat("Found", nrow(snp_list), "unique SNPs in the input file\n")
print(head(snp_list, 10))  # Show first 10 SNPs

# Initialize results storage
results_summary <- data.frame(
  CHR = character(),
  POS = numeric(),
  REF = character(),
  ALT = character(),
  R2_null = numeric(),
  R2_full = numeric(),
  Delta_R2 = numeric(),
  Percent_Variance = numeric(),
  ANOVA_F = numeric(),
  ANOVA_p = numeric(),
  stringsAsFactors = FALSE
)

# Create output directory for plots
dir.create("SNP_plots", showWarnings = FALSE)

# Loop through each SNP
for(i in 1:nrow(snp_list)) {
  
  focal_chr <- snp_list$chr[i]
  focal_snp <- snp_list$pos[i]
  
  cat("Processing SNP", i, ":", focal_chr, ":", focal_snp, "\n")
  
  # Filter data for current SNP
  long_data <- data %>%
    pivot_longer(cols = -c(V1, V2, V3, V4), names_to = "sampleID", values_to = "Genotype") %>%
    filter(V1 == focal_chr & V2 == focal_snp)
  
  # Skip if no data for this SNP
  if(nrow(long_data) == 0) {
    cat("No data found for", focal_chr, ":", focal_snp, "- skipping\n")
    next
  }
  
  # Merge with phenotype data
  long_data <- long_data %>%
    left_join(pheno %>% select(sampleID, EggspotNumberTheis, PCA_Group, SL, pc1, pc2), by = "sampleID")
  
  colnames(long_data) <- c("CHR", "POS", "REF", "ALT", "SAMPLEID", "GENOTYPE", "PHENO", "PCA_GROUP", "SL", "PC1", "PC2")
  
  # Convert genotype to string
  long_data$GENOTYPE <- ifelse(
    long_data$GENOTYPE == "0/0", paste0(long_data$REF, long_data$REF),
    ifelse(long_data$GENOTYPE == "0/1", paste0(long_data$REF, long_data$ALT), paste0(long_data$ALT, long_data$ALT))
  )
  
  # Assign ecotype labels
  long_data$Group_Label <- recode(long_data$PCA_GROUP,
                                  "High" = "Benthic",
                                  "Middle" = "Intermediate",
                                  "Low" = "Littoral")
  
  # Convert variables
  long_data$PHENO <- as.numeric(long_data$PHENO)
  long_data$GENOTYPE <- factor(long_data$GENOTYPE, levels = sort(unique(long_data$GENOTYPE)))
  
  # Genotype dosage coding
  long_data$GENO_NUM <- case_when(
    long_data$GENOTYPE == paste0(long_data$REF[1], long_data$REF[1]) ~ 0,
    long_data$GENOTYPE == paste0(long_data$REF[1], long_data$ALT[1]) ~ 1,
    long_data$GENOTYPE == paste0(long_data$ALT[1], long_data$REF[1]) ~ 1,
    long_data$GENOTYPE == paste0(long_data$ALT[1], long_data$ALT[1]) ~ 2,
    TRUE ~ NA_real_
  )
  
  # Remove incomplete cases
  long_data <- long_data %>% drop_na(PHENO, GENO_NUM, SL, PC1, PC2)
  
  # Skip if insufficient data
  if(nrow(long_data) < 10) {
    cat("Insufficient data for", focal_chr, ":", focal_snp, "- skipping\n")
    next
  }
  
  # Fit models
  model_null <- lm(PHENO ~ SL + PC1 + PC2, data = long_data)
  model_full <- lm(PHENO ~ GENO_NUM + SL + PC1 + PC2, data = long_data)
  
  # Calculate variance explained
  r2_null <- summary(model_null)$r.squared
  r2_full <- summary(model_full)$r.squared
  delta_r2 <- r2_full - r2_null
  
  # ANOVA
  anova_result <- anova(model_null, model_full)
  
  # Store results
  results_summary <- rbind(results_summary, data.frame(
    CHR = focal_chr,
    POS = focal_snp,
    REF = long_data$REF[1],
    ALT = long_data$ALT[1],
    R2_null = r2_null,
    R2_full = r2_full,
    Delta_R2 = delta_r2,
    Percent_Variance = round(delta_r2 * 100, 2),
    ANOVA_F = anova_result$F[2],
    ANOVA_p = anova_result$`Pr(>F)`[2]
  ))
  
  # Generate comparison list for plot
  comparison_list <- combn(levels(long_data$GENOTYPE), 2, simplify = FALSE)
  
  # Create plot using your preferred style
  p <- ggplot(long_data, aes(x = GENOTYPE, y = PHENO)) +
    
    # Violin plot
    geom_violin(
      width = 0.8,
      alpha = 0.3,
      color = "black",
      linewidth = 0.4,
      trim = FALSE
    ) +
    
    # Boxplot overlay (narrow)
    geom_boxplot(
      width = 0.1,
      alpha = 0.6,
      outlier.shape = NA,
      color = "black",
      linewidth = 0.4
    ) +
    
    # Jittered data points
    geom_jitter(
      aes(color = Group_Label),
      width = 0.07,
      height = 0,
      size = 2,
      alpha = 0.2,
      stroke = 0.7
    ) +
    
    # Stats
    stat_compare_means(
      method = "t.test",
      comparisons = comparison_list,
      label = "p.signif",
      hide.ns = TRUE
    ) +
    
    # Custom color scales for ecotype
    scale_color_manual(
      values = c("Benthic" = "steelblue2", "Intermediate" = "grey", "Littoral" = "goldenrod2"),
      name = "Ecotype"
    ) +
    
    # Labels and theme
    labs(
      x = paste("Genotype", focal_chr, ":", focal_snp),
      y = Phenotype,
      title = ""
    ) +
    theme_classic(base_size = 16) +
    theme(
      legend.position = "right",
      axis.text = element_text(color = "black"),
      axis.line = element_line(color = "black", linewidth = 0.8),
      axis.ticks = element_line(color = "black")
    )
  
  # Save plot
  plot_filename <- paste0("SNP_plots/SNP_", sprintf("%02d", i), "_", focal_chr, "_", focal_snp, "_plot.pdf")
  ggsave(plot_filename, plot = p, width = 6, height = 5, dpi = 300)
  
  cat("Variance explained by SNP after accounting for covariates:", round(delta_r2 * 100, 2), "%\n")
  cat("Plot saved as:", plot_filename, "\n\n")
}

# Save results summary
write.csv(results_summary, "SNP_analysis_results_summary.csv", row.names = FALSE)
