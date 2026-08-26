library(tidyverse)
data<- read.table("Masoko_biallelic_pass_snps_mean6600_DP30_Q20_GQ30_MM10.GenProbsBin.covar", header=F)
pheno<-read.csv("Masoko.csv", header=TRUE) %>% select(PCA_Group) %>% filter(PCA_Group!="")
data <- cbind (pheno, data)

# Example with prcomp()
pca_result <- prcomp(data [, 2:649], scale. = TRUE)

# Get percentage variance explained
pca_var <- pca_result$sdev^2
pca_var_perc <- round(100 * pca_var / sum(pca_var), 1)

# Prepare dataframe for plotting
pca_df <- as.data.frame(pca_result$x)
pca_df$PCA_Group <- data$PCA_Group

# Recode PCA_Group to readable labels
pca_df$Group_Label <- recode(pca_df$PCA_Group,
                             "High" = "Benthic",
                             "Middle" = "Intermediate",
                             "Low" = "Littoral")

# Set factor levels in desired order (Littoral, Intermediate, Benthic)
pca_df$Group_Label <- factor(pca_df$Group_Label, 
                             levels = c("Littoral", "Intermediate", "Benthic"))

# Plot PC1 vs PC2 with percentage and group labels
p <- ggplot(pca_df, aes(x = -(PC1), y = -(PC2), color = Group_Label)) +
  geom_point(size = 2, alpha = 0.8) +
  labs(
    title = "",
    x = paste0("PC1 (", pca_var_perc[1], "%)"),
    y = paste0("PC2 (", pca_var_perc[2], "%)"),
    color = "Ecotype"
  ) +
  theme_classic(base_size = 16) +
  theme(
    axis.text = element_text(color = "black"),
    axis.line = element_line(color = "black", linewidth = 0.8),å
    axis.ticks = element_line(color = "black")
  ) +
  scale_color_manual(
    values = c("Benthic" = "steelblue2", "Intermediate" = "grey", "Littoral" = "goldenrod2"),
    name = "Ecotype"
  )

# Save high-resolution output (publication quality)
ggsave("PCA_publication_Plot_648.pdf", plot = p, width = 6, height = 5, dpi = 300)
