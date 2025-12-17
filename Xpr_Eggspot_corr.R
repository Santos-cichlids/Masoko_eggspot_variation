# Load required libraries
library(ggplot2)
library(dplyr)
library(ggpubr)
library(tidyr)

# Load and transpose depth data
dep <- read.table("HUE_topGenes_analfin_cmp.txt", header = TRUE)
rownames(dep) <- dep$GeneID
dep$GeneID <- NULL
dep_transposed <- as.data.frame(t(dep))
dep_transposed$sampleID <- rownames(dep_transposed)
rownames(dep_transposed) <- NULL
# Define gene of interest
gene_id <- "ENSACLG00000017033" ###TMEM254
################################################################################
################################################################################
################################################################################
# Load phenotype data and recode PCA_Group
# Read the CSV without a header
pheno_raw <- read.csv("../Masoko_GEMMA/Masoko_meta_master_20230331.csv", header = FALSE, stringsAsFactors = FALSE)

# Set second row as column names, remove first two rows
colnames(pheno_raw) <- pheno_raw[2, ]
pheno <- pheno_raw[-c(1, 2), ]

# Continue with selection and recoding
pheno <- pheno %>%
  select(sampleID, PCA_Group, EggspotNumberTheis, spot.hue) %>%
  mutate(PCA_Group = recode(PCA_Group,
                            "High" = "Benthic",
                            "Middle" = "Intermediate",
                            "Low" = "Littoral"),
         spot.hue = as.numeric(spot.hue),
         EggspotNumberTheis = as.numeric(EggspotNumberTheis))
################################################################################
################################################################################
################################################################################
# Merge phenotype with depth data
merged_data <- merge(pheno, dep_transposed, by = "sampleID")

# Drop rows with NA values in key variables
merged_data <- merged_data %>%
  filter(!is.na(spot.hue) & !is.na(.data[[gene_id]]))

# Define color scheme
type_colors <- c("Benthic" = "steelblue2",
                 "Intermediate" = "grey",
                 "Littoral" = "goldenrod2")

# Define population pairwise comparisons
comparisons <- list(
  c("Benthic", "Intermediate"),
  c("Intermediate", "Littoral"),
  c("Benthic", "Littoral")
)

# Plot expression boxplot with all data
ggplot(merged_data, aes(x = PCA_Group, y = .data[[gene_id]], fill = PCA_Group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, color = "black", width = 0.6) +
  geom_jitter(aes(color = PCA_Group), width = 0.15, size = 2.2, show.legend = FALSE) +
  stat_compare_means(comparisons = comparisons, method = "wilcox.test", label = "p.signif") +
  labs(
    x = "Ecotype",
    y = "Normalized expression (CPM)",
    title = "Expression of TMEM254 gene"
  ) +
  scale_fill_manual(values = type_colors) +
  scale_color_manual(values = type_colors) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none",
    axis.text = element_text(color = "black"),
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.title = element_text()
  )

# -------------------------------------------
# MODEL 1: With all samples (AXES FLIPPED)
# -------------------------------------------

gene_expr_all <- merged_data[[gene_id]]
# Now spot.hue is the response variable, gene expression is the predictor
model_all <- lm(spot.hue ~ gene_expr_all, data = merged_data)
model_all_sum <- summary(model_all)

intercept_all <- round(coef(model_all)[1], 3)
slope_all <- round(coef(model_all)[2], 3)
r2_all <- round(model_all_sum$r.squared, 3)
pval_all <- signif(model_all_sum$coefficients[2, 4], 3)

eq_label_all <- paste0("y = ", slope_all, "x + ", intercept_all, "\n",
                       "R² = ", r2_all, ", p = ", pval_all)

# Flipped axes: gene expression on x-axis, spot.hue on y-axis
ggplot(merged_data, aes(x = gene_expr_all, y = spot.hue, color = PCA_Group)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  annotate("text", 
           x = min(gene_expr_all, na.rm = TRUE), 
           y = max(merged_data$spot.hue, na.rm = TRUE), 
           label = eq_label_all, 
           hjust = -1, vjust = 1, size = 4, color = "black") +
  scale_color_manual(values = type_colors) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 5)) +
  labs(
    x = "Normalized expression (CPM)",
    y = "Egg-spot Hue",
    color = "Ecotype",
    title = "Regression with All Samples"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text = element_text(color = "black")
  )

# -------------------------------------------
# MODEL 2: Without specific outliers (AXES FLIPPED)
# -------------------------------------------

# Define the outlier sampleIDs
outlier_ids <- c("cichl187579999", "cichl187580024", "cichl187579992")

# Filter them out
Combined_filtered <- merged_data %>% filter(!sampleID %in% outlier_ids)
gene_expr_filtered <- Combined_filtered[[gene_id]]

# Fit model with flipped axes
model_filtered <- lm(spot.hue ~ gene_expr_filtered, data = Combined_filtered)
model_filtered_sum <- summary(model_filtered)

intercept_f <- round(coef(model_filtered)[1], 3)
slope_f <- round(coef(model_filtered)[2], 3)
r2_f <- round(model_filtered_sum$r.squared, 3)
pval_f <- signif(model_filtered_sum$coefficients[2, 4], 3)

eq_label_f <- paste0("y = ", slope_f, "x + ", intercept_f, "\n",
                     "R² = ", r2_f, ", p = ", pval_f)

# Flipped axes: gene expression on x-axis, spot.hue on y-axis
ggplot(Combined_filtered, aes(x = gene_expr_filtered, y = spot.hue, color = PCA_Group)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  annotate("text", 
           x = min(gene_expr_filtered, na.rm = TRUE), 
           y = max(Combined_filtered$spot.hue, na.rm = TRUE), 
           label = eq_label_f, 
           hjust = -1, vjust = 1, size = 4, color = "black") +
  scale_color_manual(values = type_colors) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 5)) +
  labs(
    x = "Normalized expression (CPM)",
    y = "Egg-spot Hue",
    color = "Ecotype",
    title = "Regression Without 3 Outliers"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text = element_text(color = "black")
  )
