#----------Plotting using Aaron's dataset---------

# Load required packages
library(DESeq2)
library(tidyverse)
library(ggpubr)
library(rstatix)

# Gene mapping - Set 1
gene_mapping <- data.frame(
  ensembl = c("ENSACLG00000019781", "ENSACLG00000017026", "ENSACLG00000016986", 
              "ENSACLG00000019614", "ENSACLG00000014805", "ENSACLG00000001747",
              "ENSACLG00000019599", "ENSACLG00000019578", "ENSACLG00000005191",
              "ENSACLG00000004021", "ENSACLG00000024080", "ENSACLG00000024228",
              "ENSACLG00000024190"),
  gene_name = c("gnas", "tmem254", "rhoq", "prelid3b", "ntrk3a", "atp8b1",
                "aurka", "slc32a1", "ptpa", "kcnh2", "herc2", "gabrg3",
                "oca2")
)

# Gene mapping - Set 2
gene_mapping_2 <- data.frame(
  ensembl = c("ENSACLG00000013104", "ENSACLG00000019449", "ENSACLG00000003533", 
              "ENSACLG00000001370", "ENSACLG00000001390", "ENSACLG00000001395",
              "ENSACLG00000003518", "ENSACLG00000015839"),
  gene_name = c("tbx3a", "fgf10a", "hoxa13b", "hoxa13a", "evx1", "hibadha",
                "hibadhb", "kita")
)

# =============================================================================
# SEX ANALYSIS WITH DESEQ2
# =============================================================================

# Load and prepare sex data
count_data_raw <- read.csv("../ACAL_gene_counts_Extracted.csv", header = TRUE)
count_data <- count_data_raw %>%
  filter(ensembl %in% gene_mapping$ensembl)

# Create count matrix
counts_mat <- count_data %>%
  select(-ensembl) %>%
  as.matrix()
rownames(counts_mat) <- count_data$ensembl

# Create metadata
sample_names <- colnames(counts_mat)
col_data <- data.frame(
  sample = sample_names,
  group = ifelse(grepl("^F", sample_names), "F", "M")
) %>%
  filter(group %in% c("F", "M"))
rownames(col_data) <- col_data$sample

# Filter count matrix
counts_mat <- counts_mat[, col_data$sample]
col_data$group <- factor(col_data$group, levels = c("F", "M"))

# DESeq2 analysis
dds_sex <- DESeqDataSetFromMatrix(countData = counts_mat,
                                  colData = col_data,
                                  design = ~ group)
dds_sex <- DESeq(dds_sex)

# Get normalized counts for log fold change calculation
normalized_counts_sex <- counts(dds_sex, normalized = TRUE)

# Calculate log fold changes for each sample relative to female mean
female_samples <- col_data$sample[col_data$group == "F"]
male_samples <- col_data$sample[col_data$group == "M"]

# Calculate mean expression for females (reference group)
female_means <- rowMeans(normalized_counts_sex[, female_samples, drop = FALSE])

# Calculate log fold changes for each sample
sex_logfc_individual <- normalized_counts_sex %>%
  as.data.frame() %>%
  rownames_to_column("ensembl") %>%
  pivot_longer(-ensembl, names_to = "sample", values_to = "normalized_count") %>%
  left_join(col_data, by = "sample") %>%
  left_join(gene_mapping, by = "ensembl") %>%
  left_join(
    data.frame(ensembl = names(female_means), female_mean = female_means),
    by = "ensembl"
  ) %>%
  mutate(
    log2_fold_change = log2((normalized_count + 1) / (female_mean + 1))
  )

# Get DESeq2 results for significance testing
res_sex <- results(dds_sex, contrast = c("group", "M", "F"))
sex_significance <- as.data.frame(res_sex) %>%
  rownames_to_column("ensembl") %>%
  mutate(
    significance = case_when(
      padj < 0.001 ~ "***",
      padj < 0.01 ~ "**", 
      padj < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  ) %>%
  select(ensembl, significance)

# =============================================================================
# DOMAIN ANALYSIS WITH DESEQ2
# =============================================================================

# Load and prepare domain data
domain_data <- read.csv("../Domain_expression.csv", header = TRUE) %>%
  filter(ensembl %in% gene_mapping$ensembl)

# Create domain count matrix
domain_counts_mat <- domain_data %>%
  select(-ensembl) %>%
  as.matrix()
rownames(domain_counts_mat) <- domain_data$ensembl

# Create domain metadata
domain_sample_names <- colnames(domain_counts_mat)
domain_col_data <- data.frame(
  sample = domain_sample_names,
  group = case_when(
    grepl("AT$", domain_sample_names) ~ "AT",
    grepl("M$", domain_sample_names) ~ "M",
    grepl("PT$", domain_sample_names) ~ "PT"
  )
) %>%
  filter(!is.na(group))
rownames(domain_col_data) <- domain_col_data$sample

# Filter domain count matrix
domain_counts_mat <- domain_counts_mat[, domain_col_data$sample]
domain_col_data$group <- factor(domain_col_data$group, levels = c("AT", "M", "PT"))

# DESeq2 analysis for domains
dds_domain <- DESeqDataSetFromMatrix(countData = domain_counts_mat,
                                     colData = domain_col_data,
                                     design = ~ group)
dds_domain <- DESeq(dds_domain)

# Get normalized counts for domain log fold change calculation
normalized_counts_domain <- counts(dds_domain, normalized = TRUE)

# Calculate log fold changes for each sample relative to AT mean (reference group)
AT_samples <- domain_col_data$sample[domain_col_data$group == "AT"]
M_samples <- domain_col_data$sample[domain_col_data$group == "M"]
PT_samples <- domain_col_data$sample[domain_col_data$group == "PT"]

# Calculate mean expression for AT (reference group)
AT_means <- rowMeans(normalized_counts_domain[, AT_samples, drop = FALSE])

# Calculate log fold changes for each sample
domain_logfc_individual <- normalized_counts_domain %>%
  as.data.frame() %>%
  rownames_to_column("ensembl") %>%
  pivot_longer(-ensembl, names_to = "sample", values_to = "normalized_count") %>%
  left_join(domain_col_data, by = "sample") %>%
  left_join(gene_mapping, by = "ensembl") %>%
  left_join(
    data.frame(ensembl = names(AT_means), AT_mean = AT_means),
    by = "ensembl"
  ) %>%
  mutate(
    log2_fold_change = log2((normalized_count + 1) / (AT_mean + 1))
  )

# Get domain pairwise comparison results for significance
res_AT_M <- results(dds_domain, contrast = c("group", "AT", "M"))
res_AT_PT <- results(dds_domain, contrast = c("group", "AT", "PT"))
res_M_PT <- results(dds_domain, contrast = c("group", "M", "PT"))

# Create significance lookup
domain_significance <- list(
  "AT-M" = as.data.frame(res_AT_M) %>% rownames_to_column("ensembl"),
  "AT-PT" = as.data.frame(res_AT_PT) %>% rownames_to_column("ensembl"),
  "M-PT" = as.data.frame(res_M_PT) %>% rownames_to_column("ensembl")
)

# =============================================================================
# PLOTTING LOG FOLD CHANGES
# =============================================================================

# Custom theme
custom_theme <- theme_classic(base_size = 14) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "italic", size = 12),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.position = "bottom",
    legend.title = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# Color palettes
sex_colors <- c("F" = "goldenrod3", "M" = "steelblue3")
domain_colors <- c("AT" = "goldenrod4", "M" = "goldenrod2", "PT" = "goldenrod3")

# SEX BOXPLOT WITH LOG FOLD CHANGES
sex_plot_data <- sex_logfc_individual %>%
  left_join(sex_significance, by = "ensembl")

p_sex <- ggplot(sex_plot_data, aes(x = group, y = log2_fold_change, fill = group)) +
  facet_wrap(~ gene_name, scales = "free_y", ncol = 4) +
  geom_boxplot(outlier.shape = NA, color = "black", alpha = 0.8) +
  geom_jitter(width = 0.2, size = 1.5, alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", alpha = 0.7) +
  stat_compare_means(
    method = "t.test", 
    label = "p.signif",
    comparisons = list(c("F", "M")),
    size = 4
  ) +
  scale_fill_manual(values = sex_colors) +
  labs(title = "Log2 Fold Change by Sex (relative to female mean)", 
       x = "Sex", 
       y = "Log2 Fold Change") +
  custom_theme

print(p_sex)

# DOMAIN BOXPLOT WITH LOG FOLD CHANGES
p_domain <- ggplot(domain_logfc_individual, aes(x = group, y = log2_fold_change, fill = group)) +
  facet_wrap(~ gene_name, scales = "free_y", ncol = 4) +
  geom_boxplot(outlier.shape = NA, color = "black", alpha = 0.8) +
  geom_jitter(width = 0.2, size = 1.5, alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", alpha = 0.7) +
  stat_compare_means(
    method = "t.test",
    comparisons = list(c("AT", "M"), c("AT", "PT"), c("M", "PT")),
    label = "p.signif",
    tip.length = 0.01,
    step.increase = 0.1,
    size = 3
  ) +
  scale_fill_manual(values = domain_colors) +
  labs(title = "Log2 Fold Change by Domain (relative to AT mean)", 
       x = "Domain", 
       y = "Log2 Fold Change") +
  custom_theme

print(p_domain)

# =============================================================================
# LOG NORMALIZED COUNTS PLOTS (SAME FORMAT AS LOG FOLD CHANGE)
# =============================================================================

# SEX LOG NORMALIZED COUNTS PLOT
sex_lognorm_data <- sex_logfc_individual %>%
  mutate(log_norm_counts = log2(normalized_count + 1))

p_sex_lognorm <- ggplot(sex_lognorm_data, aes(x = group, y = log_norm_counts, fill = group)) +
  facet_wrap(~ gene_name, scales = "free_y", ncol = 4) +
  geom_boxplot(outlier.shape = NA, color = "black", alpha = 0.8) +
  geom_jitter(width = 0.2, size = 1.5, alpha = 0.5) +
  stat_compare_means(
    method = "t.test", 
    label = "p.signif",
    comparisons = list(c("F", "M")),
    size = 4
  ) +
  scale_fill_manual(values = sex_colors) +
  labs(title = "Log2 Normalized Counts by Sex", 
       x = "Sex", 
       y = "Log2 Normalized Counts") +
  custom_theme

print(p_sex_lognorm)

# DOMAIN LOG NORMALIZED COUNTS PLOT
domain_lognorm_data <- domain_logfc_individual %>%
  mutate(log_norm_counts = log2(normalized_count + 1))

p_domain_lognorm <- ggplot(domain_lognorm_data, aes(x = group, y = log_norm_counts, fill = group)) +
  facet_wrap(~ gene_name, scales = "free_y", ncol = 4) +
  geom_boxplot(outlier.shape = NA, color = "black", alpha = 0.8) +
  geom_jitter(width = 0.2, size = 1.5, alpha = 0.5) +
  stat_compare_means(
    method = "t.test",
    comparisons = list(c("AT", "M"), c("AT", "PT"), c("M", "PT")),
    label = "p.signif",
    tip.length = 0.01,
    step.increase = 0.1,
    size = 3
  ) +
  scale_fill_manual(values = domain_colors) +
  labs(title = "Log2 Normalized Counts by Domain", 
       x = "Domain", 
       y = "Log2 Normalized Counts") +
  custom_theme

print(p_domain_lognorm)

# =============================================================================
# GENE SET 2 ANALYSIS
# =============================================================================

cat("\n=== ANALYZING GENE SET 2 ===\n")

# =============================================================================
# SEX ANALYSIS WITH DESEQ2 - GENE SET 2
# =============================================================================

# Load and prepare sex data for gene set 2
count_data_2 <- count_data_raw %>%
  filter(ensembl %in% gene_mapping_2$ensembl)

# Create count matrix for gene set 2
counts_mat_2 <- count_data_2 %>%
  select(-ensembl) %>%
  as.matrix()
rownames(counts_mat_2) <- count_data_2$ensembl

# Filter count matrix for gene set 2
counts_mat_2 <- counts_mat_2[, col_data$sample]

# DESeq2 analysis for gene set 2
dds_sex_2 <- DESeqDataSetFromMatrix(countData = counts_mat_2,
                                    colData = col_data,
                                    design = ~ group)
dds_sex_2 <- DESeq(dds_sex_2)

# Get normalized counts for log fold change calculation - gene set 2
normalized_counts_sex_2 <- counts(dds_sex_2, normalized = TRUE)

# Calculate mean expression for females (reference group) - gene set 2
female_means_2 <- rowMeans(normalized_counts_sex_2[, female_samples, drop = FALSE])

# Calculate log fold changes for each sample - gene set 2
sex_logfc_individual_2 <- normalized_counts_sex_2 %>%
  as.data.frame() %>%
  rownames_to_column("ensembl") %>%
  pivot_longer(-ensembl, names_to = "sample", values_to = "normalized_count") %>%
  left_join(col_data, by = "sample") %>%
  left_join(gene_mapping_2, by = "ensembl") %>%
  left_join(
    data.frame(ensembl = names(female_means_2), female_mean = female_means_2),
    by = "ensembl"
  ) %>%
  mutate(
    log2_fold_change = log2((normalized_count + 1) / (female_mean + 1))
  )

# Get DESeq2 results for significance testing - gene set 2
res_sex_2 <- results(dds_sex_2, contrast = c("group", "M", "F"))
sex_significance_2 <- as.data.frame(res_sex_2) %>%
  rownames_to_column("ensembl") %>%
  mutate(
    significance = case_when(
      padj < 0.001 ~ "***",
      padj < 0.01 ~ "**", 
      padj < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  ) %>%
  select(ensembl, significance)

# =============================================================================
# DOMAIN ANALYSIS WITH DESEQ2 - GENE SET 2
# =============================================================================

# Load and prepare domain data for gene set 2
domain_data_2 <- read.csv("../Domain_expression.csv", header = TRUE) %>%
  filter(ensembl %in% gene_mapping_2$ensembl)

# Create domain count matrix for gene set 2
domain_counts_mat_2 <- domain_data_2 %>%
  select(-ensembl) %>%
  as.matrix()
rownames(domain_counts_mat_2) <- domain_data_2$ensembl

# Filter domain count matrix for gene set 2
domain_counts_mat_2 <- domain_counts_mat_2[, domain_col_data$sample]

# DESeq2 analysis for domains - gene set 2
dds_domain_2 <- DESeqDataSetFromMatrix(countData = domain_counts_mat_2,
                                       colData = domain_col_data,
                                       design = ~ group)
dds_domain_2 <- DESeq(dds_domain_2)

# Get normalized counts for domain log fold change calculation - gene set 2
normalized_counts_domain_2 <- counts(dds_domain_2, normalized = TRUE)

# Calculate mean expression for AT (reference group) - gene set 2
AT_means_2 <- rowMeans(normalized_counts_domain_2[, AT_samples, drop = FALSE])

# Calculate log fold changes for each sample - gene set 2
domain_logfc_individual_2 <- normalized_counts_domain_2 %>%
  as.data.frame() %>%
  rownames_to_column("ensembl") %>%
  pivot_longer(-ensembl, names_to = "sample", values_to = "normalized_count") %>%
  left_join(domain_col_data, by = "sample") %>%
  left_join(gene_mapping_2, by = "ensembl") %>%
  left_join(
    data.frame(ensembl = names(AT_means_2), AT_mean = AT_means_2),
    by = "ensembl"
  ) %>%
  mutate(
    log2_fold_change = log2((normalized_count + 1) / (AT_mean + 1))
  )

# Get domain pairwise comparison results for significance - gene set 2
res_AT_M_2 <- results(dds_domain_2, contrast = c("group", "AT", "M"))
res_AT_PT_2 <- results(dds_domain_2, contrast = c("group", "AT", "PT"))
res_M_PT_2 <- results(dds_domain_2, contrast = c("group", "M", "PT"))

# =============================================================================
# PLOTTING GENE SET 2
# =============================================================================

# SEX BOXPLOT WITH LOG FOLD CHANGES - GENE SET 2
sex_plot_data_2 <- sex_logfc_individual_2 %>%
  left_join(sex_significance_2, by = "ensembl")

p_sex_2 <- ggplot(sex_plot_data_2, aes(x = group, y = log2_fold_change, fill = group)) +
  facet_wrap(~ gene_name, scales = "free_y", ncol = 4) +
  geom_boxplot(outlier.shape = NA, color = "black", alpha = 0.8) +
  geom_jitter(width = 0.2, size = 1.5, alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", alpha = 0.7) +
  stat_compare_means(
    method = "t.test", 
    label = "p.signif",
    comparisons = list(c("F", "M")),
    size = 4
  ) +
  scale_fill_manual(values = sex_colors) +
  labs(title = "Log2 Fold Change by Sex - Gene Set 2 (relative to female mean)", 
       x = "Sex", 
       y = "Log2 Fold Change") +
  custom_theme

print(p_sex_2)

# DOMAIN BOXPLOT WITH LOG FOLD CHANGES - GENE SET 2
p_domain_2 <- ggplot(domain_logfc_individual_2, aes(x = group, y = log2_fold_change, fill = group)) +
  facet_wrap(~ gene_name, scales = "free_y", ncol = 4) +
  geom_boxplot(outlier.shape = NA, color = "black", alpha = 0.8) +
  geom_jitter(width = 0.2, size = 1.5, alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", alpha = 0.7) +
  stat_compare_means(
    method = "t.test",
    comparisons = list(c("AT", "M"), c("AT", "PT"), c("M", "PT")),
    label = "p.signif",
    tip.length = 0.01,
    step.increase = 0.1,
    size = 3
  ) +
  scale_fill_manual(values = domain_colors) +
  labs(title = "Log2 Fold Change by Domain - Gene Set 2 (relative to AT mean)", 
       x = "Domain", 
       y = "Log2 Fold Change") +
  custom_theme

print(p_domain_2)

# SEX LOG NORMALIZED COUNTS PLOT - GENE SET 2
sex_lognorm_data_2 <- sex_logfc_individual_2 %>%
  mutate(log_norm_counts = log2(normalized_count + 1))

p_sex_lognorm_2 <- ggplot(sex_lognorm_data_2, aes(x = group, y = log_norm_counts, fill = group)) +
  facet_wrap(~ gene_name, scales = "free_y", ncol = 4) +
  geom_boxplot(outlier.shape = NA, color = "black", alpha = 0.8) +
  geom_jitter(width = 0.2, size = 1.5, alpha = 0.5) +
  stat_compare_means(
    method = "t.test", 
    label = "p.signif",
    comparisons = list(c("F", "M")),
    size = 4
  ) +
  scale_fill_manual(values = sex_colors) +
  labs(title = "Log2 Normalized Counts by Sex - Gene Set 2", 
       x = "Sex", 
       y = "Log2 Normalized Counts") +
  custom_theme

print(p_sex_lognorm_2)

# DOMAIN LOG NORMALIZED COUNTS PLOT - GENE SET 2
domain_lognorm_data_2 <- domain_logfc_individual_2 %>%
  mutate(log_norm_counts = log2(normalized_count + 1))

p_domain_lognorm_2 <- ggplot(domain_lognorm_data_2, aes(x = group, y = log_norm_counts, fill = group)) +
  facet_wrap(~ gene_name, scales = "free_y", ncol = 4) +
  geom_boxplot(outlier.shape = NA, color = "black", alpha = 0.8) +
  geom_jitter(width = 0.2, size = 1.5, alpha = 0.5) +
  stat_compare_means(
    method = "t.test",
    comparisons = list(c("AT", "M"), c("AT", "PT"), c("M", "PT")),
    label = "p.signif",
    tip.length = 0.01,
    step.increase = 0.1,
    size = 3
  ) +
  scale_fill_manual(values = domain_colors) +
  labs(title = "Log2 Normalized Counts by Domain - Gene Set 2", 
       x = "Domain", 
       y = "Log2 Normalized Counts") +
  custom_theme

print(p_domain_lognorm_2)

# =============================================================================
# SUMMARY TABLES
# =============================================================================

# Sex comparison results
cat("\n=== SEX COMPARISON (M vs F) - DESeq2 Results ===\n")
res_sex <- results(dds_sex, contrast = c("group", "M", "F"))
sex_results <- as.data.frame(res_sex) %>%
  rownames_to_column("ensembl") %>%
  left_join(gene_mapping, by = "ensembl") %>%
  select(gene_name, log2FoldChange, pvalue, padj) %>%
  mutate(
    log2FoldChange = round(log2FoldChange, 3),
    pvalue = round(pvalue, 4),
    padj = round(padj, 4),
    significance = case_when(
      padj < 0.001 ~ "***",
      padj < 0.01 ~ "**", 
      padj < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  ) %>%
  arrange(padj)

print(sex_results)

# Domain pairwise comparison results
cat("\n=== DOMAIN PAIRWISE COMPARISONS - DESeq2 Results ===\n")

# AT vs M
AT_M_results <- as.data.frame(res_AT_M) %>%
  rownames_to_column("ensembl") %>%
  left_join(gene_mapping, by = "ensembl") %>%
  select(gene_name, log2FoldChange, pvalue, padj) %>%
  mutate(comparison = "AT_vs_M")

# AT vs PT  
AT_PT_results <- as.data.frame(res_AT_PT) %>%
  rownames_to_column("ensembl") %>%
  left_join(gene_mapping, by = "ensembl") %>%
  select(gene_name, log2FoldChange, pvalue, padj) %>%
  mutate(comparison = "AT_vs_PT")

# M vs PT
M_PT_results <- as.data.frame(res_M_PT) %>%
  rownames_to_column("ensembl") %>%
  left_join(gene_mapping, by = "ensembl") %>%
  select(gene_name, log2FoldChange, pvalue, padj) %>%
  mutate(comparison = "M_vs_PT")

# Combine all domain results
domain_results <- bind_rows(AT_M_results, AT_PT_results, M_PT_results) %>%
  mutate(
    log2FoldChange = round(log2FoldChange, 3),
    pvalue = round(pvalue, 4),
    padj = round(padj, 4),
    significance = case_when(
      padj < 0.001 ~ "***",
      padj < 0.01 ~ "**", 
      padj < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  ) %>%
  arrange(comparison, padj)

print(domain_results)

# =============================================================================
# SUMMARY TABLES - GENE SET 2
# =============================================================================

# Sex comparison results - Gene Set 2
cat("\n=== SEX COMPARISON (M vs F) - DESeq2 Results - GENE SET 2 ===\n")
sex_results_2 <- as.data.frame(res_sex_2) %>%
  rownames_to_column("ensembl") %>%
  left_join(gene_mapping_2, by = "ensembl") %>%
  select(gene_name, log2FoldChange, pvalue, padj) %>%
  mutate(
    log2FoldChange = round(log2FoldChange, 3),
    pvalue = round(pvalue, 4),
    padj = round(padj, 4),
    significance = case_when(
      padj < 0.001 ~ "***",
      padj < 0.01 ~ "**", 
      padj < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  ) %>%
  arrange(padj)

print(sex_results_2)

# Domain pairwise comparison results - Gene Set 2
cat("\n=== DOMAIN PAIRWISE COMPARISONS - DESeq2 Results - GENE SET 2 ===\n")

# AT vs M - Gene Set 2
AT_M_results_2 <- as.data.frame(res_AT_M_2) %>%
  rownames_to_column("ensembl") %>%
  left_join(gene_mapping_2, by = "ensembl") %>%
  select(gene_name, log2FoldChange, pvalue, padj) %>%
  mutate(comparison = "AT_vs_M")

# AT vs PT - Gene Set 2
AT_PT_results_2 <- as.data.frame(res_AT_PT_2) %>%
  rownames_to_column("ensembl") %>%
  left_join(gene_mapping_2, by = "ensembl") %>%
  select(gene_name, log2FoldChange, pvalue, padj) %>%
  mutate(comparison = "AT_vs_PT")

# M vs PT - Gene Set 2
M_PT_results_2 <- as.data.frame(res_M_PT_2) %>%
  rownames_to_column("ensembl") %>%
  left_join(gene_mapping_2, by = "ensembl") %>%
  select(gene_name, log2FoldChange, pvalue, padj) %>%
  mutate(comparison = "M_vs_PT")

# Combine all domain results - Gene Set 2
domain_results_2 <- bind_rows(AT_M_results_2, AT_PT_results_2, M_PT_results_2) %>%
  mutate(
    log2FoldChange = round(log2FoldChange, 3),
    pvalue = round(pvalue, 4),
    padj = round(padj, 4),
    significance = case_when(
      padj < 0.001 ~ "***",
      padj < 0.01 ~ "**", 
      padj < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  ) %>%
  arrange(comparison, padj)

print(domain_results_2)
