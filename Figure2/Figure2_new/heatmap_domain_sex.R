# Load required packages
library(DESeq2)
library(tidyverse)
library(RColorBrewer)
library(ggplot2)
library(grid)
library(ggh4x)   # for facet_nested() - nested strip labels on domain heatmap
# install with: install.packages("ggh4x")

# Gene mapping - Set 1
# FIX: was fully commented out but is used in both left_join() calls below.
# Add your ensembl-to-gene-name mappings here. Genes not listed will display
# their Ensembl ID as the label instead (handled by the ifelse below).
gene_mapping <- data.frame(
  ensembl = c("ENSACLG00000017026", "ENSACLG00000026836", "ENSACLG00000026955",
              "ENSACLG00000024190", "ENSACLG00000025445", "ENSACLG00000016986"),
  gene_name = c("tmem254", "aldh9a1a.a", "tmco1", "oca2", "plxdc1", "rhoq")
)

# Function to perform k-means clustering based on group mean expression patterns
get_kmeans_cluster_order <- function(data, group_info, num_clusters = 3) {
  mat <- as.matrix(data %>% dplyr::select(-ensembl, -GeneName))
  rownames(mat) <- data$GeneName
  
  # Calculate mean expression per group (this is what we cluster on)
  group_levels <- levels(factor(group_info))
  group_means <- sapply(group_levels, function(g) {
    rowMeans(mat[, group_info == g, drop = FALSE])
  })
  rownames(group_means) <- rownames(mat)
  
  # Scale by row for clustering
  group_means_scaled <- t(scale(t(group_means)))
  
  # Remove genes with NA or infinite values
  valid_genes <- complete.cases(group_means_scaled) & 
    apply(group_means_scaled, 1, function(x) all(is.finite(x)))
  group_means_scaled <- group_means_scaled[valid_genes, ]
  
  # Perform k-means clustering
  set.seed(123)  # For reproducibility
  kmeans_result <- kmeans(group_means_scaled, centers = num_clusters, nstart = 25, iter.max = 100)
  
  # Get cluster assignments
  clusters <- kmeans_result$cluster
  
  # Calculate mean expression for each cluster across groups to order them
  cluster_profiles <- sapply(1:num_clusters, function(c) {
    colMeans(group_means_scaled[clusters == c, , drop = FALSE])
  })
  
  # Order clusters by their first group mean (or by variance)
  cluster_variance <- apply(cluster_profiles, 2, var)
  cluster_order <- order(cluster_variance, decreasing = TRUE)
  
  # Reassign cluster numbers based on ordering
  cluster_labels <- match(clusters, cluster_order)
  names(cluster_labels) <- names(clusters)
  
  # Within each cluster, order genes by their variance across groups
  gene_vars <- apply(group_means_scaled, 1, var)
  
  gene_order <- c()
  for (i in 1:num_clusters) {
    genes_in_cluster <- names(cluster_labels)[cluster_labels == i]
    cluster_gene_vars <- gene_vars[genes_in_cluster]
    ordered_genes <- genes_in_cluster[order(cluster_gene_vars, decreasing = TRUE)]
    gene_order <- c(gene_order, ordered_genes)
  }
  
  return(list(
    order = gene_order, 
    clusters = cluster_labels,
    cluster_sizes = table(cluster_labels),
    group_means = group_means[names(cluster_labels), ]
  ))
}

###########################################################################
########################### SEX DATA ######################################
###########################################################################

# Load and prepare sex data - ALL genes
count_data <- read.csv("../ACAL_gene_counts_Extracted.csv", header = TRUE)
counts_mat <- count_data %>% dplyr::select(-ensembl) %>% as.matrix()
rownames(counts_mat) <- count_data$ensembl
sample_names <- colnames(counts_mat)
col_data <- data.frame(sample = sample_names, group = ifelse(grepl("^F", sample_names), "F", "M")) %>%
  filter(group %in% c("F", "M"))
rownames(col_data) <- col_data$sample
counts_mat <- counts_mat[, col_data$sample]
col_data$group <- factor(col_data$group, levels = c("F", "M"))

# Filter genes: keep only genes with mean count >10 across samples
gene_means <- rowMeans(counts_mat)
keep_genes <- gene_means > 10
cat(sprintf("Sex data: Before filtering: %d genes\n", nrow(counts_mat)))
cat(sprintf("Sex data: Genes with mean count >10: %d genes\n", sum(keep_genes)))
cat(sprintf("Sex data: Filtering out %d genes\n", sum(!keep_genes)))
counts_mat <- counts_mat[keep_genes, ]
cat(sprintf("Sex data: After filtering: %d genes retained\n", nrow(counts_mat)))

# DESeq2 analysis for sex - use FILTERED dataset for proper normalization
dds_sex <- DESeqDataSetFromMatrix(countData = counts_mat, colData = col_data, design = ~ group)
dds_sex <- DESeq(dds_sex)

# Get DESeq2 results and filter for significantly differentially expressed genes
res_sex <- results(dds_sex, contrast = c("group", "M", "F"))
sig_genes_sex <- rownames(res_sex)[which(res_sex$padj < 0.05 & abs(res_sex$log2FoldChange) > 1)]
cat(sprintf("Sex data: Significantly DE genes (padj < 0.05, |log2FC| > 1): %d genes\n", length(sig_genes_sex)))

# Filter to keep only significant genes
dds_sex <- dds_sex[sig_genes_sex, ]

# Get rlog transformed counts (similar to your example)
rlogcounts_sex <- assay(rlog(dds_sex, blind = FALSE))

# Prepare data for heatmap - convert to long format
sex_heatmap_data <- as.data.frame(rlogcounts_sex) %>%
  rownames_to_column(var = "ensembl") %>%
  left_join(gene_mapping, by = "ensembl") %>%
  mutate(GeneName = ifelse(is.na(gene_name), ensembl, gene_name)) %>%
  dplyr::select(ensembl, GeneName, everything(), -gene_name) %>%
  pivot_longer(cols = -c(ensembl, GeneName),
               names_to = "sample",
               values_to = "expression")

# Reshape to wide format for clustering
heatmap_sex_wide <- sex_heatmap_data %>%
  pivot_wider(names_from = sample, values_from = expression)

# Apply k-means clustering and extract gene order
cluster_result_sex <- get_kmeans_cluster_order(heatmap_sex_wide, 
                                               group_info = col_data$group,
                                               num_clusters = 3)

# Print cluster sizes
print(cluster_result_sex$cluster_sizes)

# Add gene order and cluster info as factors to the long data format
heatmap_sex_long <- sex_heatmap_data %>%
  mutate(GeneName = factor(GeneName, levels = cluster_result_sex$order),
         Cluster = factor(cluster_result_sex$clusters[match(GeneName, names(cluster_result_sex$clusters))],
                          levels = c(1, 2, 3),
                          labels = paste0(1:3, "\nn = ", as.numeric(cluster_result_sex$cluster_sizes))))

# Create CSV with cluster assignments for sex
sex_cluster_csv <- data.frame(
  ensembl = heatmap_sex_wide$ensembl[match(names(cluster_result_sex$clusters), heatmap_sex_wide$GeneName)],
  GeneName = names(cluster_result_sex$clusters),
  Cluster = cluster_result_sex$clusters,
  cluster_result_sex$group_means
)
write_csv(sex_cluster_csv, "Sig_Kmeans_sex_gene_clusters.csv")
print("Sex cluster assignments saved to: Sig_Kmeans_sex_gene_clusters.csv")

# Calculate z-scores for each gene across samples
heatmap_sex_long <- heatmap_sex_long %>%
  group_by(ensembl) %>%
  mutate(zscore = scale(expression)[,1]) %>%
  ungroup()

# Define color palette (blue to red diverging scale)
color_palette <- colorRampPalette(c("#0571b0", "#92c5de", "#f7f7f7", "#fdb863", "#ca0020"))(100)

# Add Sex info
heatmap_sex_long <- heatmap_sex_long %>%
  mutate(Sex = factor(ifelse(grepl("^F", sample), "F", "M"), levels = c("F", "M")))

#      and strip.placement = "outside" so they sit cleanly below the tiles.
sex_heatmap_plot <- ggplot(heatmap_sex_long, aes(x = sample, y = GeneName, fill = zscore)) +
  geom_tile(color = NA) +
  scale_fill_gradientn(colours = color_palette, 
                       na.value = "grey90",
                       limits = c(-4, 4),
                       oob = scales::squish,
                       name = "Z-score") +
  facet_grid(Cluster ~ Sex,
             scales = "free",
             space = "free",
             switch = "x") +           # moves Sex strip labels to the bottom
  labs(title = "Gene Expression - Sex Comparison (K-means Clustered)",
       x = "",
       y = "") +
  theme_minimal() +
  theme(axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        panel.spacing = unit(0.5, "lines"),
        strip.placement = "outside",             # places strip outside the axis area
        strip.text.x = element_text(size = 14, face = "bold"),
        strip.text.y = element_text(size = 12, face = "bold", angle = 0),
        strip.background = element_rect(fill = "grey90", color = NA),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
        legend.position = "top",
        legend.key.width = unit(2, "cm"),
        legend.key.height = unit(0.5, "cm"))

ggsave("Sig_Kmeans_heatmap_all_genes_sex.pdf", plot = sex_heatmap_plot, width = 8, height = 14, units = "in")
print(sex_heatmap_plot)

###########################################################################
########################### DOMAIN DATA ###################################
###########################################################################

# Load and prepare domain data - ALL genes
domain_data <- read.csv("../Domain_expression.csv", header = TRUE)
domain_counts_mat <- domain_data %>% dplyr::select(-ensembl) %>% as.matrix()
rownames(domain_counts_mat) <- domain_data$ensembl
domain_sample_names <- colnames(domain_counts_mat)

# FIX: Added diagnostic warning if any samples don't match AT/M/PT patterns,
#      so silently dropped samples are visible in the console output.
# NOTE: Sample names are assumed to encode both sex and domain, e.g. "F_AT1", "M_PT2".
#       Sex is extracted from the leading F/M prefix; domain from the AT/M/PT suffix.
domain_col_data <- data.frame(
  sample = domain_sample_names,
  group = case_when(
    grepl("AT$", domain_sample_names) ~ "AT",
    grepl("M$",  domain_sample_names) ~ "M",
    grepl("PT$", domain_sample_names) ~ "PT"
  ),
  sex = ifelse(grepl("^F", domain_sample_names), "F", "M")
)
unmatched <- domain_col_data$sample[is.na(domain_col_data$group)]
if (length(unmatched) > 0) {
  warning(sprintf("Domain data: %d sample(s) did not match AT/M/PT and will be dropped: %s",
                  length(unmatched), paste(unmatched, collapse = ", ")))
}
domain_col_data <- domain_col_data %>% filter(!is.na(group))
rownames(domain_col_data) <- domain_col_data$sample
domain_counts_mat <- domain_counts_mat[, domain_col_data$sample]
domain_col_data$group <- factor(domain_col_data$group, levels = c("AT", "M", "PT"))
domain_col_data$sex   <- factor(domain_col_data$sex,   levels = c("F", "M"))

# Filter genes: keep only genes with mean count >10 across samples
gene_means_domain <- rowMeans(domain_counts_mat)
keep_genes_domain <- gene_means_domain > 10
cat(sprintf("Domain data: Before filtering: %d genes\n", nrow(domain_counts_mat)))
cat(sprintf("Domain data: Genes with mean count >10: %d genes\n", sum(keep_genes_domain)))
cat(sprintf("Domain data: Filtering out %d genes\n", sum(!keep_genes_domain)))
domain_counts_mat <- domain_counts_mat[keep_genes_domain, ]
cat(sprintf("Domain data: After filtering: %d genes retained\n", nrow(domain_counts_mat)))

# DESeq2 analysis for domains - use FILTERED dataset
# Design includes sex as a covariate so domain DE results are sex-corrected,
# and sex effects are properly partitioned rather than absorbed into residuals.
dds_domain <- DESeqDataSetFromMatrix(countData = domain_counts_mat, colData = domain_col_data, design = ~ sex + group)
dds_domain <- DESeq(dds_domain)

# Get DESeq2 results and filter for significantly differentially expressed genes
# For multi-group comparison, get all pairwise comparisons
res_domain_AT_vs_M  <- results(dds_domain, contrast = c("group", "AT", "M"))
res_domain_AT_vs_PT <- results(dds_domain, contrast = c("group", "AT", "PT"))
res_domain_M_vs_PT  <- results(dds_domain, contrast = c("group", "M",  "PT"))

# Get genes significant in ANY comparison
sig_genes_AT_M  <- rownames(res_domain_AT_vs_M)[which(res_domain_AT_vs_M$padj  < 0.05 & abs(res_domain_AT_vs_M$log2FoldChange)  > 1)]
sig_genes_AT_PT <- rownames(res_domain_AT_vs_PT)[which(res_domain_AT_vs_PT$padj < 0.05 & abs(res_domain_AT_vs_PT$log2FoldChange) > 1)]
sig_genes_M_PT  <- rownames(res_domain_M_vs_PT)[which(res_domain_M_vs_PT$padj  < 0.05 & abs(res_domain_M_vs_PT$log2FoldChange)  > 1)]

sig_genes_domain <- unique(c(sig_genes_AT_M, sig_genes_AT_PT, sig_genes_M_PT))
cat(sprintf("Domain data: Significantly DE genes in any comparison (padj < 0.05, |log2FC| > 1): %d genes\n", length(sig_genes_domain)))

# Filter to keep only significant genes
dds_domain <- dds_domain[sig_genes_domain, ]

# Get rlog transformed counts
rlogcounts_domain <- assay(rlog(dds_domain, blind = FALSE))

# Prepare data for heatmap - convert to long format
domain_heatmap_data <- as.data.frame(rlogcounts_domain) %>%
  rownames_to_column(var = "ensembl") %>%
  left_join(gene_mapping, by = "ensembl") %>%
  mutate(GeneName = ifelse(is.na(gene_name), ensembl, gene_name)) %>%
  dplyr::select(ensembl, GeneName, everything(), -gene_name) %>%
  pivot_longer(cols = -c(ensembl, GeneName),
               names_to = "sample",
               values_to = "expression")

# Reshape to wide format for clustering
heatmap_domain_wide <- domain_heatmap_data %>%
  pivot_wider(names_from = sample, values_from = expression)

# Apply k-means clustering and extract gene order
cluster_result_domain <- get_kmeans_cluster_order(heatmap_domain_wide, 
                                                  group_info = domain_col_data$group,
                                                  num_clusters = 3)

# Print cluster sizes
print("Domain clustering results:")
print(cluster_result_domain$cluster_sizes)

# Add gene order and cluster info as factors to the long data format
heatmap_domain_long <- domain_heatmap_data %>%
  mutate(GeneName = factor(GeneName, levels = cluster_result_domain$order),
         Cluster = factor(cluster_result_domain$clusters[match(GeneName, names(cluster_result_domain$clusters))],
                          levels = c(1, 2, 3),
                          labels = paste0(1:3, "\nn = ", as.numeric(cluster_result_domain$cluster_sizes))))

# Create CSV with cluster assignments for domain
domain_cluster_csv <- data.frame(
  ensembl = heatmap_domain_wide$ensembl[match(names(cluster_result_domain$clusters), heatmap_domain_wide$GeneName)],
  GeneName = names(cluster_result_domain$clusters),
  Cluster = cluster_result_domain$clusters,
  cluster_result_domain$group_means
)
write_csv(domain_cluster_csv, "Sig_Kmeans_domain_gene_clusters.csv")

# Calculate z-scores for each gene across samples
heatmap_domain_long <- heatmap_domain_long %>%
  group_by(ensembl) %>%
  mutate(zscore = scale(expression)[,1]) %>%
  ungroup()

# Add Domain and Sex info for nested strip labelling
heatmap_domain_long <- heatmap_domain_long %>%
  mutate(
    Domain = factor(case_when(
      grepl("AT$", sample) ~ "AT",
      grepl("M$",  sample) ~ "M",
      grepl("PT$", sample) ~ "PT"
    ), levels = c("AT", "M", "PT")),
    Sex = factor(ifelse(grepl("^F", sample), "F", "M"), levels = c("F", "M"))
  )

# Create heatmap with nested x-axis strips: F/M (inner) nested within Domain (outer).
# facet_nested() from ggh4x handles the two-level column strip automatically.
# switch = "x" moves both strip levels to the bottom, matching the reference figure layout.
domain_heatmap_plot <- ggplot(heatmap_domain_long, aes(x = sample, y = GeneName, fill = zscore)) +
  geom_tile(color = NA) +
  scale_fill_gradientn(colours = color_palette,
                       na.value = "grey90",
                       limits = c(-4, 4),
                       oob = scales::squish,
                       name = "Z-score") +
  facet_nested(Cluster ~ Domain + Sex,
               scales = "free",
               space = "free",
               switch = "x",                  # both strip levels move to the bottom
               nest_line = element_line(linewidth = 0.2, color = "grey60")) +
  labs(title = "Gene Expression - Domain Comparison (K-means Clustered)",
       x = "",
       y = "") +
  theme_minimal() +
  theme(axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        panel.spacing = unit(0.1, "lines"),
        strip.placement = "outside",           # strips sit below the tile area
        # Inner strip (F / M)
        strip.text.x.bottom = element_text(size = 12, face = "bold"),
        # Outer strip (AT / M-domain / PT) - slightly larger to match reference
        strip.text.x = element_text(size = 14, face = "bold"),
        strip.text.y = element_text(size = 12, face = "bold", angle = 0),
        strip.background = element_rect(fill = "grey90", color = NA),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
        legend.position = "top",
        legend.key.width = unit(2, "cm"),
        legend.key.height = unit(0.5, "cm"))

ggsave("Sig_Kmeans_heatmap_all_genes_domain.pdf", plot = domain_heatmap_plot, width = 14, height = 14, units = "in")
print(domain_heatmap_plot)
