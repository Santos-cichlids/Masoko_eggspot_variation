# Load required packages
library(DESeq2)
library(tidyverse)
library(ggrepel)

# Gene mapping - Set 1
gene_mapping <- data.frame(
  ensembl = c("ENSACLG00000017026","ENSACLG00000026836", "ENSACLG00000026955",
              "ENSACLG00000024190", "ENSACLG00000025445", "ENSACLG00000016986"),
  gene_name = c("tmem254", "aldh9a1a.a", "tmco1", "oca2", "plxdc1", "rhoq")
)

# Load and prepare sex data - keep ALL genes for proper normalization
count_data <- read.csv("../ACAL_gene_counts_Extracted.csv", header = TRUE)
counts_mat <- count_data %>% select(-ensembl) %>% as.matrix()
rownames(counts_mat) <- count_data$ensembl
sample_names <- colnames(counts_mat)
col_data <- data.frame(sample = sample_names, group = ifelse(grepl("^F", sample_names), "F", "M")) %>%
  filter(group %in% c("F", "M"))
rownames(col_data) <- col_data$sample
counts_mat <- counts_mat[, col_data$sample]
col_data$group <- factor(col_data$group, levels = c("F", "M"))

# DESeq2 analysis for sex - use FULL dataset for proper normalization
dds_sex <- DESeqDataSetFromMatrix(countData = counts_mat, colData = col_data, design = ~ group)
dds_sex <- DESeq(dds_sex)
normalized_counts_sex <- counts(dds_sex, normalized = TRUE)

# Now subset AFTER normalization
sex_norm_data <- normalized_counts_sex %>%
  as.data.frame() %>%
  rownames_to_column("ensembl") %>%
  filter(ensembl %in% gene_mapping$ensembl) %>%  # Subset here
  pivot_longer(-ensembl, names_to = "sample", values_to = "normalized_count") %>%
  left_join(col_data, by = "sample") %>%
  left_join(gene_mapping, by = "ensembl") %>%
  mutate(log_norm_counts = log2(normalized_count + 1))

# Calculate mean expression for sex
sex_mean_data <- sex_norm_data %>%
  group_by(gene_name, group) %>%
  summarise(mean_log_counts = mean(log_norm_counts, na.rm = TRUE), .groups = "drop") %>%
  mutate(group_num = as.numeric(group))

# Load and prepare domain data - keep ALL genes for proper normalization
domain_data <- read.csv("../Domain_expression.csv", header = TRUE)  # Don't filter yet!
domain_counts_mat <- domain_data %>% select(-ensembl) %>% as.matrix()
rownames(domain_counts_mat) <- domain_data$ensembl
domain_sample_names <- colnames(domain_counts_mat)
domain_col_data <- data.frame(
  sample = domain_sample_names,
  group = case_when(
    grepl("AT$", domain_sample_names) ~ "AT",
    grepl("M$", domain_sample_names) ~ "M",
    grepl("PT$", domain_sample_names) ~ "PT")
) %>% filter(!is.na(group))
rownames(domain_col_data) <- domain_col_data$sample
domain_counts_mat <- domain_counts_mat[, domain_col_data$sample]
domain_col_data$group <- factor(domain_col_data$group, levels = c("AT", "M", "PT"))

# DESeq2 analysis for domains - use FULL dataset for proper normalization
dds_domain <- DESeqDataSetFromMatrix(countData = domain_counts_mat, colData = domain_col_data, design = ~ group)
dds_domain <- DESeq(dds_domain)
normalized_counts_domain <- counts(dds_domain, normalized = TRUE)

# Now subset AFTER normalization
domain_norm_data <- normalized_counts_domain %>%
  as.data.frame() %>%
  rownames_to_column("ensembl") %>%
  filter(ensembl %in% gene_mapping$ensembl) %>%  # Subset here
  pivot_longer(-ensembl, names_to = "sample", values_to = "normalized_count") %>%
  left_join(domain_col_data, by = "sample") %>%
  left_join(gene_mapping, by = "ensembl") %>%
  mutate(log_norm_counts = log2(normalized_count + 1))

# Calculate mean expression for domains
domain_mean_data <- domain_norm_data %>%
  group_by(gene_name, group) %>%
  summarise(mean_log_counts = mean(log_norm_counts, na.rm = TRUE), .groups = "drop") %>%
  mutate(group_num = as.numeric(group))

# Get labels for the end points only
sex_labels <- sex_mean_data %>%
  group_by(gene_name) %>%
  filter(group == "M") %>%
  ungroup()

# Create sex plot
sex_plot <- ggplot(sex_mean_data, aes(x = group, y = mean_log_counts, group = gene_name, color = gene_name)) +
  geom_line(size = 1.2) +
  geom_point(size = 3.5) +
  geom_text_repel(data = sex_labels,
                  aes(label = gene_name), 
                  size = 4,
                  nudge_x = 0.3,
                  hjust = 0,
                  direction = "y",
                  segment.size = 0.3,
                  segment.color = "grey50") + scale_color_manual(values = c(
    "tmem254" = "#1B9E77",
    "aldh9a1a.a" = "#D95F02",
    "tmco1" = "#7570B3",
    "oca2" = "#E7298A",
    "plxdc1" = "#66A61E",
    "rhoq" = "#E6AB02"  
  )) +
labs(x = "Sex", y = "Mean Log2 Normalized Counts", 
       title = "Gene Expression Across Sex") +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none",
    axis.text = element_text(color = "black"),
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.title = element_text(face = "bold")
  ) +
  coord_cartesian(clip = "off") +
  theme(plot.margin = margin(5, 50, 5, 5))

# Get labels for the end points only
domain_labels <- domain_mean_data %>%
  group_by(gene_name) %>%
  filter(group == "PT") %>%
  ungroup()

# Create domain plot
domain_plot <- ggplot(domain_mean_data, aes(x = group, y = mean_log_counts, group = gene_name, color = gene_name)) +
  geom_line(size = 1.2) +
  geom_point(size = 3.5) +
  geom_text_repel(data = domain_labels,
                  aes(label = gene_name), 
                  size = 4,
                  nudge_x = 0.3,
                  hjust = 0,
                  direction = "y",
                  segment.size = 0.3,
                  segment.color = "grey50") +
  scale_color_manual(values = c(
    "tmem254" = "#1B9E77",
    "aldh9a1a.a" = "#D95F02",
    "tmco1" = "#7570B3",
    "oca2" = "#E7298A",
    "plxdc1" = "#66A61E",
    "rhoq" = "#E6AB02"
  ))+labs(x = "Domain", y = "Mean (Log Normalized Counts)", 
       title = "Gene Expression Across Domains") +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none",
    axis.text = element_text(color = "black"),
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.title = element_text(face = "bold")
  ) +
  coord_cartesian(clip = "off") +
  theme(plot.margin = margin(5, 50, 5, 5))

# Save plots
ggsave("sex_mean_expression_all_genes.pdf", plot = sex_plot, width = 8, height = 6, units = "in")
ggsave("domain_mean_expression_all_genes.pdf", plot = domain_plot, width = 10, height = 8, units = "in")

# Display plots
print(sex_plot)
print(domain_plot)
