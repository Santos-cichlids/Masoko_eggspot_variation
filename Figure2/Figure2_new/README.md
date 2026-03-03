# Figure 2

## Pipeline Overview

### 1. QC
### 2. Trimming
### 3. Annotation and Counts
### 4. Plotting

---

## Datasets & Scripts

### Adult Dataset

**Data**
- lncRNA dataset from adult anal fins, gonads and muscles: `TPM_counts_lncRNA.txt`
- Small RNA dataset from anal fins, gonads and muscles: `Masoko_NormCounts_smallRNA.txt`

**Scripts**
- Use `Plot_violin_corr_line.R` to plot Violin plots, Line plots and correlations for the Genes of interest.

---

### Domain & Sex Dataset

**Data**
- `ACAL_gene_counts_Extracted.csv`
- `Domain_expression.csv`

**Scripts**
- Use `line_plot_domain_divergence.R` to plot the line plots for genes of interest.
- Use `heatmap_domain_sex.R` to plot heatmap for domain and sex at the time of trait initiation and divergence.

---

### Single Cell Dataset

**Data**
- `cds_object.rds`

**Scripts**
- Use `1to1_ortho_Umap_zebra.R` to plot the pigment cell expression for each gene orthologous to our gene of interest.
