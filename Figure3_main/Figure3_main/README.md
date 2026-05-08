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
- lncRNA dataset from adult anal fins, gonads and muscles: `longRNA_Raw.txt`
- Small RNA dataset from anal fins, gonads and muscles: `SmallRNA_Raw.txt`

**Scripts**
- Use `deseq_fromRawcounts.R` to plot Violin plots, Line plots and correlations for the Genes of interest.

---

### Domain Dataset

**Data**
- `Domain_expression.csv`

**Scripts**
- Use `deseq_fromRawcounts.R`

---

### Single Cell Dataset

**Data**
- `cds_object.rds`

**Scripts**
- Use `1to1_ortho_Umap_zebra.R` to plot the pigment cell expression for each gene orthologous to our gene of interest.
