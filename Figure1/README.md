## Figure 2 Panels - Masoko Eggspot Paper

The following sections correspond to all the panels in Figure 2:

**Panels:**
- **A.** [GWAS](#a-gwas-for-eggspot-theis-count--eggspot-hue) for egg spot number (Theis count)  
- **B.** Zoom plot for lncRNA (chr22: ENSACLG00000029630)  
- **C.** Zoom plot for (chr23: oca2)  
- **D.** Genotype-Phenotype violin plots (chr22:12607876)  
- **E.** Genotype-Phenotype violin plots (chr23:12211546)  
- **F.** Potential binding sites for miR10 miRNA  
- **G.** [GWAS](#a-gwas-for-eggspot-theis-count--eggspot-hue) for egg spot Hue  
- **H.** Zoom plot for (chr13: tmem254-RHOQ)  
- **I.** Zoom plot for (chr23: aldh9a1-tmco1)  
- **K.** Genotype-Phenotype violin plots (chr13:14886065)  

---

### A & G. GWAS for Eggspot (Theis Count & Eggspot Hue) [example - eggspot number]

**Perform GWAS using GEMMA:**  
Follow these steps to perform GWAS using GEMMA (assumes QC steps are already done):

#### i. Input Files
- **samples_for_theis.txt** – includes all samples used for the analyses (no header).  
- **phenotype.txt** – formatted as:

```
cichl187579993  cichl187579993  3.5
cichl187579994  cichl187579994  6
cichl187579995  cichl187579995  6.5
cichl187579996  cichl187579996  4.5
cichl187579997  cichl187579997  4.5
cichl187579999  cichl187579999  6.5
cichl187580000  cichl187580000  5.5
cichl187580001  cichl187580001  8
cichl187580002  cichl187580002  3.5
```

#### ii. Convert VCF to BIM/BAM
**Before converting, make sure GEMMA is installed:**
```bash
conda create -n GEMMA -c bioconda gemma -y
```
**Run conversion:**
```bash
bash VCF_to_BIM_BAM.bash &
```
This outputs `.bed`, `.fam`, etc. Check that the last column in `.fam` (phenotype) is read correctly (0's are read as NA `-9` in GEMMA).

#### iii. Compute Relatedness / Kinship Matrix
```bash
sbatch compute_relatedness.job
```
Output: `.cXX.txt` file in your output folder (kinship matrix).

#### iv. Perform GWAS
```bash
sbatch perform_GWAS.job
```
Covariate (e.g., standard length) should be formatted as:

```
1  88.8
1  88.49
1  78.5
1  98.24
1  58.57
1  53.7
1  86.07
1  96.01
1  87.3
```

#### v. Plot Data & Extract Top SNPs
```R
Rscript Plot.R
```

---

### B. Zoom Plot for lncRNA (chr22: ENSACLG00000029630)
```R
Rscript OCA2_Zoomplot.R
```

### C. Zoom Plot for chr23: oca2
```R
Rscript lncRNA_Zoomplot.R
```

### H & I. Zoom Plots
- Same as B & C but update `"GWAS"`, `"output"`, `"chrom"`, `"start"`, `"end"` in the script.

---

### D, E & K. Genotype-Phenotype Plots
- Requires a GT file:
```bash
bash extract_GT.bash &
```
- Generate plots:
```R
Rscript GenoXPheno.R
```

---

### F. Potential Binding Sites for miR10 miRNA

