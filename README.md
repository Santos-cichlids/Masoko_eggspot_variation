# Masoko Dataset Analyses

This folder contains all scripts used to analyze the **Masoko** dataset.

## Analyses Overview

1. [GWAS](#1-gwas)  
   - `egg_spot_Area`  
   - `egg_spot_size`  
   - `egg_spot_number`  
   - `egg_spot_number_Theis`  
   - `egg_spot_Hue`  

2. [Selective Scans](#2-selective-scans)  
   - FST  
   - DXY  
   - Tajima's D  
   - Pi  
   - xp-EHH  

3. [Expression Analyses](#3-expression-analyses)  
   - RNAseq data from wild-caught individuals  
   - RNAseq data from Aaron's experiment
      
4. [3'UTR and mi-RNA binding](#4-3utr-and-mi-rna-binding) 
   - Extracting 3'UTR from vcf 
---

## 1. GWAS  

_Metadata: [Google Sheet](https://docs.google.com/spreadsheets/d/1rJda7cw5H4N52hMqD5NrAQfIXbaVSL6f/edit?gid=1587235546#gid=1587235546)_

All GWAS were carried out using **GEMMA** and **ANGSD**.  
📎 [Pipeline Reference](https://github.com/Santos-cichlids/Methods-and-tools/tree/main/Bioinformatics/GWAS)

### a. `egg_spot_Area`
- **N** = 288 males  
- **Phenotype column**: `Aligned_Total_Eggspot_area_PropOfFin`  
- **Directory**:  
  `/rds/project/rds-8b3VcZwY7rY/projects/cichlid/suhaib/masoko/GWAS/masoko_eggspot_Area`

**Filtering Criteria:**
```r
subset(masoko, 
       Sex == "Male" & 
       seq_depth > 10 & 
       Aligned_SpotSize_PropOfFin != "NA" & 
       SL != "NA")
```

### b. `egg_spot_size`
- **N** = 288 males  
- **Phenotype column**: `Aligned_SpotSize_PropOfFin`  
- **Directory**:  
  `/rds/project/rds-8b3VcZwY7rY/projects/cichlid/suhaib/masoko/GWAS/masoko_eggspot_size`

**Filtering Criteria:**
```r
subset(masoko, 
       Sex == "Male" & 
       seq_depth > 10 & 
       Aligned_SpotSize_PropOfFin != "NA" & 
       SL != "NA")
```

### c. `egg_spot_number`
- **N** = 260 males  
- **Phenotype column**: `EggspotNumber`  
- **Directory**:  
  `/rds/project/rds-8b3VcZwY7rY/projects/cichlid/suhaib/masoko/GWAS/masoko_GWAS_eggsot_number`

**Filtering Criteria:**
```r
subset(masoko, 
       Sex == "Male" & 
       seq_depth > 10 & 
       EggspotNumber != "NA" & 
       SL != "NA")
```

### d. `egg_spot_number_Theis`
- **N** = 243 males  
- **Phenotype column**: `EggspotNumberTheis`  
- **Directory**:  
  `/rds/project/rds-8b3VcZwY7rY/projects/cichlid/suhaib/masoko/GWAS/masoko_GWAS_eggSpot_theis`

**Filtering Criteria:**
```r
subset(masoko, 
       Sex == "Male" & 
       seq_depth > 10 & 
       EggspotNumberTheis != "NA" & 
       SL != "NA")
```

### e. `egg_spot_Hue`
- **N** = 159 males  
- **Phenotype column**: `spot.hue`  
- **Directory**:  
  `/rds/project/rds-8b3VcZwY7rY/projects/cichlid/suhaib/masoko/GWAS/masoko_HUE`

**Filtering Criteria:**
```r
subset(masoko, 
       Sex == "Male" & 
       seq_depth > 10 & 
       spot.hue != "NA" & 
       SL != "NA")
```

---

## 2. Selective Scans

### Filtering:

- VCFs split by ecology: **Benthic**, **Littoral**, and **Middle**
- Sites with >30% missingness per group removed

### Directories:

- VCFs with monomorphic sites:  
  `/rds/project/rds-8b3VcZwY7rY/projects/cichlid/suhaib/masoko/genome_scans_ALL/BENTHIC_LITTORAL/`

- VCFs with biallelic sites only:  
  `/rds/project/rds-8b3VcZwY7rY/projects/cichlid/suhaib/masoko/genome_scans/`

### a. FST  
- **Directory**:  
  `/rds/project/rds-8b3VcZwY7rY/projects/cichlid/suhaib/masoko/genome_scans/FST_BIALLELICS/`

### b. DXY  
- Uses VCF with monomorphic sites included.

### c. Tajima’s D  
- **Directory**:  
  `/rds/project/rds-8b3VcZwY7rY/projects/cichlid/suhaib/masoko/genome_scans_ALL/BENTHIC_LITTORAL/PER_SNP`

- Calculated using `VCFTOOLS` with **10kb non-overlapping windows**

**Run command:**
```bash
bash Calc_Pi_Tajima.bash
# Runs Calc_Pi_Tajima.job over all chromosomes
```

### d. Pi  
- **Directory**:  
  `/rds/project/rds-8b3VcZwY7rY/projects/cichlid/suhaib/masoko/genome_scans_ALL/BENTHIC_LITTORAL/PER_SNP`

- Calculated with `VCFTOOLS` using **10kb window** and **5kb step**

**Run command:**
```bash
bash Calc_Pi_Tajima.bash
# Runs Calc_Pi_Tajima.job over all chromosomes
```

### e. xp-EHH  
> _Requires phased data (phased to fAstCal1.2 reference)_

- 98 **benthic** + 187 **littoral** individuals  
- Script: `Calc_IHH_xpEHH.R`  
- Outputs: IHH (population-specific) + xp-EHH

- **Directory**:  
  `/rds/project/rds-8b3VcZwY7rY/projects/cichlid/suhaib/masoko/genome_scans_ALL/BENTHIC_LITTORAL/EHH_IHH_Selection/biallelics/HIGH_SEQ_depth/`

---

## 3. Expression Analyses

### a. RNAseq from 58 wild-caught individuals  
- Analyzed using: `Plot_Xpr_by_eco_top_genes.R`

### b. RNAseq from Aaron’s experiment  
- Analyzed using: `Domain_Sex_deseq_Aaron.R`
## 4. 3'UTR and mi-RNA binding
- Extracting 3'UTR from vcf...
