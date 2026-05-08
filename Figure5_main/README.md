# Figure 5 main

This folder contains scripts used to generate the analyses and plots for the corresponding figure in the manuscript.

---

#### removing sites that are missing in more that 30% samples
- a loop can was run over all the chr files 
for i in {1..20};
 do sbatch 1_extract_failed_sites.job $i;
done

for i in {22..23};
 do sbatch 1_extract_failed_sites.job $i;
done
- this gives you all the failed sites per population, then use a similar loop to remove those sites
for i in {1..20};
 do sbatch 2_removeFailedSites.job $i;
done

for i in {22..23};
 do sbatch 2_removeFailedSites.job $i;
done

---

## xpEHH Calculation

### `Calc_bootstrap_xpEHH.R`
Use this script **only when the two populations have unequal sample sizes**.

- Iteratively runs the standard xpEHH analysis on the larger population  
- In each iteration, a random subset of individuals (equal in size to the smaller population) is sampled  
- This approach reduces bias introduced by sample size differences

### `Calc_xpEHH_IHH.R`
Use this script when **both populations have roughly equal sample sizes**.

- Implements the standard xpEHH / IHH calculation workflow  
- No bootstrapping or subsampling is performed

---

## Threshold Calculation

### `Calculate_thresholds.R`
Standalone script used to calculate statistical thresholds for downstream analyses.

---

## Diversity Statistics

### `Calc_Tajima.bash`
Wrapper script to run `Calc_Pi_Tajima.job` across all chromosome files in parallel.

---

## Plotting

### `plot_xpEHH_IHS_Taji_GWAS.R`
Generates genome-wide plots combining:

- xpEHH  
- iHS  
- Tajima’s D  
- GWAS results

---
## Lake Malawi-wide phylogeny

### Plot Allele frequencies



