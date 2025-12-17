#!/bin/bash

# === Step 1: Extract genotype calls from VCF ===
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' ../Masoko_SUBSET_FILTERED_for_THEIS.vcf.gz > genotype_data_for_plotting.txt

# === Step 2: Extract chr and ps columns from candidate loci file ===
awk 'NR==1 {print "chr\tps"} NR>1 {print $1"\t"$3}' ../Only_SL_as_Covar/top_candidate_loci.txt > chr_ps_extracted.txt

# === Step 3: Subset genotype data for those loci ===
awk 'FNR==NR && NR>1 {keys["chr"$1"\t"$2]; next} ($1"\t"$2) in keys' chr_ps_extracted.txt genotype_data_for_plotting.txt > data_for_plotting_THEIS_SL.txt

