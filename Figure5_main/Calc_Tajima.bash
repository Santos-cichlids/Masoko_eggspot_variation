#!/bin/bash
for i in {1..20} {22..23};
 do sbatch calculate_Pi_PER_SNP_TAJIMA.job $i;
done
