#!/bin/bash
#SBATCH --job-name=LD_blocks
#SBATCH -A XXXXXXXXX
#SBATCH --array=1-23
#SBATCH --output=logs/LD_blocks_%A_%a.out
#SBATCH --error=logs/LD_blocks_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=03:00:00
#SBATCH --partition=XXXXXXXX

set -euo pipefail

module load plink-1.9-gcc-5.4.0-sm3ojoi
mkdir -p logs ld_results/blocks

CHR=${SLURM_ARRAY_TASK_ID}

# chr21 does not exist — skip
if [[ "${CHR}" -eq 21 ]]; then exit 0; fi

VCF="../masoko_mono_biall.GQ35.depth.filtered.chr${CHR}.vcf.gz"
OUT="ld_results/blocks/cichlid_chr${CHR}"

if [[ ! -f "${VCF}" ]]; then
    echo "ERROR: VCF not found: ${VCF}" >&2
    exit 1
fi

# Step 1: QC filter and make BED
plink \
    --vcf "${VCF}" \
    --allow-extra-chr \
    --chr "${CHR}" \
    --double-id \
    --set-missing-var-ids @:# \
    --maf 0.05 \
    --geno 0.1 \
    --mind 0.5 \
    --make-bed \
    --threads "${SLURM_CPUS_PER_TASK}" \
    --out "${OUT}_qc"

# Step 2: Define LD blocks
awk '{print $1, $2, 1}' "${OUT}_qc.fam" > /tmp/dummy_pheno_chr${CHR}.txt

plink \
    --bfile "${OUT}_qc" \
    --allow-extra-chr \
    --allow-no-sex \
    --blocks \
    --blocks-max-kb 100 \
    --pheno /tmp/dummy_pheno_chr${CHR}.txt \
    --1 \
    --threads "${SLURM_CPUS_PER_TASK}" \
    --out "${OUT}"

rm -f /tmp/dummy_pheno_chr${CHR}.txt
rm -f "${OUT}_qc".{bed,bim,fam,log,nosex}