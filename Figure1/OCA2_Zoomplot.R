#!/usr/bin/env Rscript

# Load libraries (explicit for namespace conflicts)
library(rtracklayer)
library(dplyr)

# Input files
GWAS <- "/rds/project/rds-8b3VcZwY7rY/projects/cichlid/suhaib/masoko/GWAS/masoko_GWAS_eggSpot_theis/Only_SL_as_Covar/output/Masoko_GWAS_with_covariates_SL_covar.assoc.txt"
gtf_file <- "/rds/project/rds-8b3VcZwY7rY/projects/cichlid/Bettina/data/Suhaib_expression_20250422/fAstCal1.2.111.gtf.gz"
output <- "OCA2_"

# Region of interest
chrom <- "chr23"
start <- 12149100
end <- 12315000

# Thresholds
total_snps <- 1000000
bonf <- -log10(0.05 / total_snps)
suggestive <- -log10(1e-5)
 

# Read GWAS
gwas <- read.table(GWAS, header = TRUE)
gwas$chr <- paste0("chr", gwas$chr)
gwas$neglogp <- -log10(gwas$p_score)
gwas_zoom <- dplyr::filter(gwas, chr == chrom, ps >= start, ps <= end)

# Read GTF
gtf <- as.data.frame(rtracklayer::import(gtf_file))
gtf_sub <- dplyr::filter(gtf, seqnames == sub("chr", "", chrom), start >= start, end <= end)

# Get unique genes and assign levels to avoid overlap
genes <- dplyr::filter(gtf_sub, type == "gene") %>%
  dplyr::select(gene_id, gene_name, start, end, strand) %>%
  dplyr::distinct() %>%
  dplyr::arrange(start)

# Assign plotting levels to genes to avoid overlap
margin <- 5000
genes$level <- 1
genes$t_start <- genes$start
genes$t_end <- genes$end

for (i in 1:nrow(genes)) {
  tx <- dplyr::filter(gtf_sub, gene_id == genes$gene_id[i], type == "transcript")
  if (nrow(tx) > 0) {
    genes$t_start[i] <- min(tx$start)
    genes$t_end[i] <- max(tx$end)
  }
  for (lvl in 1:3) {
    overlap <- any(genes$level[1:(i - 1)] == lvl &
                     genes$t_start[i] <= genes$t_end[1:(i - 1)] + margin &
                     genes$t_end[i] >= genes$t_start[1:(i - 1)] - margin)
    if (!overlap) {
      genes$level[i] <- lvl
      break
    }
  }
}

# Plot
max_lvl <- max(genes$level)
ann_height <- 0.4 + 0.6 * max_lvl
pdf(paste0(output, "zoom_plot_with_thresh.pdf"), width = 8 + ann_height, height = 8 + ann_height)

layout(matrix(1:2, ncol = 1), heights = c(4, ann_height))
par(mar = c(0.5, 5, 3, 2), oma = c(4, 0, 0, 0))

# Plot GWAS
plot(gwas_zoom$ps, gwas_zoom$neglogp, type = "n", xlab = "", ylab = "-log10(p-value)", 
     xlim = c(start, end), ylim = range(c(gwas_zoom$neglogp, bonf)), las = 1,
     xaxt = "n", main = "")

abline(h = bonf, col = "red", lty = 1)
abline(h = suggestive, col = "blue", lty = 2)

#colz <- ifelse(gwas_zoom$neglogp >= bonf, "red",
               #ifelse(gwas_zoom$neglogp >= suggestive, "blue", "black"))
colz <- "black"
cexz <- 0.6
points(gwas_zoom$ps, gwas_zoom$neglogp, col = colz, pch = 16, cex = cexz)

# Plot genes
par(mar = c(4, 5, 0.5, 2))
plot(0, type = "n", xlab = "Genomic Position (bp)", ylab = "", xlim = c(start, end), ylim = c(0, 1.1), axes = FALSE)
axis(1, at = pretty(c(start, end)), labels = format(pretty(c(start, end)), big.mark = ","))

for (i in 1:nrow(genes)) {
  lvl <- genes$level[i]
  center_y <- (max_lvl - lvl + 0.5) / max_lvl
  y_top <- center_y + 0.1
  y_bot <- center_y - 0.1
  segments(genes$t_start[i], center_y, genes$t_end[i], center_y, lwd = 2)
  exons <- dplyr::filter(gtf_sub, gene_id == genes$gene_id[i], type == "exon")
  if (nrow(exons) > 0) {
    apply(exons, 1, function(e) {
      rect(as.numeric(e["start"]), y_bot, as.numeric(e["end"]), y_top, col = "black", border = NA)
    })
  }
  name <- ifelse(!is.na(genes$gene_name[i]) && genes$gene_name[i] != "", genes$gene_name[i], genes$gene_id[i])
  arrow <- ifelse(genes$strand[i] == "+", paste0(name, " >"), paste0("< ", name))
  text((genes$t_start[i] + genes$t_end[i]) / 2, y_top + 0.05, labels = arrow, cex = 0.6)
}

dev.off()
