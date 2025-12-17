(base) [sfy23@login-q-3 HIGH_SEQ_depth]$ cat tmem_rhoq_plot_xpEHH_IHS.R
library(dplyr)
library(rtracklayer)

# ---------- Input Files ----------
xpihs_file <- "HighSEQ_BENTHvsLITT_chr13_IHH_xpEHH_combined.tsv"
taj_file   <- "/rds/project/rds-8b3VcZwY7rY/projects/cichlid/suhaib/masoko/genome_scans_ALL/BENTHIC_LITTORAL/OVERALL_TAJI_PI_DXY/OVERALL_Tajima_chr13.Tajima.D"
gwas_file  <- "/rds/project/rds-8b3VcZwY7rY/projects/cichlid/suhaib/masoko/GWAS/masoko_HUE/output/Masoko_SUBSET_FILTERED_160_HUE_SL.assoc.txt"
gtf_file   <- "/rds/project/rds-8b3VcZwY7rY/projects/cichlid/Bettina/data/Suhaib_expression_20250422/fAstCal1.2.111.gtf.gz"

# ---------- Region of interest ----------
chrom <- "chr13"
start <- 4280000
end   <- 4850000

# TWO HIGHLIGHT REGIONS
highlight_start1 <- 4570415 - 10000
highlight_end1   <- 4570415 + 10000
highlight_start2 <- 4549161 - 10000
highlight_end2   <- 4549161 + 10000

# ---------- Bonferroni thresholds ----------
total_snps <- 1000000
bonf <- -log10(0.05 / total_snps)
suggestive <- -log10(1e-5)

# ---------- Read and process data ----------
xpihs_raw <- read.table(xpihs_file, header = TRUE)

freq_bins <- seq(0.05, 0.95, by = 0.05)

xpihs_raw$freq_bin_benth <- cut(xpihs_raw$freq_maj.benth,
                                breaks = c(0, freq_bins, 1),
                                labels = FALSE,
                                include.lowest = TRUE)
xpihs_raw <- xpihs_raw %>%
  group_by(freq_bin_benth) %>%
  mutate(unihs_std_benth = (unihs - mean(unihs, na.rm = TRUE)) / sd(unihs, na.rm = TRUE)) %>%
  ungroup()

xpihs_raw$freq_bin_litt <- cut(xpihs_raw$freq_maj.litt,
                               breaks = c(0, freq_bins, 1),
                               labels = FALSE,
                               include.lowest = TRUE)
xpihs_raw <- xpihs_raw %>%
  group_by(freq_bin_litt) %>%
  mutate(unihs_std_litt = (unihs_litt - mean(unihs_litt, na.rm = TRUE)) / sd(unihs_litt, na.rm = TRUE)) %>%
  ungroup()

xpihs <- xpihs_raw %>%
  transmute(CHR = gsub("^chr", "", chr),
            POSITION = position,
            xpEHH = logpvalue,
            unihs_benth = unihs_std_benth,
            unihs_litt = unihs_std_litt) %>%
  filter(CHR == gsub("^chr", "", chrom),
         POSITION >= start, POSITION <= end)

# Tajima's D
taj <- read.table(taj_file, header = TRUE) %>%
  transmute(CHR = gsub("^chr", "", CHROM),
            POSITION = BIN_START + 500,
            TajimaD = TajimaD) %>%
  filter(CHR == gsub("^chr", "", chrom),
         POSITION >= start, POSITION <= end)

# GWAS
gwas <- read.table(gwas_file, header = TRUE)
gwas$chr <- paste0("chr", gwas$chr)
gwas$neglogp <- -log10(gwas$p_score)
gwas_zoom <- dplyr::filter(gwas, chr == chrom, ps >= start, ps <= end)

# ---------- Gene annotation ----------
gtf <- as.data.frame(rtracklayer::import(gtf_file))
gtf_sub <- dplyr::filter(gtf, seqnames == sub("chr", "", chrom),
                         start >= start, end <= end)

genes <- dplyr::filter(gtf_sub, type == "gene") %>%
  dplyr::select(gene_id, gene_name, start, end, strand) %>%
  dplyr::distinct() %>%
  dplyr::arrange(start)

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
  for (lvl in 1:4) {
    overlap <- any(genes$level[1:(i - 1)] == lvl &
                     genes$t_start[i] <= genes$t_end[1:(i - 1)] + margin &
                     genes$t_end[i] >= genes$t_start[1:(i - 1)] - margin)
    if (!overlap) {
      genes$level[i] <- lvl
      break
    }
  }
}

max_lvl <- max(genes$level)
ann_height <- 0.4 + 0.6 * max_lvl

# ---------- UPDATED HIGHLIGHT FUNCTION ----------
add_highlight <- function() {
  rect(highlight_start1, par("usr")[3], highlight_end1, par("usr")[4],
       col = rgb(0.8, 0.8, 0.8, 0.3), border = NA)
  rect(highlight_start2, par("usr")[3], highlight_end2, par("usr")[4],
       col = rgb(0.8, 0.8, 0.8, 0.3), border = NA)
}

# ---------- Output PDF ----------
pdf("tmem_rhoq_IHS_xpEHH_taji.pdf", width = 11, height = 15 + ann_height)
layout(matrix(c(1, 2, 3, 4, 5), ncol = 1),
       heights = c(3, 3, 3, 3, ann_height))

# ---------- PANEL 1: Standardized IHS ----------
par(mar = c(0.5, 5, 2.5, 2))
ylims_ihs <- range(c(xpihs$unihs_benth, xpihs$unihs_litt), na.rm = TRUE)
plot(xpihs$POSITION, xpihs$unihs_benth, type = "n",
     xlab = "", ylab = "Standardized iHS",
     xlim = c(start, end), ylim = ylims_ihs, xaxt = "n", las = 1, cex.lab=1)
add_highlight()
lines(xpihs$POSITION, xpihs$unihs_benth, col = "steelblue2", lwd = 1.8)
lines(xpihs$POSITION, xpihs$unihs_litt, col = "goldenrod2", lwd = 1.8)
legend("topright", legend = c("Benthics", "Littorals"),
       col = c("steelblue2", "goldenrod2"), lwd = 1.8, bty = "n")

# ---------- PANEL 2: xpEHH ----------
par(mar = c(0.5, 5, 2.5, 2))
ylims1 <- range(xpihs$xpEHH, na.rm = TRUE)
plot(xpihs$POSITION, xpihs$xpEHH, type = "n",
     xlab = "", ylab = "xpEHH (logpvalue)",
     xlim = c(start, end), ylim = ylims1, xaxt = "n", las = 1, cex.lab=1)
add_highlight()
lines(xpihs$POSITION, xpihs$xpEHH, col = "black", lwd = 1.8)

# ---------- PANEL 3: Tajima's D ----------
par(mar = c(0.5, 5, 2.5, 2))
ylims2 <- range(taj$TajimaD, na.rm = TRUE)
plot(taj$POSITION, taj$TajimaD, type = "n",
     xlab = "", ylab = "Tajima's D",
     xlim = c(start, end), ylim = ylims2, xaxt = "n", las = 1, cex.lab=1)
add_highlight()
lines(taj$POSITION, taj$TajimaD, col = "black", lwd = 1.8)

# ---------- PANEL 4: GWAS ----------
par(mar = c(0.5, 5, 2.5, 2))
plot(gwas_zoom$ps, gwas_zoom$neglogp, type = "n",
     xlab = "", ylab = expression(-log[10](p)),
     xlim = c(start, end),
     ylim = range(c(gwas_zoom$neglogp, bonf)),
     xaxt = "n", las = 1, cex.lab= 1.5)
add_highlight()
abline(h = bonf, col = "black", lty = 1)
abline(h = suggestive, col = "black", lty = 2)
points(gwas_zoom$ps, gwas_zoom$neglogp, col = "black", pch = 16, cex = 0.9)

# ---------- PANEL 5: Gene Annotation ----------
par(mar = c(4, 5, 1, 2))
plot(0, type = "n", xlab = "Genomic Position (bp)", ylab = "",
     xlim = c(start, end), ylim = c(0, 1.1), axes = FALSE)
axis(1, at = pretty(c(start, end)),
     labels = format(pretty(c(start, end)), big.mark = ","))

for (i in 1:nrow(genes)) {
  lvl <- genes$level[i]
  center_y <- (max_lvl - lvl + 0.5) / max_lvl
  y_top <- center_y + 0.1
  y_bot <- center_y - 0.1
  
  segments(genes$t_start[i], center_y, genes$t_end[i], center_y, lwd = 2)

  exons <- dplyr::filter(gtf_sub, gene_id == genes$gene_id[i], type == "exon")
  if (nrow(exons) > 0) {
    apply(exons, 1, function(e) {
      rect(as.numeric(e["start"]), y_bot, as.numeric(e["end"]), y_top,
           col = "black", border = NA)
    })
  }

  name <- ifelse(!is.na(genes$gene_name[i]) && genes$gene_name[i] != "",
                 genes$gene_name[i], genes$gene_id[i])
  arrow <- ifelse(genes$strand[i] == "+",
                  paste0(name, " >"),
                  paste0("< ", name))
  text((genes$t_start[i] + genes$t_end[i]) / 2, y_top + 0.05,
       labels = arrow, cex = 1)
}

dev.off()