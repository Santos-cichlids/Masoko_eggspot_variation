library(dplyr)
library(rtracklayer)

# ---------- Input Files ----------
xpihs_file <- "/rds/project/rds-8b3VcZwY7rY/projects/cichlid/suhaib/masoko/genome_scans_ALL/BENTHIC_LITTORAL/EHH_IHH_Selection/biallelics/HIGH_SEQ_depth/Downsampled_xpEHH_IHS_POLARIZED/HighSEQ_BENTHvsLITT_chr13_IHH_xpEHH_combined_downsampled_POLARIZED.tsv"
taj_file   <- "/rds/project/rds-8b3VcZwY7rY/projects/cichlid/suhaib/masoko/genome_scans_ALL/BENTHIC_LITTORAL/OVERALL_TAJI_PI_DXY/OVERALL_Tajima_10K_chr23.Tajima.DY"
gwas_file  <- "/rds/project/rds-8b3VcZwY7rY/projects/cichlid/suhaib/masoko/GWAS/masoko_HUE/output/Masoko_SUBSET_FILTERED_160_HUE_SL.assoc.txt"
#gtf_file   <- "/rds/project/rds-8b3VcZwY7rY/projects/cichlid/suhaib/masoko/genome_scans_ALL/BENTHIC_LITTORAL/EHH_IHH_Selection/biallelics/HIGH_SEQ_depth/Astatotilapia_calliptera.fAstCal1.3.115.gtf.gz" (used only for lncRNA annotation)
gtf_file   <- "/rds/project/rds-8b3VcZwY7rY/projects/cichlid/Bettina/data/Suhaib_expression_20250422/fAstCal1.2.111.gtf.gz"


# ---------- Region of interest ----------
chrom <- "chr13"
start <- 4350000
end   <- 4750000

# TWO HIGHLIGHT REGIONS
highlight_start1 <- 4570415 - 10000
highlight_end1   <- 4570415 + 10000
highlight_start2 <- 4549161 - 10000
highlight_end2   <- 4549161 + 10000

# ---------- Bonferroni thresholds ----------
total_snps <- 1691704 #[from GWAS file]
bonf <- -log10(0.05 / total_snps)
suggestive <- -log10(1e-5)

# ---------- xpEHH Thresholds ----------
xpehh_upper_975 <- 1.965812
xpehh_lower_25  <- -2.072641
xpehh_upper_995 <- 2.743514
xpehh_lower_05  <- -3.064567

# ---------- Read and process data ----------
xpihs_raw <- read.table(xpihs_file, header = TRUE)

# The downsampled polarized file has different column names
# Extract relevant columns: position, xpehh_benth_litt, xpehh_benth_litt_sd, unihs, unihs_sd, unihs_litt, unihs_litt_sd

xpihs <- xpihs_raw %>%
  transmute(CHR = gsub("^chr", "", chr),
            POSITION = position,
            xpEHH = xpehh_benth_litt,
            xpEHH_sd = xpehh_benth_litt_sd,
            unihs_benth = unihs,
            unihs_litt = unihs_litt) %>%
  filter(CHR == gsub("^chr", "", chrom),
         POSITION >= start, POSITION <= end)

# Tajima's D
taj <- read.table(taj_file, header = TRUE) %>%
  transmute(CHR = gsub("^chr", "", CHROM),
            POSITION = BIN_START + 5000,
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
pdf("StandardisedPanel_tmem_rhoq_IHS_xpEHH_taji.pdf", width = 11, height = 15 + ann_height)
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

# ---------- PANEL 2: xpEHH (FIXED Y-AXIS WITH SD) ----------
par(mar = c(0.5, 5, 2.5, 2))
plot(xpihs$POSITION, xpihs$xpEHH, type = "n",
     xlab = "", ylab = "xpEHH",
     xlim = c(start, end), ylim = c(-3.5, 3.5), xaxt = "n", las = 1, cex.lab=1)
add_highlight()
# Add threshold lines
abline(h = 0, col = "gray50", lty = 2, lwd = 1)
abline(h = xpehh_upper_975, col = "blue", lty = 2, lwd = 1.2)
abline(h = xpehh_lower_25, col = "blue", lty = 2, lwd = 1.2)
abline(h = xpehh_upper_995, col = "red", lty = 2, lwd = 1.2)
abline(h = xpehh_lower_05, col = "red", lty = 2, lwd = 1.2)
# Add SD shaded region (mean ± SD)
xpihs_clean <- xpihs[!is.na(xpihs$xpEHH) & !is.na(xpihs$xpEHH_sd), ]
if(nrow(xpihs_clean) > 0) {
  polygon(c(xpihs_clean$POSITION, rev(xpihs_clean$POSITION)),
          c(xpihs_clean$xpEHH - xpihs_clean$xpEHH_sd, 
            rev(xpihs_clean$xpEHH + xpihs_clean$xpEHH_sd)),
          col = rgb(0.5, 0.5, 0.5, 0.3), border = NA)
}
# Plot mean line
lines(xpihs$POSITION, xpihs$xpEHH, col = "black", lwd = 1.8)

# ---------- PANEL 3: Tajima's D (FIXED Y-AXIS) ----------
par(mar = c(0.5, 5, 2.5, 2))
plot(taj$POSITION, taj$TajimaD, type = "n",
     xlab = "", ylab = "Tajima's D",
     xlim = c(start, end), ylim = c(-3.5, 3.5), xaxt = "n", las = 1, cex.lab=1)
add_highlight()
# Add horizontal line at 0
abline(h = 0, col = "gray50", lty = 2, lwd = 1)
# Plot data
lines(taj$POSITION, taj$TajimaD, col = "black", lwd = 1.8)

# ---------- PANEL 4: GWAS ----------
par(mar = c(0.5, 5, 2.5, 2))
plot(gwas_zoom$ps, gwas_zoom$neglogp, type = "n",
     xlab = "", ylab = expression(-log[10](p)),
     xlim = c(start, end),
     ylim = c(0, 9),,
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