#!/usr/bin/env Rscript

  library(data.table)
  library(ggplot2)

BLOCKS_DIR <- "ld_results/blocks"
CHROMS     <- setdiff(1:23, 21)
OUT_DIR    <- "plots"
dir.create(OUT_DIR, showWarnings = FALSE)

# Load all chromosomes
blocks <- rbindlist(lapply(CHROMS, function(chr) {
  f <- file.path(BLOCKS_DIR, paste0("cichlid_chr", chr, ".blocks.det"))
  if (!file.exists(f)) { warning("Missing: ", f); return(NULL) }
  dt <- fread(f, header = TRUE)
  if (nrow(dt) == 0) { message("Skipping chr", chr, " — empty"); return(NULL) }
  setnames(dt, trimws(names(dt)))
  dt[, .(CHR = chr, block_size_kb = (BP2 - BP1) / 1000)]
}))

# Summary stats
mean_kb   <- round(mean(blocks$block_size_kb), 2)
median_kb <- round(median(blocks$block_size_kb), 2)
cat("Mean block size:  ", mean_kb, "kb\n")
cat("Median block size:", median_kb, "kb\n")

# Per-chromosome means
chr_means <- blocks[, .(mean_kb = mean(block_size_kb)), by = CHR]

# Plot
p <- ggplot(blocks, aes(x = factor(CHR), y = block_size_kb)) +
  geom_boxplot(
    outlier.shape = NA,
    fill          = "#63B995",
    colour        = "#2d6e52",
    linewidth     = 0.4,
    coef          = 1.5
  ) +
  geom_point(
    data   = chr_means,
    aes(x  = factor(CHR), y = mean_kb),
    shape  = 16,
    size   = 2,
    colour = "#085041"
  ) +
  annotate("text",
    x = 1, y = Inf,
    label  = paste0("Mean = ", mean_kb, " kb\nMedian = ", median_kb, " kb"),
    hjust  = 0, vjust = 1.4, size = 3.5, colour = "#2d6e52"
  ) +
  coord_cartesian(ylim = c(0, quantile(blocks$block_size_kb, 0.99))) +
  scale_x_discrete(name = "Chromosome") +
  scale_y_continuous(name = "Block size (kb)") +
  labs(
    title    = "LD block size per chromosome",
    subtitle = "Boxes = IQR, whiskers = 1.5×IQR, point = mean"
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title         = element_text(face = "bold", size = 11),
    plot.subtitle      = element_text(colour = "grey50", size = 8),
    panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.3)
  )

ggsave(file.path(OUT_DIR, "TMP2_LD_block_sizes_per_chr.pdf"),
       plot = p, width = 12, height = 5, device = cairo_pdf)
ggsave(file.path(OUT_DIR, "TMP2_LD_block_sizes_per_chr.png"),
       plot = p, width = 12, height = 5, dpi = 300)