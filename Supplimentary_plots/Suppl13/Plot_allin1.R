#!/usr/bin/env Rscript
# ============================================================
#  Regional LD triangle + GWAS + CLR + gene annotation
#  Statistic: R²
#
#  Wrapped in a function so it can be run once per region
#  (chr13 / chr22 / chr23) from a single script. Only the
#  "PARAMETERS" block differs between regions; everything else
#  (triangle building, axis stacking, gene track, plotting) is
#  identical across regions.
# ============================================================

library(data.table)
library(ggplot2)
library(dplyr)
library(scales)

# ============================================================
#  GLOBAL DEFAULTS (shared across all regions)
# ============================================================
GWAS_CHR_COL   <- "chr"
GWAS_POS_COL   <- "ps"
GWAS_PVAL_COL  <- "p_score"
GWAS_SNP_COL   <- "rs"

# ============================================================
#  PER-REGION CONFIGS — edit / add regions here
# ============================================================
GTF_MASOKO_DEFAULT <- "PATH_TO_GTF_MASOKO_DEFAULT"
GWAS_HUE_SL        <- "PATH_TO_GWAS_HUE_SL"
CLR_DIR_GENOMEWIDE <- "PATH_TO_CLR_DIR_GENOMEWIDE"
CLR_COLOURS_2POP   <- c(benthic = "steelblue2", littoral = "goldenrod2")

GWAS_NUMBER <- "PATH_TO_GWAS_NUMBER"

regions <- list(

  chr13 = list(
    LD_FILE          = "ld_results/heatmap_chr13.ld.gz",
    GTF_FILE         = GTF_MASOKO_DEFAULT,
    GWAS_FILE        = GWAS_HUE_SL,
    CHR              = 13,
    START_BP         = 4350000,
    END_BP           = 4750000,
    OUT_FILE         = "plots/LD_GWAS_CLR_chr13_4.35_4.75Mb.pdf",
    TITLE            = "Chr13: 4.35\u20134.75Mb",
    CLR_DIR          = CLR_DIR_GENOMEWIDE,
    CLR_REGION_LABEL = "chr13",
    CLR_POPS         = c("benthic", "littoral"),
    CLR_COLOURS      = CLR_COLOURS_2POP
  ),

  chr23 = list(
    LD_FILE          = "ld_results/heatmap_chr23.ld.gz",
    GTF_FILE         = GTF_MASOKO_DEFAULT,
    GWAS_FILE        = GWAS_NUMBER,
    CHR              = 23,
    START_BP         = 12100000,
    END_BP           = 12500000,
    OUT_FILE         = "plots/LD_GWAS_CLR_chr23_12_12.5Mb.pdf",
    TITLE            = "Chr23: 12\u201312.5Mb",
    CLR_DIR          = CLR_DIR_GENOMEWIDE,
    CLR_REGION_LABEL = "chr23",
    CLR_POPS         = c("benthic", "littoral"),
    CLR_COLOURS      = CLR_COLOURS_2POP
  ),

  chr22 = list(
    LD_FILE          = "ld_results/heatmap_chr22.ld.gz",
    GTF_FILE         = "PATH_TO_GTF_CHR22",
    GWAS_FILE        = GWAS_NUMBER,
    CHR              = 22,
    START_BP         = 12500000,
    END_BP           = 12900000,
    OUT_FILE         = "plots/LD_GWAS_CLR_chr22_12.5_12.9Mb.pdf",
    TITLE            = "Chr22: 12.5\u201312.9Mb",
    CLR_DIR          = CLR_DIR_GENOMEWIDE,
    CLR_REGION_LABEL = "chr22",
    CLR_POPS         = c("benthic", "littoral"),
    CLR_COLOURS      = CLR_COLOURS_2POP
  )

)

# ============================================================
#  MAIN PLOTTING FUNCTION
# ============================================================
plot_region <- function(cfg) {

  LD_FILE   <- cfg$LD_FILE
  GTF_FILE  <- cfg$GTF_FILE
  GWAS_FILE <- cfg$GWAS_FILE
  CHR       <- cfg$CHR
  START_BP  <- cfg$START_BP
  END_BP    <- cfg$END_BP
  OUT_FILE  <- cfg$OUT_FILE
  TITLE     <- cfg$TITLE

  CLR_DIR          <- cfg$CLR_DIR
  CLR_REGION_LABEL <- cfg$CLR_REGION_LABEL
  CLR_POPS         <- cfg$CLR_POPS
  CLR_COLOURS      <- cfg$CLR_COLOURS

  # ---------- 1. Read + filter LD ----------
  ld <- fread(cmd = paste0("zcat ", LD_FILE), header = TRUE)

  ld_reg <- ld %>%
    filter(CHR_A == CHR, CHR_B == CHR,
           BP_A >= START_BP, BP_A <= END_BP,
           BP_B >= START_BP, BP_B <= END_BP)

  if (nrow(ld_reg) == 0) stop("No SNPs found in this region. Check START_BP / END_BP.")

  snp_pos <- sort(unique(c(ld_reg$BP_A, ld_reg$BP_B)))
  n       <- length(snp_pos)

  mat <- matrix(NA, nrow = n, ncol = n, dimnames = list(snp_pos, snp_pos))
  for (i in seq_len(nrow(ld_reg))) {
    ra <- as.character(ld_reg$BP_A[i])
    rb <- as.character(ld_reg$BP_B[i])
    v  <- ld_reg$R2[i]
    mat[ra, rb] <- v
    mat[rb, ra] <- v
  }
  diag(mat) <- 1

  df_list <- vector("list", n * n)
  k <- 0
  for (row_i in seq_len(n)) {
    for (col_j in seq_len(n)) {
      if (col_j < row_i) next
      k <- k + 1
      df_list[[k]] <- data.frame(
        x  = (col_j + row_i) / 2,
        y  = (col_j - row_i) / 2,
        r2 = mat[row_i, col_j]
      )
    }
  }
  df_tri  <- bind_rows(df_list[seq_len(k)])
  tri_top <- (n - 1) / 2

  bp_to_idx <- approxfun(snp_pos, seq_along(snp_pos), rule = 2)
  idx_to_bp <- approxfun(seq_along(snp_pos), snp_pos, rule = 2)

  idx_breaks <- seq(1, n, length.out = min(8, n))
  bp_labels  <- function(idx) paste0(round(idx_to_bp(idx) / 1000), " kb")

  track_axis_layers <- function(tick_y, tick_labels, axis_title) {
    list(
      annotate("segment", x = 1, xend = 1,
               y = min(tick_y), yend = max(tick_y), color = "black", linewidth = 0.3),
      annotate("text", x = 1 - (n - 1) * 0.025, y = tick_y,
               label = tick_labels, hjust = 1, size = 2.8),
      annotate("text", x = 1 - (n - 1) * 0.075, y = mean(tick_y),
               label = axis_title, angle = 90, size = 3.2)
    )
  }

  # ---------- 2. Read GWAS + find top hit ----------
  gwas <- fread(GWAS_FILE, header = TRUE)

  gwas_reg <- gwas %>%
    filter(as.character(.data[[GWAS_CHR_COL]]) == as.character(CHR),
           .data[[GWAS_POS_COL]] >= START_BP, .data[[GWAS_POS_COL]] <= END_BP) %>%
    transmute(
      snp  = .data[[GWAS_SNP_COL]],
      pos  = .data[[GWAS_POS_COL]],
      pval = .data[[GWAS_PVAL_COL]],
      logp = -log10(.data[[GWAS_PVAL_COL]]),
      x    = bp_to_idx(.data[[GWAS_POS_COL]])
    )

  top_hit <- NULL
  x_top   <- NULL
  if (nrow(gwas_reg) > 0) {
    top_hit <- gwas_reg %>% slice_min(pval, n = 1, with_ties = FALSE)
    x_top   <- top_hit$x

    lead_bp <- as.character(top_hit$pos)
    if (lead_bp %in% rownames(mat)) {
      gwas_reg$r2_to_lead <- sapply(as.character(gwas_reg$pos), function(p) {
        if (p %in% colnames(mat)) mat[lead_bp, p] else NA_real_
      })
    } else {
      gwas_reg$r2_to_lead <- NA_real_
    }
  }

  # ---------- 3. Position the GWAS track above the triangle ----------
  gwas_gap    <- tri_top * 0.15
  gwas_height <- tri_top * 0.6
  logp_max    <- if (nrow(gwas_reg) > 0) max(gwas_reg$logp, na.rm = TRUE) else 1
  gwas_scale  <- gwas_height / (logp_max * 1.15)

  if (nrow(gwas_reg) > 0) {
    gwas_reg$track_y <- tri_top + gwas_gap + gwas_reg$logp * gwas_scale
    top_hit$track_y  <- tri_top + gwas_gap + top_hit$logp  * gwas_scale
  }

  gwas_breaks <- pretty(c(0, logp_max))
  gwas_breaks <- gwas_breaks[gwas_breaks >= 0 & gwas_breaks <= logp_max * 1.15]
  gwas_tick_y <- tri_top + gwas_gap + gwas_breaks * gwas_scale
  gwas_top    <- tri_top + gwas_gap + gwas_height

  # ---------- 4. Read CLR (SweepFinder2), position track above GWAS ----------
  clr_list <- lapply(CLR_POPS, function(pop) {
    f <- file.path(CLR_DIR, pop, paste0(CLR_REGION_LABEL, ".txt"))
    if (!file.exists(f)) return(NULL)
    d <- fread(f, header = TRUE)
    d %>%
      filter(location >= START_BP, location <= END_BP) %>%
      transmute(pop = pop, location = location, LR = LR, x = bp_to_idx(location))
  })
  clr_df <- bind_rows(Filter(Negate(is.null), clr_list))

  clr_gap    <- tri_top * 0.15
  clr_height <- tri_top * 0.6
  clr_max    <- if (nrow(clr_df) > 0) max(clr_df$LR, na.rm = TRUE) else 1
  clr_scale  <- clr_height / (clr_max * 1.15)
  clr_bottom <- gwas_top + clr_gap

  if (nrow(clr_df) > 0) clr_df$track_y <- clr_bottom + clr_df$LR * clr_scale

  clr_breaks <- pretty(c(0, clr_max))
  clr_breaks <- clr_breaks[clr_breaks >= 0 & clr_breaks <= clr_max * 1.15]
  clr_tick_y <- clr_bottom + clr_breaks * clr_scale
  clr_top    <- clr_bottom + clr_height

  # ---------- 5. Read GTF, position gene track below the triangle ----------
  gtf <- fread(cmd = paste0("zcat ", GTF_FILE, " | grep -v '^#'"), header = FALSE, sep = "\t", quote = "",
               col.names = c("seqname", "source", "feature", "start", "end",
                             "score", "strand", "frame", "attribute"))

  chr_candidates <- c(as.character(CHR), paste0("chr", CHR))
  chr_use <- chr_candidates[chr_candidates %in% unique(as.character(gtf$seqname))][1]
  if (is.na(chr_use)) stop("Could not match CHR to any seqname in the GTF — check chromosome naming style.")

  genes <- gtf %>%
    filter(as.character(seqname) == chr_use, feature == "gene") %>%
    filter(start <= END_BP, end >= START_BP) %>%
    mutate(
      gene_id = sub('.*gene_id "([^"]+)".*', "\\1", attribute),
      gene_name = ifelse(grepl("gene_name", attribute),
                          sub('.*gene_name "([^"]+)".*', "\\1", attribute),
                          gene_id),
      plot_start      = pmax(start, START_BP),
      plot_end        = pmin(end, END_BP),
      truncated_left  = start < START_BP,
      truncated_right = end   > END_BP,
      gene_label      = paste0(ifelse(truncated_left, "\u2026", ""), gene_name,
                                ifelse(truncated_right, "\u2026", ""))
    )

  exons <- data.frame()
  if (nrow(genes) > 0) {
    exons_all <- gtf %>%
      filter(as.character(seqname) == chr_use, feature == "exon") %>%
      filter(start <= END_BP, end >= START_BP) %>%
      mutate(
        gene_id       = sub('.*gene_id "([^"]+)".*', "\\1", attribute),
        transcript_id = sub('.*transcript_id "([^"]+)".*', "\\1", attribute),
        is_canonical  = grepl('tag "Ensembl_canonical"', attribute)
      ) %>%
      filter(gene_id %in% genes$gene_id)

    if (nrow(exons_all) > 0) {
      exons <- exons_all %>%
        group_by(gene_id) %>%
        group_modify(~ {
          if (any(.x$is_canonical)) .x %>% filter(is_canonical)
          else .x %>% filter(transcript_id == transcript_id[1])
        }) %>%
        ungroup() %>%
        mutate(plot_start = pmax(start, START_BP), plot_end = pmin(end, END_BP))
    }
  }

  gene_gap        <- tri_top * 0.05
  gene_row_height <- tri_top * 0.12
  exon_half       <- gene_row_height * 0.3
  chevrons        <- data.frame()

  if (nrow(genes) > 0) {
    genes <- genes %>%
      mutate(x1 = bp_to_idx(plot_start), x2 = bp_to_idx(plot_end)) %>%
      arrange(x1)

    row_end   <- c()
    rows      <- integer(nrow(genes))
    label_pad <- (n - 1) * 0.02
    for (i in seq_len(nrow(genes))) {
      placed <- FALSE
      reach  <- max(genes$x2[i], genes$x1[i] + nchar(genes$gene_label[i]) * 0.35)
      for (r in seq_along(row_end)) {
        if (genes$x1[i] > row_end[r] + label_pad) {
          rows[i]     <- r
          row_end[r]  <- reach
          placed      <- TRUE
          break
        }
      }
      if (!placed) {
        row_end <- c(row_end, reach)
        rows[i] <- length(row_end)
      }
    }
    genes$row     <- rows
    genes$track_y <- -(gene_gap + (genes$row - 0.5) * gene_row_height)

    if (nrow(exons) > 0) {
      exons <- exons %>%
        mutate(x1 = bp_to_idx(plot_start), x2 = bp_to_idx(plot_end)) %>%
        left_join(genes %>% select(gene_id, track_y), by = "gene_id")
    }

    arrow_spacing <- (n - 1) / 25
    chevrons <- genes %>%
      rowwise() %>%
      reframe({
        width    <- x2 - x1
        n_arrows <- max(1, min(6, round(width / arrow_spacing)))
        xs       <- seq(x1, x2, length.out = n_arrows + 2)
        xs       <- xs[-c(1, length(xs))]
        data.frame(x = xs, track_y = track_y, label = ifelse(strand == "-", "\u2039", "\u203A"))
      })
  }

  n_rows      <- if (nrow(genes) == 0) 1 else max(genes$row)
  gene_bottom <- -(gene_gap + n_rows * gene_row_height)

  # ---------- 6. Build the single combined plot ----------
  clr_legend_y <- clr_top - (clr_height * 0.09) * (seq_along(CLR_POPS) - 1)

  p <- ggplot() +
    geom_tile(data = df_tri, aes(x = x, y = y, fill = r2)) +
    {if (nrow(genes) > 0) geom_segment(data = genes, aes(x = x1, xend = x2, y = track_y, yend = track_y),
                                        color = "black", linewidth = 0.4)} +
    {if (nrow(exons) > 0) geom_rect(data = exons,
                                      aes(xmin = x1, xmax = x2, ymin = track_y - exon_half, ymax = track_y + exon_half),
                                      fill = "black", color = NA)} +
    {if (nrow(chevrons) > 0) geom_text(data = chevrons, aes(x = x, y = track_y, label = label),
                                         size = 2.6, color = "black", fontface = "bold")} +
    {if (nrow(genes) > 0) geom_text(data = genes,
                                      aes(x = (x1 + x2) / 2, y = track_y + exon_half + gene_row_height * 0.5, label = gene_label),
                                      size = 2.8, fontface = "italic")} +
    {if (nrow(genes) == 0) annotate("text", x = n / 2, y = -(gene_gap + gene_row_height),
                                      label = "no genes in window", size = 3, color = "grey50")} +
    {if (nrow(gwas_reg) > 0) geom_point(data = gwas_reg, aes(x = x, y = track_y, fill = r2_to_lead),
                                          shape = 21, color = "grey30", stroke = 0.2, size = 2)} +
    {if (!is.null(top_hit)) geom_point(data = top_hit, aes(x = x, y = track_y), color = "purple", size = 3, shape = 18)} +
    {if (!is.null(top_hit)) geom_text(data = top_hit, aes(x = x, y = track_y, label = snp),
                                        vjust = -0.8, size = 3, fontface = "bold", color = "purple")} +
    track_axis_layers(gwas_tick_y, gwas_breaks, "-log10(p) GWAS") +
    {if (!is.null(x_top)) geom_vline(xintercept = x_top, linetype = "dashed", color = "grey30")} +
    {if (nrow(clr_df) > 0) lapply(CLR_POPS, function(pop) {
      d <- filter(clr_df, pop == !!pop)
      if (nrow(d) == 0) return(NULL)
      geom_line(data = d, aes(x = x, y = track_y), color = CLR_COLOURS[pop], linewidth = 0.7)
    })} +
    {if (nrow(clr_df) > 0) lapply(seq_along(CLR_POPS), function(i) {
      annotate("text", x = n, y = clr_legend_y[i], hjust = 1, vjust = 1,
               size = 2.6, color = CLR_COLOURS[CLR_POPS[i]], label = CLR_POPS[i], fontface = "bold")
    })} +
    track_axis_layers(clr_tick_y, clr_breaks, "CLR") +
    scale_fill_gradientn(
      colors   = c("white", "yellow", "red"),
      values   = rescale(c(0, 0.5, 1)),
      limits   = c(0, 1),
      na.value = "grey70",
      name     = expression(R^2)
    ) +
    scale_x_continuous(limits = c(1 - (n - 1) * 0.1, n), breaks = idx_breaks, labels = bp_labels, expand = c(0.01, 0)) +
    scale_y_continuous(limits = c(gene_bottom - gene_row_height, clr_top + clr_gap)) +
    coord_fixed(ratio = 1) +
    labs(title = TITLE, x = paste0("Chr", CHR, " position"), y = NULL) +
    theme_minimal(base_size = 11) +
    theme(
      axis.text.y     = element_blank(),
      axis.ticks.y    = element_blank(),
      panel.grid      = element_blank(),
      plot.title      = element_text(face = "bold", hjust = 0.5),
      legend.position = "right",
      plot.margin     = margin(10, 10, 5, 10)
    )

  # ---------- 7. Save ----------
  dir.create(dirname(OUT_FILE), showWarnings = FALSE, recursive = TRUE)
  ggsave(OUT_FILE, plot = p, width = 10, height = 11 + 0.4 * n_rows, device = "pdf")
}

# ============================================================
#  RUN ALL REGIONS
# ============================================================
for (region_name in names(regions)) {
  tryCatch(
    plot_region(regions[[region_name]]),
    error = function(e) message("FAILED (", region_name, "): ", conditionMessage(e))
  )
}