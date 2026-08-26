library(ggplot2)
library(dplyr)
library(patchwork)

# =============================================================================
# ── USER PARAMETERS ──────────────────────────────────────────────────────────
# =============================================================================

sf2_dir  <- "sf2_results_genomewide"
out_dir  <- "PDFs"
out_file <- "SF2_manhattan_benthic_littoral.pdf"

CHROMS <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,22,23)

chrom_sizes <- c(
  `1`=41162407,  `2`=38215696,  `3`=52512415,  `4`=33433079,
  `5`=38678279,  `6`=41428020,  `7`=68397778,  `8`=25827708,
  `9`=35913139,  `10`=34074121, `11`=36676067, `12`=38669361,
  `13`=33145951, `14`=39985354, `15`=40222765, `16`=35801853,
  `17`=39448915, `18`=35852588, `19`=30863130, `20`=31467755,
  `22`=35011464, `23`=43921804
)

# Highlight region
highlight_chr   <- 23
highlight_start <- 12200000
highlight_end   <- 12320000

# =============================================================================
# ── LOAD DATA ────────────────────────────────────────────────────────────────
# =============================================================================

load_pop <- function(pop) {
  bind_rows(lapply(CHROMS, function(chr) {
    f <- file.path(sf2_dir, pop, paste0("chr", chr, ".txt"))
    if (!file.exists(f)) return(NULL)
    df <- tryCatch(
      read.table(f, header = TRUE, sep = "\t"),
      error = function(e) NULL
    )
    if (is.null(df) || nrow(df) == 0) return(NULL)
    names(df)[names(df) == "location"] <- "pos"
    df$chr <- chr
    df
  }))
}


# =============================================================================
# ── THRESHOLDS ───────────────────────────────────────────────────────────────
# =============================================================================

# Standard for SF2 Manhattan: empirical 95th and 99th percentile genome-wide
# (Nielsen et al. 2005 original SF paper; also used in Pavlidis et al. 2013)
t_ben <- quantile(ben$LR, c(0.975, 0.995), na.rm = TRUE)
t_lit <- quantile(lit$LR, c(0.975, 0.995), na.rm = TRUE)


# =============================================================================
# ── CUMULATIVE POSITIONS ─────────────────────────────────────────────────────
# =============================================================================

offsets <- c(0, cumsum(chrom_sizes[as.character(CHROMS[-length(CHROMS)])]))
names(offsets) <- as.character(CHROMS)

add_cumpos <- function(df, col1, col2) {
  df$cum_pos   <- df$pos + offsets[as.character(df$chr)]
  df$chrom_col <- ifelse(match(df$chr, CHROMS) %% 2 == 1, col1, col2)
  df
}

# Steelblue shades for benthic, goldenrod shades for littoral
ben <- add_cumpos(ben, "#1A5276", "#5DADE2")   # dark/light steelblue
lit <- add_cumpos(lit, "#9A6B00", "#F0B429")   # dark/light goldenrod

chrom_mids <- offsets + chrom_sizes[as.character(CHROMS)] / 2

highlight_cum_start <- highlight_start + offsets[as.character(highlight_chr)]
highlight_cum_end   <- highlight_end   + offsets[as.character(highlight_chr)]

x_max <- max(ben$cum_pos, na.rm = TRUE)

# =============================================================================
# ── PLOT FUNCTION ─────────────────────────────────────────────────────────────
# =============================================================================

make_manhattan <- function(df, title, title_col, t95, t99) {

  y_max <- pmin(
    quantile(df$LR, 0.9999, na.rm = TRUE) * 1.15,
    max(df$LR, na.rm = TRUE) * 1.05
  )

  ggplot(df, aes(x = cum_pos, y = LR, colour = chrom_col)) +
    annotate("rect",
             xmin = highlight_cum_start, xmax = highlight_cum_end,
             ymin = -Inf, ymax = Inf,
             fill = "#FFDD99", alpha = 0.5) +
    geom_point(size = 0.3, alpha = 0.6) +
    geom_hline(yintercept = t95, linetype = "dashed",
               colour = "#F4A261", linewidth = 0.6) +
    geom_hline(yintercept = t99, linetype = "dashed",
               colour = "#E63946", linewidth = 0.6) +
    annotate("text", x = x_max, y = t99,
             label = paste0("99th (", round(t99, 1), ")"),
             hjust = 1, vjust = -0.4,
             size = 2.8, colour = "#E63946") +
    annotate("text", x = x_max, y = t95,
             label = paste0("95th (", round(t95, 1), ")"),
             hjust = 1, vjust = -0.4,
             size = 2.8, colour = "#F4A261") +
    scale_colour_identity() +
    scale_x_continuous(
      breaks = chrom_mids,
      labels = CHROMS,
      expand = expansion(mult = 0.01)
    ) +
    scale_y_continuous(
      limits = c(0, y_max),
      expand = expansion(mult = c(0, 0.02))
    ) +
    labs(title = title, x = NULL, y = "CLR") +
    theme_classic(base_size = 10) +
    theme(
      axis.text.x     = element_text(size = 7),
      plot.title      = element_text(face = "bold",
                                     colour = title_col, size = 11),
      legend.position = "none",
      plot.margin     = margin(4, 8, 4, 8)
    )
}

# =============================================================================
# ── ASSEMBLE AND SAVE ─────────────────────────────────────────────────────────
# =============================================================================

p_ben <- make_manhattan(ben, "Benthic",  "steelblue",  t_ben[1], t_ben[2])
p_lit <- make_manhattan(lit, "Littoral", "goldenrod3", t_lit[1], t_lit[2])

combined <- (p_ben / p_lit) +
  plot_annotation(
    title    = "SweepFinder2 — Composite Likelihood Ratio",
    subtitle = "Dashed lines: empirical 97.5th and 99.5th genome-wide percentiles | yellow = chr23 peak region",
    theme    = theme(
      plot.title    = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(size = 9, colour = "grey40")
    )
  )

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
out_path <- file.path(out_dir, out_file)

ggsave(out_path, combined,
       width = 18, height = 7,
       device = cairo_pdf, limitsize = FALSE)


