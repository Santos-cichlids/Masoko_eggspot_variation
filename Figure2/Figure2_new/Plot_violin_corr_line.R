# =============================================================================
# Combined Expression Plotting Script
# Produces: (1) violin + box plots per gene (all 18 genes + miR10)
#           (2) line plots (mean F vs M) per tissue for focal genes
# Data: pre-normalised TPM counts — no additional normalisation applied
# Significance: Wilcoxon raw p-value, converted to stars (ns, *, **, ***)
# =============================================================================

library(tidyverse)
library(ggrepel)

# =============================================================================
# GENE MAPPINGS
# violin plots: all 18 genes + miR10
# line plots:   5 focal genes only
# =============================================================================

gene_mapping_violin <- data.frame(
  ensembl = c("ENSACLG00000019781", "ENSACLG00000017026",
              "ENSACLG00000034115", "ENSACLG00000026836", "ENSACLG00000026955",
              "ENSACLG00000019614", "ENSACLG00000014805", "ENSACLG00000001747",
              "ENSACLG00000019599", "ENSACLG00000019578", "ENSACLG00000005191",
              "ENSACLG00000004021", "ENSACLG00000024080", "ENSACLG00000024228",
              "ENSACLG00000024190", "ENSACLG00000025445", "ENSACLG00000029630"),
  gene_name = c("gnas", "tmem254a", "rhoq", "aldh9a1", "tmco1",
                "prelid3b", "ntrk3a", "atp8b1", "aurka", "slc32a1", "ptpa",
                "kcnh2", "herc2", "gabrg3", "oca2", "plxdc1", "lncRNA"),
  stringsAsFactors = FALSE
)
#plxdc1=ENSACLG00000025445
gene_mapping_line <- data.frame(
  ensembl = c("ENSACLG00000017026", "ENSACLG00000024190", "ENSACLG00000034115",
              "ENSACLG00000010391", "ENSACLG00000029630"),
  gene_name = c("tmem254", "oca2", "rhoq", "miR10", "lncRNA"),
  stringsAsFactors = FALSE
)

# colours for line plot genes
gene_colours <- c(
  "tmem254" = "#1B9E77",
  "oca2"    = "#E7298A",
  "rhoq"    = "#E6AB02",
  "miR10"   = "#7570B3",
  "lncRNA"  = "#A6761D"
)

# =============================================================================
# COLOURS — non-stereotypical sex colours
# =============================================================================

sex_colors    <- c("F" = "#2E8B57", "M" = "yellow3")  
tissue_colors <- c("Anal Fin" = "#D55E00", "Muscle" = "#0072B2", "Gonad" = "#009E73")

# =============================================================================
# HELPER: parse metadata from sample column names
# masoko_{tissue_code}_{sex_code}_rep{n}
# =============================================================================

parse_sample_meta <- function(sample_names) {
  tissue_map <- c(anf = "Anal Fin", gon = "Gonad", mus = "Muscle")
  sex_map    <- c(fem = "F", mal = "M")
  parts      <- strsplit(sample_names, "_")
  data.frame(
    sample     = sample_names,
    tissue     = factor(tissue_map[sapply(parts, `[`, 2)],
                        levels = c("Anal Fin", "Muscle", "Gonad")),
    sex        = factor(sex_map[sapply(parts, `[`, 3)],
                        levels = c("F", "M")),
    individual = sapply(parts, `[`, 4),
    stringsAsFactors = FALSE
  )
}

# =============================================================================
# HELPER: p-value to significance stars
# =============================================================================

pval_to_stars <- function(p) {
  if (is.na(p))  return(NA_character_)
  if (p < 0.001) return("***")
  if (p < 0.01)  return("**")
  if (p < 0.05)  return("*")
  return("ns")
}

# =============================================================================
# LOAD DATA
# =============================================================================

fc_norm <- read.table("./TPM_counts_lncRNA.txt", header = TRUE,
                      sep = "\t", check.names = FALSE)
colnames(fc_norm)[1] <- "ensembl"

count_data_miR10 <- read.table("Masoko_NormCounts_smallRNA.txt", header = TRUE,
                               sep = "\t", check.names = FALSE)
colnames(count_data_miR10)[1] <- "gene"

# =============================================================================
# METADATA
# =============================================================================

meta_fc    <- parse_sample_meta(colnames(fc_norm)[-1])
meta_miR10 <- parse_sample_meta(colnames(count_data_miR10)[-1])

cat("fc_norm samples    :", nrow(meta_fc),    "\n")
cat("miR10 samples      :", nrow(meta_miR10), "\n")
cat("Shared samples     :", length(intersect(meta_fc$sample, meta_miR10$sample)), "\n\n")

# =============================================================================
# THEME
# =============================================================================

custom_theme <- theme_classic(base_size = 14) + theme(
  axis.text       = element_text(color = "black"),
  plot.title      = element_text(face = "bold", hjust = 0.5),
  axis.text.x     = element_text(angle = 45, hjust = 1),
  legend.position = "right"
)

# =============================================================================
# BUILD LONG DATA — fc_norm genes (violin + line)
# drop zeros, log2
# =============================================================================

fc_long <- fc_norm %>%
  dplyr::filter(ensembl %in% gene_mapping_violin$ensembl) %>%
  tidyr::pivot_longer(-ensembl, names_to = "sample", values_to = "raw_count") %>%
  dplyr::left_join(gene_mapping_violin, by = "ensembl") %>%
  dplyr::left_join(meta_fc, by = "sample") %>%
  tidyr::drop_na(tissue) %>%
  dplyr::filter(raw_count > 0) %>%
  dplyr::mutate(
    norm_count = log2(raw_count),
    tissue     = factor(as.character(tissue), levels = c("Anal Fin", "Muscle", "Gonad")),
    sex        = factor(as.character(sex),    levels = c("F", "M"))
  ) %>%
  dplyr::ungroup() %>%
  as.data.frame()

# miR10 long data
miR10_long <- count_data_miR10 %>%
  dplyr::filter(gene == "ENSACLG00000010391") %>%
  tidyr::pivot_longer(-gene, names_to = "sample", values_to = "raw_count") %>%
  dplyr::left_join(meta_miR10, by = "sample") %>%
  tidyr::drop_na(tissue) %>%
  dplyr::filter(raw_count > 0) %>%
  dplyr::mutate(
    norm_count = log2(raw_count),
    gene_name  = "miR10",
    tissue     = factor(as.character(tissue), levels = c("Anal Fin", "Muscle", "Gonad")),
    sex        = factor(as.character(sex),    levels = c("F", "M"))
  ) %>%
  dplyr::ungroup() %>%
  as.data.frame()

# =============================================================================
# SIGNIFICANCE 
# =============================================================================

# -- Tissue-level significance (violin only, per gene, pooled across sex) --
tissue_sig_df <- do.call(rbind, lapply(
  split(fc_long, fc_long$ensembl), function(sub) {
    data.frame(
      ensembl  = sub$ensembl[1],
      sig_AF_M = pval_to_stars(tryCatch(
        wilcox.test(sub$norm_count[sub$tissue == "Anal Fin"],
                    sub$norm_count[sub$tissue == "Muscle"])$p.value,  error = \(e) NA_real_)),
      sig_AF_G = pval_to_stars(tryCatch(
        wilcox.test(sub$norm_count[sub$tissue == "Anal Fin"],
                    sub$norm_count[sub$tissue == "Gonad"])$p.value,   error = \(e) NA_real_)),
      sig_M_G  = pval_to_stars(tryCatch(
        wilcox.test(sub$norm_count[sub$tissue == "Muscle"],
                    sub$norm_count[sub$tissue == "Gonad"])$p.value,   error = \(e) NA_real_)),
      stringsAsFactors = FALSE
    )
  }
))

# -- Sex-level significance (both plot types, per gene x tissue) --
# fc_norm genes
sex_sig_fc <- do.call(rbind, lapply(
  split(fc_long, list(fc_long$ensembl, as.character(fc_long$tissue))), function(sub) {
    if (nrow(sub) == 0) return(NULL)
    data.frame(
      ensembl      = as.character(sub$ensembl[1]),
      gene_name    = as.character(sub$gene_name[1]),
      tissue       = as.character(sub$tissue[1]),
      significance = pval_to_stars(tryCatch(
        wilcox.test(sub$norm_count[sub$sex == "F"],
                    sub$norm_count[sub$sex == "M"])$p.value, error = \(e) NA_real_)),
      stringsAsFactors = FALSE
    )
  }
))

# miR10 sex significance
sex_sig_miR10 <- do.call(rbind, lapply(
  split(miR10_long, as.character(miR10_long$tissue)), function(sub) {
    if (nrow(sub) == 0) return(NULL)
    data.frame(
      ensembl      = "ENSACLG00000010391",
      gene_name    = "miR10",
      tissue       = as.character(sub$tissue[1]),
      significance = pval_to_stars(tryCatch(
        wilcox.test(sub$norm_count[sub$sex == "F"],
                    sub$norm_count[sub$sex == "M"])$p.value, error = \(e) NA_real_)),
      stringsAsFactors = FALSE
    )
  }
))

# Combined sex sig for line plots
sex_sig_all <- rbind(sex_sig_fc, sex_sig_miR10)
sex_sig_all <- sex_sig_all[sex_sig_all$gene_name %in% gene_mapping_line$gene_name, ]

print(sex_sig_fc[sex_sig_fc$ensembl == "ENSACLG00000025445", ])

# plot_data for violin plots (tissue sig joined, sex sig passed separately)
plot_data <- merge(fc_long, tissue_sig_df, by = "ensembl", all.x = TRUE)
plot_data$tissue <- factor(as.character(plot_data$tissue), levels = c("Anal Fin", "Muscle", "Gonad"))
plot_data$sex    <- factor(as.character(plot_data$sex),    levels = c("F", "M"))

# =============================================================================
# VIOLIN PLOT FUNCTION
# =============================================================================

make_violin_plot <- function(df, gene_label, sex_sig_lookup,
                             y_label = "Normalised Counts (log2)") {
  
  ts_row   <- as.data.frame(df[, c("sig_AF_M", "sig_AF_G", "sig_M_G")])[1, ]
  sig_AF_M <- as.character(ts_row$sig_AF_M)
  sig_AF_G <- as.character(ts_row$sig_AF_G)
  sig_M_G  <- as.character(ts_row$sig_M_G)
  
  get_sex_sig <- function(tis) {
    rows <- sex_sig_lookup$significance[sex_sig_lookup$tissue == tis]
    if (length(rows) == 0 || all(is.na(rows))) return("ns")
    as.character(rows[1])
  }
  
  df$tissue <- factor(as.character(df$tissue), levels = c("Anal Fin", "Muscle", "Gonad"))
  df$sex    <- factor(as.character(df$sex),    levels = c("F", "M"))
  y_max <- max(df$norm_count, na.rm = TRUE)
  y_pos <- y_max + 0.5
  
  p <- ggplot(df, aes(x = tissue, y = norm_count, fill = sex)) +
    geom_violin(trim = FALSE, color = "black", alpha = 0.15,
                position = position_dodge(0.8), width = 1) +
    geom_boxplot(width = 0.12, outlier.shape = NA, alpha = 0.5,
                 color = "black", position = position_dodge(0.8)) +
    geom_point(aes(color = sex),
               position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.8),
               size = 3, alpha = 0.85) +
    scale_fill_manual(values  = sex_colors) +
    scale_color_manual(values = sex_colors) +
    labs(x = "Tissue", y = y_label, title = gene_label) +
    custom_theme +
    guides(fill  = guide_legend(title = "Sex"),
           color = guide_legend(title = "Sex"))
  
  # Tissue comparison brackets
  if (!is.na(sig_AF_M) && sig_AF_M != "ns") {
    p <- p +
      annotate("segment", x = 1, xend = 2, y = y_pos, yend = y_pos,
               color = "black", linewidth = 0.8) +
      annotate("segment", x = 1, xend = 1, y = y_pos - 0.08, yend = y_pos,
               color = "black", linewidth = 0.8) +
      annotate("segment", x = 2, xend = 2, y = y_pos - 0.08, yend = y_pos,
               color = "black", linewidth = 0.8) +
      annotate("text", x = 1.5, y = y_pos + 0.15, label = sig_AF_M,
               color = "black", size = 4, fontface = "bold")
    y_pos <- y_pos + 0.6
  }
  if (!is.na(sig_AF_G) && sig_AF_G != "ns") {
    p <- p +
      annotate("segment", x = 1, xend = 3, y = y_pos, yend = y_pos,
               color = "black", linewidth = 0.8) +
      annotate("segment", x = 1, xend = 1, y = y_pos - 0.08, yend = y_pos,
               color = "black", linewidth = 0.8) +
      annotate("segment", x = 3, xend = 3, y = y_pos - 0.08, yend = y_pos,
               color = "black", linewidth = 0.8) +
      annotate("text", x = 2, y = y_pos + 0.15, label = sig_AF_G,
               color = "black", size = 4, fontface = "bold")
    y_pos <- y_pos + 0.6
  }
  if (!is.na(sig_M_G) && sig_M_G != "ns") {
    p <- p +
      annotate("segment", x = 2, xend = 3, y = y_pos, yend = y_pos,
               color = "black", linewidth = 0.8) +
      annotate("segment", x = 2, xend = 2, y = y_pos - 0.08, yend = y_pos,
               color = "black", linewidth = 0.8) +
      annotate("segment", x = 3, xend = 3, y = y_pos - 0.08, yend = y_pos,
               color = "black", linewidth = 0.8) +
      annotate("text", x = 2.5, y = y_pos + 0.15, label = sig_M_G,
               color = "black", size = 4, fontface = "bold")
  }
  
  # Within-tissue sex brackets
  tissue_xpos <- c("Anal Fin" = 1, "Muscle" = 2, "Gonad" = 3)
  for (tis in names(tissue_xpos)) {
    sig_val  <- get_sex_sig(tis)
    xi       <- tissue_xpos[tis]
    tis_vals <- df$norm_count[as.character(df$tissue) == tis]
    if (length(tis_vals) == 0 || all(is.na(tis_vals))) next
    y_sex <- max(tis_vals, na.rm = TRUE) + 0.15
    if (!is.na(sig_val) && sig_val != "ns") {
      p <- p +
        annotate("segment", x = xi - 0.2, xend = xi + 0.2,
                 y = y_sex, yend = y_sex, color = "grey30", linewidth = 0.5) +
        annotate("text", x = xi, y = y_sex + 0.15,
                 label = sig_val, color = "grey30", size = 4, fontface = "bold")
    }
  }
  p
}

# =============================================================================
# GENERATE VIOLIN PLOTS (fc_norm genes + miR10)
# =============================================================================

dir.create("./MIG_TPM", showWarnings = FALSE)

# fc_norm genes
plot_list   <- base::split(plot_data, plot_data$gene_name)
violin_plots <- lapply(names(plot_list), function(gname) {
  df          <- plot_list[[gname]]
  ensembl_id  <- as.character(df$ensembl[1])
  sex_sig_sub <- sex_sig_fc[sex_sig_fc$ensembl == ensembl_id,
                            c("tissue", "significance")]
  make_violin_plot(df, gname, sex_sig_sub)
})
names(violin_plots) <- names(plot_list)

# miR10 violin
miR10_ts_sig <- data.frame(
  sig_AF_M = pval_to_stars(tryCatch(
    wilcox.test(miR10_long$norm_count[miR10_long$tissue == "Anal Fin"],
                miR10_long$norm_count[miR10_long$tissue == "Muscle"])$p.value,  error = \(e) NA_real_)),
  sig_AF_G = pval_to_stars(tryCatch(
    wilcox.test(miR10_long$norm_count[miR10_long$tissue == "Anal Fin"],
                miR10_long$norm_count[miR10_long$tissue == "Gonad"])$p.value,   error = \(e) NA_real_)),
  sig_M_G  = pval_to_stars(tryCatch(
    wilcox.test(miR10_long$norm_count[miR10_long$tissue == "Muscle"],
                miR10_long$norm_count[miR10_long$tissue == "Gonad"])$p.value,   error = \(e) NA_real_)),
  stringsAsFactors = FALSE
)
miR10_plot_data <- miR10_long
miR10_plot_data$sig_AF_M <- miR10_ts_sig$sig_AF_M
miR10_plot_data$sig_AF_G <- miR10_ts_sig$sig_AF_G
miR10_plot_data$sig_M_G  <- miR10_ts_sig$sig_M_G

violin_plots[["miR10"]] <- make_violin_plot(
  miR10_plot_data, "miR10",
  sex_sig_miR10[, c("tissue", "significance")]
)

# Save individual violin PDFs
for (gname in names(violin_plots)) {
  clean_name <- gsub("[^A-Za-z0-9_]", "_", gname)
  ggsave(paste0("./MIG_TPM/violin_", clean_name, ".pdf"),
         violin_plots[[gname]], width = 10, height = 7)
}
pdf("./MIG_TPM/all_violins.pdf", width = 10, height = 7)
for (p in violin_plots) print(p)
dev.off()

# =============================================================================
# BUILD DATA FOR LINE PLOTS
# uses fc_long (already built) + miR10_long, filtered to focal genes
# =============================================================================

# fc_norm focal genes (subset from fc_long, rename to match line gene_mapping)
line_fc_long <- fc_long[fc_long$ensembl %in%
                          gene_mapping_line$ensembl[gene_mapping_line$gene_name != "miR10"], ]
# remap gene names to line plot names (tmem254a -> tmem254, rhoq stays rhoq etc.)
line_fc_long <- merge(line_fc_long[, !colnames(line_fc_long) %in% "gene_name"],
                      gene_mapping_line, by = "ensembl")

line_all <- rbind(
  line_fc_long[, c("gene_name", "sex", "tissue", "norm_count")],
  miR10_long[,  c("gene_name", "sex", "tissue", "norm_count")]
)
colnames(line_all)[colnames(line_all) == "norm_count"] <- "log_norm_counts"

# Mean per gene x sex x tissue
sex_mean_data <- do.call(rbind, lapply(
  split(line_all, list(line_all$gene_name,
                       as.character(line_all$sex),
                       as.character(line_all$tissue))), function(sub) {
                         if (nrow(sub) == 0) return(NULL)
                         data.frame(
                           gene_name       = as.character(sub$gene_name[1]),
                           sex             = as.character(sub$sex[1]),
                           tissue          = as.character(sub$tissue[1]),
                           mean_log_counts = mean(sub$log_norm_counts, na.rm = TRUE),
                           stringsAsFactors = FALSE
                         )
                       }
))
sex_mean_data$sex    <- factor(sex_mean_data$sex,    levels = c("F", "M"))
sex_mean_data$tissue <- factor(sex_mean_data$tissue, levels = c("Anal Fin", "Muscle", "Gonad"))

# Significance label positions for line plots
sig_positions <- do.call(rbind, lapply(
  split(sex_mean_data, list(sex_mean_data$gene_name,
                            as.character(sex_mean_data$tissue))), function(sub) {
                              if (nrow(sub) == 0) return(NULL)
                              gn  <- as.character(sub$gene_name[1])
                              tis <- as.character(sub$tissue[1])
                              sig_row <- sex_sig_all[sex_sig_all$gene_name == gn &
                                                       sex_sig_all$tissue   == tis, ]
                              sig_lbl <- if (nrow(sig_row) > 0) as.character(sig_row$significance[1]) else "ns"
                              data.frame(
                                gene_name  = gn,
                                tissue     = tis,
                                mid_x      = 1.5,
                                mid_y      = mean(sub$mean_log_counts, na.rm = TRUE),
                                sig_label  = sig_lbl,
                                stringsAsFactors = FALSE
                              )
                            }
))
sig_positions$tissue <- factor(sig_positions$tissue, levels = c("Anal Fin", "Muscle", "Gonad"))

# =============================================================================
# LINE PLOT FUNCTION (one per tissue)
# =============================================================================

make_line_plot <- function(tiss) {
  
  tiss_data   <- sex_mean_data[as.character(sex_mean_data$tissue) == tiss, ]
  tiss_jitter <- line_all[as.character(line_all$tissue) == tiss, ]
  tiss_sig    <- sig_positions[as.character(sig_positions$tissue) == tiss, ]
  tiss_labels <- tiss_data[as.character(tiss_data$sex) == "M", ]
  
  ggplot() +
    geom_jitter(
      data = tiss_jitter,
      aes(sex, log_norm_counts, color = gene_name),
      width = 0.06, size = 2, alpha = 0.3, show.legend = FALSE
    ) +
    geom_line(
      data = tiss_data,
      aes(sex, mean_log_counts, group = gene_name, color = gene_name),
      linewidth = 1.2
    ) +
    geom_point(
      data = tiss_data,
      aes(sex, mean_log_counts, color = gene_name),
      size = 3.5
    ) +
    geom_text(
      data = tiss_sig,
      aes(mid_x, mid_y, label = sig_label, color = gene_name),
      size = 5, vjust = -0.5, show.legend = FALSE
    ) +
    geom_text_repel(
      data = tiss_labels,
      aes(sex, mean_log_counts, label = gene_name, color = gene_name),
      size = 4, nudge_x = 0.3, hjust = 0,
      segment.size = 0.3, show.legend = FALSE
    ) +
    scale_color_manual(values = gene_colours) +
    scale_x_discrete(labels = c("F" = "Female", "M" = "Male")) +
    labs(
      x     = "Sex",
      y     = expression("Mean log"[2]~"Normalised Expression"),
      title = tiss
    ) +
    theme_classic(base_size = 14) +
    theme(
      legend.position = "none",
      plot.title      = element_text(face = "bold", hjust = 0.5),
      plot.margin     = margin(5, 70, 5, 5)
    ) +
    coord_cartesian(clip = "off")
}

# =============================================================================
# GENERATE LINE PLOTS
# =============================================================================

tissues    <- c("Anal Fin", "Muscle", "Gonad")
line_plots <- lapply(tissues, make_line_plot)
names(line_plots) <- tissues

for (tiss in tissues) {
  ggsave(
    paste0("./MIG_TPM/line_", gsub(" ", "_", tiss), ".pdf"),
    line_plots[[tiss]], width = 8, height = 6
  )
}
pdf("./MIG_TPM/all_lines.pdf", width = 8, height = 6)
for (p in line_plots) print(p)
dev.off()

# =============================================================================
# CORRELATION PLOTS (oca2 vs lncRNA vs miR10)
# =============================================================================

fc_corr <- fc_long[fc_long$ensembl %in% c("ENSACLG00000024190", "ENSACLG00000029630"),
                   c("sample", "ensembl", "norm_count")]
fc_corr      <- reshape(fc_corr, idvar = "sample", timevar = "ensembl", direction = "wide")
colnames(fc_corr) <- gsub("norm_count\\.", "", colnames(fc_corr))
colnames(fc_corr)[colnames(fc_corr) == "ENSACLG00000024190"] <- "oca2"
colnames(fc_corr)[colnames(fc_corr) == "ENSACLG00000029630"] <- "lncRNA"

miR10_corr        <- miR10_long[, c("sample", "norm_count")]
colnames(miR10_corr)[2] <- "miR10"

corr_data <- merge(fc_corr, miR10_corr, by = "sample")
corr_data <- merge(corr_data, meta_fc[, c("sample", "tissue", "sex")], by = "sample")
corr_data$tissue <- factor(as.character(corr_data$tissue), levels = c("Anal Fin", "Muscle", "Gonad"))
corr_data$sex    <- factor(as.character(corr_data$sex),    levels = c("F", "M"))

cat("Correlation samples:", nrow(corr_data), "\n")

make_corr_plot <- function(df, x_col, y_col, x_lab, y_lab) {
  r_all     <- round(cor(df[[x_col]], df[[y_col]], use = "complete.obs"), 3)
  p_all     <- cor.test(df[[x_col]], df[[y_col]])$p.value
  p_all_lab <- ifelse(p_all < 0.001, "p < 0.001", paste0("p = ", round(p_all, 3)))
  
  r_labels <- do.call(rbind, lapply(
    split(df, as.character(df$tissue)), function(sub) {
      r <- round(cor(sub[[x_col]], sub[[y_col]], use = "complete.obs"), 3)
      p <- tryCatch(cor.test(sub[[x_col]], sub[[y_col]])$p.value, error = \(e) NA_real_)
      data.frame(
        tissue = as.character(sub$tissue[1]),
        x = min(sub[[x_col]], na.rm = TRUE),
        y = max(sub[[y_col]], na.rm = TRUE),
        label = paste0("r = ", r, "\n",
                       ifelse(is.na(p), "p = NA",
                              ifelse(p < 0.001, "p < 0.001",
                                     paste0("p = ", round(p, 3))))),
        stringsAsFactors = FALSE
      )
    }
  ))
  r_labels$tissue <- factor(r_labels$tissue, levels = c("Anal Fin", "Muscle", "Gonad"))
  
  ggplot(df, aes(x = .data[[x_col]], y = .data[[y_col]],
                 color = tissue, shape = sex)) +
    geom_smooth(method = "lm", se = TRUE, color = "black",
                linewidth = 0.7, alpha = 0.15) +
    geom_point(size = 3.5, alpha = 0.85) +
    scale_color_manual(values = tissue_colors) +
    geom_text(data = r_labels, aes(x = x, y = y, label = label),
              inherit.aes = FALSE, hjust = 0, vjust = 1,
              size = 3.5, fontface = "italic") +
    facet_wrap(~ tissue, scales = "free") +
    labs(x = x_lab, y = y_lab,
         title    = paste(x_lab, "vs", y_lab),
         subtitle = paste0("Overall: r = ", r_all, ", ", p_all_lab)) +
    theme_classic(base_size = 13) +
    theme(plot.title    = element_text(face = "bold",  hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5, color = "grey40"),
          strip.text    = element_text(face = "bold"),
          legend.position = "right") +
    guides(color = guide_legend(title = "Tissue"),
           shape = guide_legend(title = "Sex"))
}

p_lnc_oca2   <- make_corr_plot(corr_data, "lncRNA", "oca2",   "lncRNA (log2)", "oca2 (log2)")
p_lnc_miR10  <- make_corr_plot(corr_data, "lncRNA", "miR10",  "lncRNA (log2)", "miR10 (log2)")
p_oca2_miR10 <- make_corr_plot(corr_data, "oca2",   "miR10",  "oca2 (log2)",   "miR10 (log2)")

ggsave("./MIG_TPM/corr_lncRNA_vs_oca2.pdf",   p_lnc_oca2,   width = 12, height = 5)
ggsave("./MIG_TPM/corr_lncRNA_vs_miR10.pdf",  p_lnc_miR10,  width = 12, height = 5)
ggsave("./MIG_TPM/corr_oca2_vs_miR10.pdf",    p_oca2_miR10, width = 12, height = 5)

pdf("./MIG_TPM/all_correlations.pdf", width = 12, height = 5)
print(p_lnc_oca2); print(p_lnc_miR10); print(p_oca2_miR10)
dev.off()


# =============================================================================
# SUMMARY TABLE
# =============================================================================
write.csv(sex_mean_data, "./MIG_TPM/sex_mean_expression_summary.csv", row.names = FALSE)
