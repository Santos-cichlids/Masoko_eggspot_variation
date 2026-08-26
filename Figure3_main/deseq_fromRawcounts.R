# =============================================================================
# DESeq2 + Plotting Script [ribodepleted long RNA + small RNA]
# =============================================================================

library(DESeq2)
library(tidyverse)
library(ggrepel)

# =============================================================================
# OUTPUT DIRECTORIES
# =============================================================================

base_dir    <- "./deseq_rawCounts"
deseq_dir   <- file.path(base_dir, "deseq_results")
norm_dir    <- file.path(base_dir, "norm_counts")
violin_dir  <- file.path(base_dir, "violin_plots")
line_dir    <- file.path(base_dir, "line_plots")
corr_dir    <- file.path(base_dir, "correlations")
volcano_dir <- file.path(base_dir, "volcano_plots")

for (d in c(base_dir, deseq_dir, norm_dir, violin_dir, line_dir, corr_dir, volcano_dir))
  dir.create(d, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# GENE MAPPINGS
# =============================================================================

gene_mapping_violin <- data.frame(
  ensembl = c(
    "ENSACLG00000019781", "ENSACLG00000017026", "ENSACLG00000034115",
    "ENSACLG00000026836", "ENSACLG00000026955", "ENSACLG00000019614",
    "ENSACLG00000014805", "ENSACLG00000001747", "ENSACLG00000019599",
    "ENSACLG00000019578", "ENSACLG00000005191", "ENSACLG00000004021",
    "ENSACLG00000024080", "ENSACLG00000024228", "ENSACLG00000024190",
    "ENSACLG00000025445", "ENSACLG00000029630"
  ),
  gene_name = c(
    "gnas", "tmem254a", "rhoq", "aldh9a1", "tmco1", "prelid3b",
    "ntrk3a", "atp8b1", "aurka", "slc32a1", "ptpa", "kcnh2",
    "herc2", "gabrg3", "oca2", "plxdc1", "lncRNA"
  ),
  stringsAsFactors = FALSE
)

gene_mapping_line <- data.frame(
  ensembl   = c(
    "ENSACLG00000017026", "ENSACLG00000024190", "ENSACLG00000034115",
    "ENSACLG00000010391", "ENSACLG00000029630", "ENSACLG00000010387"
  ),
  gene_name = c("tmem254", "oca2", "rhoq", "miR10", "lncRNA", "let7"),
  stringsAsFactors = FALSE
)

gene_mapping_line_core <- data.frame(
  ensembl   = c("ENSACLG00000017026", "ENSACLG00000024190", "ENSACLG00000034115"),
  gene_name = c("tmem254", "oca2", "rhoq"),
  stringsAsFactors = FALSE
)

gene_mapping_line_noncoding <- data.frame(
  ensembl   = c("ENSACLG00000029630", "ENSACLG00000010391", "ENSACLG00000010387"),
  gene_name = c("lncRNA", "miR10", "let7"),
  stringsAsFactors = FALSE
)

gene_name_lookup <- rbind(
  gene_mapping_violin,
  data.frame(
    ensembl   = c("ENSACLG00000010391", "ENSACLG00000010387"),
    gene_name = c("miR10", "let7"),
    stringsAsFactors = FALSE
  )
)

gene_colours <- c(
  tmem254 = "#1B9E77", oca2   = "#E7298A", rhoq  = "#E6AB02",
  miR10   = "#7570B3", lncRNA = "#A6761D", let7  = "#66A61E"
)
sex_colors    <- c("F" = "#2E8B57", "M" = "yellow3")
tissue_colors <- c("Anal Fin" = "#D55E00", "Muscle" = "#0072B2", "Gonad" = "#009E73")

# =============================================================================
# HELPERS
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

pval_to_stars <- function(p) {
  if (is.na(p))  return(NA_character_)
  if (p < 0.001) return("***")
  if (p < 0.01)  return("**")
  if (p < 0.05)  return("*")
  "ns"
}

prefilter_any_tissue <- function(count_mat, meta,
                                 min_counts = 10, min_samples = 3) {
  passes <- vapply(levels(meta$tissue), function(tiss) {
    samps <- meta$sample[meta$tissue == tiss]
    rowSums(count_mat[, samps, drop = FALSE] >= min_counts) >= min_samples
  }, logical(nrow(count_mat)))
  keep <- rowSums(passes) >= 1
  cat("  Features before pre-filter:", nrow(count_mat), "\n")
  cat("  Features after  pre-filter:", sum(keep), "\n")
  count_mat[keep, , drop = FALSE]
}

# =============================================================================
# LOAD RAW COUNTS
# =============================================================================

raw_long  <- read.table("longRNA_Raw.txt",  header = TRUE, sep = "\t", check.names = FALSE)
raw_small <- read.table("SmallRNA_RAW.txt", header = TRUE, sep = "\t", check.names = FALSE)
colnames(raw_long)[1]  <- "ensembl"
colnames(raw_small)[1] <- "ensembl"

shared_samples <- intersect(colnames(raw_long)[-1], colnames(raw_small)[-1])

count_matrix_longRNA            <- round(as.matrix(raw_long[,  shared_samples]))
rownames(count_matrix_longRNA)  <- raw_long$ensembl

count_matrix_smallRNA           <- round(as.matrix(raw_small[, shared_samples]))
rownames(count_matrix_smallRNA) <- raw_small$ensembl

# =============================================================================
# SAMPLE METADATA
# =============================================================================

meta           <- parse_sample_meta(shared_samples)
rownames(meta) <- meta$sample
meta$tissue    <- factor(meta$tissue, levels = c("Anal Fin", "Muscle", "Gonad"))
meta$sex       <- factor(meta$sex,    levels = c("F", "M"))

tissue_contrasts <- list(
  c("tissue", "Muscle", "Anal Fin"),
  c("tissue", "Gonad",  "Anal Fin"),
  c("tissue", "Gonad",  "Muscle")
)

# =============================================================================
# SHARED DE PIPELINE
# Returns named list of per-tissue sex DE result data frames.
# =============================================================================

run_deseq_pipeline <- function(count_mat, meta, tissue_contrasts,
                               gene_name_lookup, deseq_dir, prefix,
                               apply_prefilter = FALSE) {
  
  if (apply_prefilter) {
    cat("\nPre-filtering", prefix, "matrix:\n")
    count_mat_de <- prefilter_any_tissue(count_mat, meta)
  } else {
    count_mat_de <- count_mat
  }
  
  # Tissue DE
  dds_tissue <- DESeqDataSetFromMatrix(
    countData = count_mat_de[, meta$sample],
    colData   = meta,
    design    = ~ sex + tissue
  )
  dds_tissue <- dds_tissue[rowSums(counts(dds_tissue) >= 10) >= 3, ]
  dds_tissue <- DESeq(dds_tissue)
  
  for (ct in tissue_contrasts) {
    tag    <- paste0(gsub(" ", "_", ct[2]), "_vs_", gsub(" ", "_", ct[3]))
    res_df <- as.data.frame(results(dds_tissue, contrast = ct, alpha = 0.05))
    res_df$ensembl <- rownames(res_df)
    res_df <- merge(res_df, gene_name_lookup, by = "ensembl", all.x = TRUE)
    write.csv(res_df[order(res_df$padj, na.last = TRUE), ],
              file.path(deseq_dir, paste0(prefix, "_tissue_", tag, ".csv")),
              row.names = FALSE)
  }
  
  # Per-tissue sex DE — results returned for line plots and volcano plots
  sex_results <- list()
  for (tiss in levels(meta$tissue)) {
    samps    <- meta$sample[meta$tissue == tiss]
    meta_sub <- droplevels(meta[samps, ])
    
    dds_sex <- DESeqDataSetFromMatrix(
      countData = count_mat_de[, samps, drop = FALSE],
      colData   = meta_sub,
      design    = ~ sex
    )
    dds_sex <- dds_sex[rowSums(counts(dds_sex) >= 10) >= 3, ]
    dds_sex <- DESeq(dds_sex)
    
    res_df <- as.data.frame(results(dds_sex, contrast = c("sex", "M", "F"), alpha = 0.05))
    res_df$ensembl <- rownames(res_df)
    res_df <- merge(res_df, gene_name_lookup, by = "ensembl", all.x = TRUE)
    res_df <- res_df[order(res_df$padj, na.last = TRUE), ]
    
    write.csv(res_df,
              file.path(deseq_dir, paste0(prefix, "_sex_", gsub(" ", "_", tiss), ".csv")),
              row.names = FALSE)
    
    sex_results[[tiss]] <- res_df
  }
  
  invisible(sex_results)
}

# =============================================================================
# PIPELINE 1 — LONG RNA-SEQ
# =============================================================================

sex_results_longRNA <- run_deseq_pipeline(
  count_matrix_longRNA, meta, tissue_contrasts,
  gene_name_lookup, deseq_dir, prefix = "longRNA",
  apply_prefilter = FALSE
)

dds_plot_long       <- DESeqDataSetFromMatrix(
  countData = count_matrix_longRNA[, meta$sample],
  colData   = meta, design = ~ sex + tissue
)
dds_plot_long       <- estimateSizeFactors(dds_plot_long)
norm_counts_longRNA <- log2(counts(dds_plot_long, normalized = TRUE) + 1)
write.csv(as.data.frame(norm_counts_longRNA),
          file.path(norm_dir, "log2_normalised_counts_longRNA.csv"))

# =============================================================================
# PIPELINE 2 — SMALL RNA-SEQ
# =============================================================================

sex_results_smallRNA <- run_deseq_pipeline(
  count_matrix_smallRNA, meta, tissue_contrasts,
  gene_name_lookup, deseq_dir, prefix = "smallRNA",
  apply_prefilter = TRUE
)

dds_plot_small        <- DESeqDataSetFromMatrix(
  countData = count_matrix_smallRNA[, meta$sample],
  colData   = meta, design = ~ sex + tissue
)
dds_plot_small        <- estimateSizeFactors(dds_plot_small)
norm_counts_smallRNA  <- log2(counts(dds_plot_small, normalized = TRUE) + 1)
write.csv(as.data.frame(norm_counts_smallRNA),
          file.path(norm_dir, "log2_normalised_counts_smallRNA.csv"))

# =============================================================================
# MERGE NORMALISED COUNTS — FOR PLOTTING ONLY
# =============================================================================

norm_df_long          <- as.data.frame(norm_counts_longRNA)
norm_df_long$ensembl  <- rownames(norm_df_long)

norm_df_small         <- as.data.frame(norm_counts_smallRNA)
norm_df_small$ensembl <- rownames(norm_df_small)

write.csv(rbind(norm_df_long, norm_df_small),
          file.path(norm_dir, "log2_normalised_counts_merged_plotOnly.csv"),
          row.names = FALSE)

# =============================================================================
# BUILD TIDY LONG DATA FOR PLOTTING
# =============================================================================

make_gene_long <- function(norm_df, ensembl_ids, gene_map, meta_df) {
  norm_df %>%
    dplyr::filter(ensembl %in% ensembl_ids) %>%
    tidyr::pivot_longer(-ensembl, names_to = "sample", values_to = "norm_count") %>%
    dplyr::left_join(gene_map, by = "ensembl") %>%
    dplyr::left_join(meta_df[, c("sample", "tissue", "sex")], by = "sample") %>%
    tidyr::drop_na(tissue) %>%
    dplyr::mutate(
      tissue = factor(as.character(tissue), levels = c("Anal Fin", "Muscle", "Gonad")),
      sex    = factor(as.character(sex),    levels = c("F", "M"))
    ) %>%
    dplyr::ungroup() %>%
    as.data.frame()
}

fc_long <- make_gene_long(norm_df_long, gene_mapping_violin$ensembl, gene_mapping_violin, meta)

miR10_long <- make_gene_long(
  norm_df_small, "ENSACLG00000010391",
  data.frame(ensembl = "ENSACLG00000010391", gene_name = "miR10", stringsAsFactors = FALSE),
  meta
)
let7_long <- make_gene_long(
  norm_df_small, "ENSACLG00000010387",
  data.frame(ensembl = "ENSACLG00000010387", gene_name = "let7", stringsAsFactors = FALSE),
  meta
)

# =============================================================================
# DESEQ2-BASED SIGNIFICANCE LOOKUP FOR PLOTS
#
# For each gene of interest x tissue, extract the padj from the appropriate
# per-tissue sex DE result (longRNA or smallRNA pipeline). This is used for
# both violin sex brackets and line plot annotations.
# =============================================================================

# Combine all per-tissue sex results into one tidy table
# longRNA results cover all violin + line coding + lncRNA genes
# smallRNA results cover miR10 and let7

compile_deseq_sex_sig <- function(sex_results, gene_map) {
  do.call(rbind, lapply(names(sex_results), function(tiss) {
    res <- sex_results[[tiss]]
    res_sub <- res[res$ensembl %in% gene_map$ensembl, ]
    if (nrow(res_sub) == 0) return(NULL)
    data.frame(
      ensembl      = res_sub$ensembl,
      tissue       = tiss,
      padj         = res_sub$padj,
      significance = sapply(res_sub$padj, pval_to_stars),
      stringsAsFactors = FALSE
    )
  }))
}

# Long RNA genes (violin + coding line + lncRNA)
deseq_sex_sig_long <- compile_deseq_sex_sig(
  sex_results_longRNA,
  gene_name_lookup  # covers all long RNA genes of interest
)

# Small RNA genes (miR10, let7)
deseq_sex_sig_small <- compile_deseq_sex_sig(
  sex_results_smallRNA,
  gene_mapping_line_noncoding[gene_mapping_line_noncoding$gene_name %in% c("miR10", "let7"), ]
)

# Single unified table with gene names attached
deseq_sex_sig_all <- rbind(
  merge(deseq_sex_sig_long,  gene_name_lookup,          by = "ensembl", all.x = TRUE),
  merge(deseq_sex_sig_small, gene_mapping_line_noncoding, by = "ensembl", all.x = TRUE)
)
deseq_sex_sig_all$tissue <- factor(deseq_sex_sig_all$tissue,
                                   levels = c("Anal Fin", "Muscle", "Gonad"))

# =============================================================================
# WILCOXON SIGNIFICANCE — tissue comparisons for violin brackets only
# (still exploratory; DESeq2 tissue contrasts are in the CSV outputs)
# =============================================================================

tissue_sig_df <- do.call(rbind, lapply(
  split(fc_long, fc_long$ensembl), function(sub) {
    wt <- function(a, b) tryCatch(
      wilcox.test(sub$norm_count[sub$tissue == a],
                  sub$norm_count[sub$tissue == b])$p.value, error = \(e) NA_real_)
    data.frame(
      ensembl  = sub$ensembl[1],
      sig_AF_M = pval_to_stars(wt("Anal Fin", "Muscle")),
      sig_AF_G = pval_to_stars(wt("Anal Fin", "Gonad")),
      sig_M_G  = pval_to_stars(wt("Muscle",   "Gonad")),
      stringsAsFactors = FALSE
    )
  }
))

make_mirna_tissue_sig <- function(long_df) {
  wt <- function(a, b) tryCatch(
    wilcox.test(long_df$norm_count[long_df$tissue == a],
                long_df$norm_count[long_df$tissue == b])$p.value, error = \(e) NA_real_)
  data.frame(
    sig_AF_M = pval_to_stars(wt("Anal Fin", "Muscle")),
    sig_AF_G = pval_to_stars(wt("Anal Fin", "Gonad")),
    sig_M_G  = pval_to_stars(wt("Muscle",   "Gonad")),
    stringsAsFactors = FALSE
  )
}

miR10_ts_sig <- make_mirna_tissue_sig(miR10_long)
let7_ts_sig  <- make_mirna_tissue_sig(let7_long)

plot_data        <- merge(fc_long, tissue_sig_df, by = "ensembl", all.x = TRUE)
plot_data$tissue <- factor(as.character(plot_data$tissue), levels = c("Anal Fin", "Muscle", "Gonad"))
plot_data$sex    <- factor(as.character(plot_data$sex),    levels = c("F", "M"))

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
# VIOLIN PLOT FUNCTION
# Sex bracket significance comes from DESeq2 padj.
# Tissue bracket significance comes from Wilcoxon (exploratory).
# =============================================================================

make_violin_plot <- function(df, gene_label, ensembl_id,
                             y_label = "Normalised Counts (log2)") {
  
  sig_AF_M <- as.character(df$sig_AF_M[1])
  sig_AF_G <- as.character(df$sig_AF_G[1])
  sig_M_G  <- as.character(df$sig_M_G[1])
  
  get_sex_sig <- function(tis) {
    row <- deseq_sex_sig_all[
      deseq_sex_sig_all$ensembl == ensembl_id &
        deseq_sex_sig_all$tissue  == tis, ]
    if (nrow(row) == 0 || is.na(row$significance[1])) return("ns")
    as.character(row$significance[1])
  }
  
  df$tissue <- factor(as.character(df$tissue), levels = c("Anal Fin", "Muscle", "Gonad"))
  df$sex    <- factor(as.character(df$sex),    levels = c("F", "M"))
  y_pos     <- max(df$norm_count, na.rm = TRUE) + 0.5
  
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
    guides(fill = guide_legend(title = "Sex"), color = guide_legend(title = "Sex"))
  
  add_bracket <- function(p, x1, x2, y, label) {
    p +
      annotate("segment", x = x1, xend = x2,   y = y,        yend = y,        linewidth = 0.8) +
      annotate("segment", x = x1, xend = x1,   y = y - 0.08, yend = y,        linewidth = 0.8) +
      annotate("segment", x = x2, xend = x2,   y = y - 0.08, yend = y,        linewidth = 0.8) +
      annotate("text", x = (x1 + x2) / 2, y = y + 0.15, label = label,
               color = "black", size = 4, fontface = "bold")
  }
  
  if (!is.na(sig_AF_M) && sig_AF_M != "ns") { p <- add_bracket(p, 1, 2, y_pos, sig_AF_M); y_pos <- y_pos + 0.6 }
  if (!is.na(sig_AF_G) && sig_AF_G != "ns") { p <- add_bracket(p, 1, 3, y_pos, sig_AF_G); y_pos <- y_pos + 0.6 }
  if (!is.na(sig_M_G)  && sig_M_G  != "ns") { p <- add_bracket(p, 2, 3, y_pos, sig_M_G)  }
  
  for (tis in c("Anal Fin", "Muscle", "Gonad")) {
    xi   <- which(c("Anal Fin", "Muscle", "Gonad") == tis)
    sig  <- get_sex_sig(tis)
    ymax <- max(df$norm_count[as.character(df$tissue) == tis], na.rm = TRUE)
    if (!is.na(sig) && sig != "ns")
      p <- p +
      annotate("segment", x = xi - 0.2, xend = xi + 0.2,
               y = ymax + 0.15, yend = ymax + 0.15, color = "grey30", linewidth = 0.5) +
      annotate("text", x = xi, y = ymax + 0.30, label = sig,
               color = "grey30", size = 4, fontface = "bold")
  }
  p
}

attach_ts_sig <- function(long_df, ts_sig) {
  long_df$sig_AF_M <- ts_sig$sig_AF_M
  long_df$sig_AF_G <- ts_sig$sig_AF_G
  long_df$sig_M_G  <- ts_sig$sig_M_G
  long_df
}

# =============================================================================
# GENERATE VIOLIN PLOTS
# =============================================================================

plot_list    <- base::split(plot_data, plot_data$gene_name)
violin_plots <- lapply(names(plot_list), function(gname) {
  df  <- plot_list[[gname]]
  eid <- as.character(df$ensembl[1])
  make_violin_plot(df, gname, eid)
})
names(violin_plots) <- names(plot_list)

if (nrow(miR10_long) > 0)
  violin_plots[["miR10"]] <- make_violin_plot(
    attach_ts_sig(miR10_long, miR10_ts_sig), "miR10", "ENSACLG00000010391"
  )

if (nrow(let7_long) > 0)
  violin_plots[["let7"]] <- make_violin_plot(
    attach_ts_sig(let7_long, let7_ts_sig), "let7", "ENSACLG00000010387"
  )

for (gname in names(violin_plots))
  ggsave(file.path(violin_dir, paste0("violin_", gsub("[^A-Za-z0-9_]", "_", gname), ".pdf")),
         violin_plots[[gname]], width = 10, height = 7)

pdf(file.path(violin_dir, "all_violins.pdf"), width = 10, height = 7)
for (p in violin_plots) print(p)
dev.off()

# =============================================================================
# BUILD LONG DATA FOR LINE PLOTS
# =============================================================================

line_fc_long <- fc_long[fc_long$ensembl %in%
                          gene_mapping_line$ensembl[!gene_mapping_line$gene_name %in% c("miR10","let7")], ]
line_fc_long <- merge(line_fc_long[, !colnames(line_fc_long) %in% "gene_name"],
                      gene_mapping_line, by = "ensembl")

line_all <- rbind(
  line_fc_long[, c("gene_name","sex","tissue","norm_count")],
  miR10_long[,  c("gene_name","sex","tissue","norm_count")],
  let7_long[,   c("gene_name","sex","tissue","norm_count")]
)
colnames(line_all)[colnames(line_all) == "norm_count"] <- "log_norm_counts"

line_all_core      <- line_all[line_all$gene_name %in% gene_mapping_line_core$gene_name, ]
line_all_noncoding <- line_all[line_all$gene_name %in% gene_mapping_line_noncoding$gene_name, ]

build_line_tables <- function(dat) {
  out <- do.call(rbind, lapply(
    split(dat, list(dat$gene_name, as.character(dat$sex), as.character(dat$tissue))),
    function(sub) {
      if (nrow(sub) == 0) return(NULL)
      data.frame(
        gene_name       = sub$gene_name[1],
        sex             = sub$sex[1],
        tissue          = sub$tissue[1],
        mean_log_counts = mean(sub$log_norm_counts, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    }
  ))
  out$sex    <- factor(out$sex,    levels = c("F", "M"))
  out$tissue <- factor(out$tissue, levels = c("Anal Fin", "Muscle", "Gonad"))
  out
}

sex_mean_data           <- build_line_tables(line_all)
sex_mean_data_core      <- build_line_tables(line_all_core)
sex_mean_data_noncoding <- build_line_tables(line_all_noncoding)

# =============================================================================
# LINE PLOT SIGNIFICANCE — DESeq2 padj
#
# For each gene x tissue combination shown on the line plot, look up the padj
# from the appropriate pipeline result and attach it as a label positioned
# between the female and male means.
# =============================================================================

build_sig_positions <- function(sex_mean_df, gene_map) {
  # gene_map tells us which ensembl ID corresponds to each gene_name
  do.call(rbind, lapply(
    split(sex_mean_df, list(sex_mean_df$gene_name, as.character(sex_mean_df$tissue))),
    function(sub) {
      if (nrow(sub) == 0) return(NULL)
      gname <- as.character(sub$gene_name[1])
      tiss  <- as.character(sub$tissue[1])
      
      ens   <- gene_map$ensembl[gene_map$gene_name == gname]
      if (length(ens) == 0) return(NULL)
      
      sig_row <- deseq_sex_sig_all[
        deseq_sex_sig_all$ensembl == ens &
          deseq_sex_sig_all$tissue  == tiss, ]
      
      sig_label <- if (nrow(sig_row) > 0 && !is.na(sig_row$significance[1])) {
        as.character(sig_row$significance[1])
      } else "ns"
      
      data.frame(
        gene_name = gname,
        tissue    = tiss,
        mid_x     = 1.5,
        mid_y     = mean(sub$mean_log_counts, na.rm = TRUE),
        sig_label = sig_label,
        stringsAsFactors = FALSE
      )
    }
  )) %>%
    dplyr::mutate(tissue = factor(tissue, levels = c("Anal Fin", "Muscle", "Gonad")))
}

sig_positions           <- build_sig_positions(sex_mean_data,           gene_mapping_line)
sig_positions_core      <- build_sig_positions(sex_mean_data_core,      gene_mapping_line_core)
sig_positions_noncoding <- build_sig_positions(sex_mean_data_noncoding, gene_mapping_line_noncoding)

# =============================================================================
# LINE PLOT FUNCTION
# =============================================================================

make_line_plot <- function(tiss, sex_mean, jitter_dat, sig_pos,
                           gene_colours_use = gene_colours) {
  td  <- sex_mean[  as.character(sex_mean$tissue)   == tiss, ]
  tj  <- jitter_dat[as.character(jitter_dat$tissue) == tiss, ]
  ts  <- sig_pos[   as.character(sig_pos$tissue)    == tiss, ]
  ts  <- ts[!is.na(ts$sig_label) & ts$sig_label != "ns", ]  # only show non-ns
  lbl <- td[as.character(td$sex) == "M", ]
  
  ggplot() +
    geom_jitter(data = tj, aes(sex, log_norm_counts, color = gene_name),
                width = 0.06, size = 2.5, alpha = 0.3, show.legend = FALSE) +
    geom_line(data = td, aes(sex, mean_log_counts, group = gene_name, color = gene_name),
              linewidth = 2.2) +
    geom_point(data = td, aes(sex, mean_log_counts, color = gene_name), size = 5) +
    geom_text(data = ts, aes(mid_x, mid_y, label = sig_label, color = gene_name),
              size = 6, vjust = -0.5, show.legend = FALSE) +
    geom_text_repel(data = lbl, aes(sex, mean_log_counts, label = gene_name, color = gene_name),
                    size = 5, nudge_x = 0.3, hjust = 0, segment.size = 0.4, show.legend = FALSE) +
    scale_color_manual(values = gene_colours_use) +
    scale_x_discrete(labels = c("F" = "Female", "M" = "Male")) +
    labs(x = "Sex", y = expression("Mean log"[2]~"Normalised Expression"), title = tiss) +
    theme_classic(base_size = 15) +
    theme(legend.position = "none",
          plot.title  = element_text(face = "bold", hjust = 0.5),
          plot.margin = margin(5, 70, 5, 5)) +
    coord_cartesian(clip = "off")
}

# =============================================================================
# GENERATE LINE PLOTS
# =============================================================================

tissues                <- c("Anal Fin", "Muscle", "Gonad")
gene_colours_core      <- gene_colours[names(gene_colours) %in% gene_mapping_line_core$gene_name]
gene_colours_noncoding <- gene_colours[names(gene_colours) %in% gene_mapping_line_noncoding$gene_name]

save_line_plots <- function(plots, prefix) {
  for (tiss in tissues)
    ggsave(file.path(line_dir, paste0(prefix, gsub(" ", "_", tiss), ".pdf")),
           plots[[tiss]], width = 8, height = 6)
  pdf(file.path(line_dir, paste0(prefix, "all.pdf")), width = 8, height = 6)
  for (p in plots) print(p)
  dev.off()
}

line_plots <- setNames(lapply(tissues, make_line_plot,
                              sex_mean = sex_mean_data, jitter_dat = line_all,
                              sig_pos = sig_positions), tissues)
save_line_plots(line_plots, "line_")

line_plots_core <- setNames(lapply(tissues, make_line_plot,
                                   sex_mean = sex_mean_data_core, jitter_dat = line_all_core,
                                   sig_pos  = sig_positions_core,
                                   gene_colours_use = gene_colours_core), tissues)
save_line_plots(line_plots_core, "line_core_")

line_plots_noncoding <- setNames(lapply(tissues, make_line_plot,
                                        sex_mean = sex_mean_data_noncoding, jitter_dat = line_all_noncoding,
                                        sig_pos  = sig_positions_noncoding,
                                        gene_colours_use = gene_colours_noncoding), tissues)
save_line_plots(line_plots_noncoding, "line_noncoding_")

# =============================================================================
# SUMMARY TABLE
# =============================================================================

write.csv(sex_mean_data,
          file.path(base_dir, "sex_mean_expression_summary.csv"),
          row.names = FALSE)







#==============================================================================
# POLY-A gene dataset for Astatotilapia calliptera anal fin domain dataset
#==================================HEATMAP=====================================#
# Load required packages
library(DESeq2)
library(tidyverse)
library(pheatmap)

###########################################################################
########################### HELPER FUNCTIONS ##############################
###########################################################################

get_sig_genes <- function(dds, contrast) {
  res <- results(dds, contrast = contrast, alpha = 0.05)
  rownames(res)[which(res$padj < 0.05 & abs(res$log2FoldChange) > 1 & !is.na(res$padj))]
}

###########################################################################
########################### ALL SAMPLES ###################################
###########################################################################

domain_data         <- read.csv("Domain_expression.csv", header = TRUE)
domain_counts_mat   <- domain_data %>% select(-ensembl) %>% as.matrix()
rownames(domain_counts_mat) <- domain_data$ensembl
domain_sample_names <- colnames(domain_counts_mat)

domain_col_data <- data.frame(
  sample  = domain_sample_names,
  group   = case_when(
    grepl("AT$", domain_sample_names) ~ "AT",
    grepl("M$",  domain_sample_names) ~ "M",
    grepl("PT$", domain_sample_names) ~ "PT"
  ),
  initial = sub("_.*", "", domain_sample_names)
) %>% filter(!is.na(group))

rownames(domain_col_data) <- domain_col_data$sample
domain_counts_mat         <- domain_counts_mat[, domain_col_data$sample]
domain_col_data$group     <- factor(domain_col_data$group, levels = c("AT", "M", "PT"))
domain_col_data$initial   <- factor(domain_col_data$initial)

# Pre-filtering: rowSums >= 10 - exactly like theirs
keep              <- rowSums(domain_counts_mat) >= 10
domain_counts_mat <- domain_counts_mat[keep, ]
cat("Genes after pre-filtering:", nrow(domain_counts_mat), "\n")

# DESeq2 - exactly like theirs
dds_domain <- DESeqDataSetFromMatrix(countData = domain_counts_mat,
                                     colData   = domain_col_data,
                                     design    = ~ initial + group)
dds_domain$group <- relevel(dds_domain$group, "AT")
dds_domain <- DESeq(dds_domain)

# Get significant genes - exactly like theirs
res_ATvM  <- results(dds_domain, contrast = c("group", "M",  "AT"))
res_ATvPT <- results(dds_domain, contrast = c("group", "PT", "AT"))
res_MvPT  <- results(dds_domain, contrast = c("group", "M",  "PT"))

# Combine all three results - rbind like theirs then filter
res_ATvM_df  <- as.data.frame(res_ATvM)  %>% rownames_to_column("ensembl")
res_ATvPT_df <- as.data.frame(res_ATvPT) %>% rownames_to_column("ensembl")
res_MvPT_df  <- as.data.frame(res_MvPT)  %>% rownames_to_column("ensembl")

Domain_only_factor <- rbind(res_ATvM_df, res_ATvPT_df, res_MvPT_df)
Domain_only_factor <- subset(Domain_only_factor, padj < 0.05 & abs(log2FoldChange) > 1)

# Remove duplicates - exactly like theirs
duplicate_rows <- duplicated(Domain_only_factor$ensembl)
DESeq_results  <- Domain_only_factor[!duplicate_rows, ]
cat("Significant DE genes (all samples):", nrow(DESeq_results), "\n")

# Get normalised counts for sig genes - DESeq2 median-of-ratios normalised counts (NOT rlog)
df <- counts(dds_domain, normalized = TRUE)[DESeq_results$ensembl, ]

# Z-score - exactly like theirs
df.z <- t(apply(df, 1, scale))
colnames(df.z) <- colnames(df)

# Order columns by domain - exactly like theirs
sample_order <- c(
  rownames(domain_col_data)[domain_col_data$group == "AT"],
  rownames(domain_col_data)[domain_col_data$group == "M"],
  rownames(domain_col_data)[domain_col_data$group == "PT"]
)
df.z_ordered <- df.z[, sample_order]

# Column annotation
col_ann <- data.frame(
  Domain = domain_col_data[sample_order, "group"],
  Sex    = ifelse(grepl("^F", sample_order), "F", "M")
)
rownames(col_ann) <- sample_order

ann_colors <- list(
  Domain = c(AT = "lightgrey", M = "lightgrey", PT = "lightgrey"),
  Sex    = c(F  = "lightyellow", M = "lightgreen")
)

#COLOUR FOR CONSISTENCY WITH other plots
my_colors      <- colorRampPalette(c("#440154", "#35b779", "#fde725"))
color_gradient <- my_colors(100)

# gaps_col - calculate from data to avoid hardcoding
gaps <- c(sum(domain_col_data$group == "AT"),
          sum(domain_col_data$group == "AT") + sum(domain_col_data$group == "M"))

# NOTE: cutree_rows / cutree(k=) now consistently set to 3 throughout,
# matching the reported clustering (3 hierarchical clusters, complete
# linkage, Euclidean distance).
N_CLUSTERS <- 3

# Close any stray open graphics devices before plotting, to avoid a blank
# PDF caused by leftover device state from the interactive screen plot below.
graphics.off()

# Plot to screen
Domain_heatmap <- pheatmap(df.z_ordered,
                           cluster_rows      = TRUE,
                           cluster_cols      = FALSE,
                           cutree_rows       = N_CLUSTERS,
                           fontsize          = 12,
                           gaps_col          = gaps,
                           show_rownames     = FALSE,
                           show_colnames     = FALSE,
                           annotation_col    = col_ann,
                           annotation_colors = ann_colors,
                           main              = "Domain DE Genes (All Samples)",
                           color             = color_gradient)

# Save - fresh device, explicit close before opening new one
graphics.off()
pdf("hier_clust/heatmap_domain_allsamples.pdf", width = 12, height = 14)
pheatmap(df.z_ordered,
         cluster_rows      = TRUE,
         cluster_cols      = FALSE,
         cutree_rows       = N_CLUSTERS,
         fontsize          = 12,
         gaps_col          = gaps,
         show_rownames     = FALSE,
         show_colnames     = FALSE,
         annotation_col    = col_ann,
         annotation_colors = ann_colors,
         main              = "Domain DE Genes (All Samples)",
         color             = color_gradient)
dev.off()

# Save cluster CSV - using hclust on same matrix as pheatmap, k matches pheatmap's cutree_rows
hc_all <- hclust(dist(df.z_ordered), method = "complete")
cluster_csv_all <- data.frame(
  ensembl = rownames(df.z_ordered),
  Cluster = cutree(hc_all, k = N_CLUSTERS)
)
write_csv(cluster_csv_all, "hier_clust/domain_gene_clusters_allsamples.csv")

###########################################################################
########################### MALES ONLY ####################################
###########################################################################

male_samples      <- domain_col_data$sample[!grepl("^F", domain_col_data$sample)]
domain_col_data_M <- domain_col_data[male_samples, ]
domain_counts_M   <- domain_counts_mat[, male_samples]

keep_M          <- rowSums(domain_counts_M) >= 10
domain_counts_M <- domain_counts_M[keep_M, ]

dds_domain_M <- DESeqDataSetFromMatrix(countData = domain_counts_M,
                                       colData   = domain_col_data_M,
                                       design    = ~ initial + group)
dds_domain_M$group <- relevel(dds_domain_M$group, "AT")
dds_domain_M <- DESeq(dds_domain_M)

res_ATvM_M  <- as.data.frame(results(dds_domain_M, contrast = c("group", "M",  "AT"))) %>% rownames_to_column("ensembl")
res_ATvPT_M <- as.data.frame(results(dds_domain_M, contrast = c("group", "PT", "AT"))) %>% rownames_to_column("ensembl")
res_MvPT_M  <- as.data.frame(results(dds_domain_M, contrast = c("group", "M",  "PT"))) %>% rownames_to_column("ensembl")

Domain_only_M <- rbind(res_ATvM_M, res_ATvPT_M, res_MvPT_M)
Domain_only_M <- subset(Domain_only_M, padj < 0.05 & abs(log2FoldChange) > 1)
DESeq_results_M <- Domain_only_M[!duplicated(Domain_only_M$ensembl), ]
cat("Significant DE genes (males only):", nrow(DESeq_results_M), "\n")

df_M <- counts(dds_domain_M, normalized = TRUE)[DESeq_results_M$ensembl, ]

df.z_M <- t(apply(df_M, 1, scale))
colnames(df.z_M) <- colnames(df_M)

sample_order_M <- c(
  rownames(domain_col_data_M)[domain_col_data_M$group == "AT"],
  rownames(domain_col_data_M)[domain_col_data_M$group == "M"],
  rownames(domain_col_data_M)[domain_col_data_M$group == "PT"]
)
df.z_M_ordered <- df.z_M[, sample_order_M]

col_ann_M <- data.frame(Domain = domain_col_data_M[sample_order_M, "group"])
rownames(col_ann_M) <- sample_order_M
ann_colors_M <- list(Domain = c(AT = "#45C4B0", M = "#9AEBA3", PT = "#13678A"))

gaps_M <- c(sum(domain_col_data_M$group == "AT"),
            sum(domain_col_data_M$group == "AT") + sum(domain_col_data_M$group == "M"))

graphics.off()

# Plot to screen
Domain_heatmap_M <- pheatmap(df.z_M_ordered,
                             cluster_rows      = TRUE,
                             cluster_cols      = FALSE,
                             cutree_rows       = N_CLUSTERS,
                             fontsize          = 12,
                             gaps_col          = gaps_M,
                             show_rownames     = FALSE,
                             show_colnames     = FALSE,
                             annotation_col    = col_ann_M,
                             annotation_colors = ann_colors_M,
                             main              = "Domain DE Genes (Males Only)",
                             color             = color_gradient)

# Save
graphics.off()
pdf("hier_clust/heatmap_domain_malesonly.pdf", width = 10, height = 14)
pheatmap(df.z_M_ordered,
         cluster_rows      = TRUE,
         cluster_cols      = FALSE,
         cutree_rows       = N_CLUSTERS,
         fontsize          = 12,
         gaps_col          = gaps_M,
         show_rownames     = FALSE,
         show_colnames     = FALSE,
         annotation_col    = col_ann_M,
         annotation_colors = ann_colors_M,
         main              = "Domain DE Genes (Males Only)",
         color             = color_gradient)
dev.off()

hc_M <- hclust(dist(df.z_M_ordered), method = "complete")
cluster_csv_M <- data.frame(
  ensembl = rownames(df.z_M_ordered),
  Cluster = cutree(hc_M, k = N_CLUSTERS)
)
write_csv(cluster_csv_M, "hier_clust/domain_gene_clusters_malesonly.csv")

##### other CSVs
write_csv(DESeq_results,   "hier_clust/DESeq_sig_genes_allsamples.csv")
write_csv(DESeq_results_M, "hier_clust/DESeq_sig_genes_malesonly.csv")

###########################################################################
########################### SUMMARY COUNTS #################################
###########################################################################

# Per-contrast counts (all samples) - before rbind/dedup
n_ATvM_all  <- sum(res_ATvM_df$padj  < 0.05 & abs(res_ATvM_df$log2FoldChange)  > 1, na.rm = TRUE)
n_ATvPT_all <- sum(res_ATvPT_df$padj < 0.05 & abs(res_ATvPT_df$log2FoldChange) > 1, na.rm = TRUE)
n_MvPT_all  <- sum(res_MvPT_df$padj  < 0.05 & abs(res_MvPT_df$log2FoldChange)  > 1, na.rm = TRUE)

# Per-contrast counts (males only) - before rbind/dedup
n_ATvM_M  <- sum(res_ATvM_M$padj  < 0.05 & abs(res_ATvM_M$log2FoldChange)  > 1, na.rm = TRUE)
n_ATvPT_M <- sum(res_ATvPT_M$padj < 0.05 & abs(res_ATvPT_M$log2FoldChange) > 1, na.rm = TRUE)
n_MvPT_M  <- sum(res_MvPT_M$padj  < 0.05 & abs(res_MvPT_M$log2FoldChange)  > 1, na.rm = TRUE)

###########################################################################
########################### Line Plots dOMAIN #################################
###########################################################################

# Load required packages
library(DESeq2)
library(tidyverse)
library(ggrepel)

# Gene mapping - Set 1
# CHANGE 1: Only keeping tmem254, oca2, rhoq (others commented out)
gene_mapping <- data.frame(
  ensembl = c("ENSACLG00000017026",
              #"ENSACLG00000026836",   # aldh9a1a.a - excluded
              #"ENSACLG00000026955",   # tmco1 - excluded
              "ENSACLG00000024190",
              #"ENSACLG00000025445",   # plxdc1 - excluded
              "ENSACLG00000016986"),
  gene_name = c("tmem254",
                #"aldh9a1a.a",         # excluded
                #"tmco1",              # excluded
                "oca2",
                #"plxdc1",             # excluded
                "rhoq")
)

# Load and prepare sex data - keep ALL genes for proper normalization
count_data <- read.csv("../ACAL_gene_counts_Extracted.csv", header = TRUE)
counts_mat <- count_data %>% select(-ensembl) %>% as.matrix()
rownames(counts_mat) <- count_data$ensembl
sample_names <- colnames(counts_mat)
col_data <- data.frame(sample = sample_names, group = ifelse(grepl("^F", sample_names), "F", "M")) %>%
  filter(group %in% c("F", "M"))
rownames(col_data) <- col_data$sample
counts_mat <- counts_mat[, col_data$sample]
col_data$group <- factor(col_data$group, levels = c("F", "M"))

# DESeq2 analysis for sex - use FULL dataset for proper normalization
# NOTE: design fixed to ~ group only. "sample" was one level per row (just an ID),
# fully confounded with group, which made the model matrix not full rank.
dds_sex <- DESeqDataSetFromMatrix(countData = counts_mat, colData = col_data, design = ~ group)
dds_sex <- DESeq(dds_sex)
normalized_counts_sex <- counts(dds_sex, normalized = TRUE)

# Now subset AFTER normalization
sex_norm_data <- normalized_counts_sex %>%
  as.data.frame() %>%
  rownames_to_column("ensembl") %>%
  filter(ensembl %in% gene_mapping$ensembl) %>%  # Subset here
  pivot_longer(-ensembl, names_to = "sample", values_to = "normalized_count") %>%
  left_join(col_data, by = "sample") %>%
  left_join(gene_mapping, by = "ensembl") %>%
  mutate(log_norm_counts = log2(normalized_count + 1))

# Calculate mean expression for sex
sex_mean_data <- sex_norm_data %>%
  group_by(gene_name, group) %>%
  summarise(mean_log_counts = mean(log_norm_counts, na.rm = TRUE), .groups = "drop") %>%
  mutate(group_num = as.numeric(group))

# Get labels for the end points only
sex_labels <- sex_mean_data %>%
  group_by(gene_name) %>%
  filter(group == "M") %>%
  ungroup()

# CHANGE 1 (colours): Only 3 genes retained; others commented out
gene_colors <- c(
  "tmem254" = "#1B9E77",
  #"aldh9a1a.a" = "#D95F02",   # excluded
  #"tmco1" = "#7570B3",        # excluded
  "oca2" = "#E7298A",
  #"plxdc1" = "#66A61E",       # excluded
  "rhoq" = "#E6AB02"
)

# Create sex plot
# CHANGE 2: Increased line size, point size, and label text size
# CHANGE 3: Added geom_jitter layer (individual points, same colour, translucent) BEFORE mean line/points
sex_plot <- ggplot(sex_mean_data, aes(x = group, y = mean_log_counts, group = gene_name, color = gene_name)) +
  # CHANGE 3: Individual jittered points from raw normalised data (alpha for translucency)
  geom_jitter(data = sex_norm_data,
              aes(x = group, y = log_norm_counts, group = gene_name, color = gene_name),
              width = 0.08, size = 3, alpha = 0.25, inherit.aes = FALSE) +
  # Mean line and points on top
  geom_line(size = 1.8) +                        
  geom_point(size = 5.5) +                        
  geom_text_repel(data = sex_labels,
                  aes(label = gene_name),
                  size = 5.5,                     
                  nudge_x = 0.3,
                  hjust = 0,
                  direction = "y",
                  segment.size = 0.4,
                  segment.color = "grey50") +
  scale_color_manual(values = gene_colors) +
  labs(x = "Sex", y = "Mean Log2 Normalized Counts",
       title = "Gene Expression Across Sex") +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none",
    axis.text = element_text(color = "black"),
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.title = element_text(face = "bold")
  ) +
  coord_cartesian(clip = "off") +
  theme(plot.margin = margin(5, 50, 5, 5))

# Load and prepare domain data - keep ALL genes for proper normalization
domain_data <- read.csv("../Domain_expression.csv", header = TRUE)  # Don't filter yet!
domain_counts_mat <- domain_data %>% select(-ensembl) %>% as.matrix()
rownames(domain_counts_mat) <- domain_data$ensembl
domain_sample_names <- colnames(domain_counts_mat)
domain_col_data <- data.frame(
  sample = domain_sample_names,
  group = case_when(
    grepl("AT$", domain_sample_names) ~ "AT",
    grepl("M$", domain_sample_names) ~ "M",
    grepl("PT$", domain_sample_names) ~ "PT")
) %>% filter(!is.na(group))
rownames(domain_col_data) <- domain_col_data$sample
domain_counts_mat <- domain_counts_mat[, domain_col_data$sample]
domain_col_data$group <- factor(domain_col_data$group, levels = c("AT", "M", "PT"))

# DESeq2 analysis for domains - use FULL dataset for proper normalization
dds_domain <- DESeqDataSetFromMatrix(countData = domain_counts_mat, colData = domain_col_data, design = ~ group)
dds_domain <- DESeq(dds_domain)
normalized_counts_domain <- counts(dds_domain, normalized = TRUE)

# Now subset AFTER normalization
domain_norm_data <- normalized_counts_domain %>%
  as.data.frame() %>%
  rownames_to_column("ensembl") %>%
  filter(ensembl %in% gene_mapping$ensembl) %>%  # Subset here
  pivot_longer(-ensembl, names_to = "sample", values_to = "normalized_count") %>%
  left_join(domain_col_data, by = "sample") %>%
  left_join(gene_mapping, by = "ensembl") %>%
  mutate(log_norm_counts = log2(normalized_count + 1))

# Calculate mean expression for domains
domain_mean_data <- domain_norm_data %>%
  group_by(gene_name, group) %>%
  summarise(mean_log_counts = mean(log_norm_counts, na.rm = TRUE), .groups = "drop") %>%
  mutate(group_num = as.numeric(group))

# Get labels for the end points only
domain_labels <- domain_mean_data %>%
  group_by(gene_name) %>%
  filter(group == "PT") %>%
  ungroup()

# Create domain plot
# CHANGE 2: Increased line size, point size, and label text size
# CHANGE 3: Added geom_jitter layer (individual points, same colour, translucent) BEFORE mean line/points
domain_plot <- ggplot(domain_mean_data, aes(x = group, y = mean_log_counts, group = gene_name, color = gene_name)) +
  # CHANGE 3: Individual jittered points from raw normalised data (alpha for translucency)
  geom_jitter(data = domain_norm_data,
              aes(x = group, y = log_norm_counts, group = gene_name, color = gene_name),
              width = 0.08, size = 3, alpha = 0.25, inherit.aes = FALSE) +
  # Mean line and points on top
  geom_line(size = 1.8) +                        
  geom_point(size = 5.5) +                      
  geom_text_repel(data = domain_labels,
                  aes(label = gene_name),
                  size = 5.5,                     
                  nudge_x = 0.3,
                  hjust = 0,
                  direction = "y",
                  segment.size = 0.4,
                  segment.color = "grey50") +
  scale_color_manual(values = gene_colors) +
  labs(x = "Domain", y = "Mean (Log Normalized Counts)",
       title = "Gene Expression Across Domains") +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none",
    axis.text = element_text(color = "black"),
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.title = element_text(face = "bold")
  ) +
  coord_cartesian(clip = "off") +
  theme(plot.margin = margin(5, 50, 5, 5))

# Save plots
ggsave("sex_mean_expression_all_genes.pdf", plot = sex_plot, width = 8, height = 6, units = "in")
ggsave("domain_mean_expression_all_genes.pdf", plot = domain_plot, width = 10, height = 8, units = "in")

# Display plots
print(sex_plot)
print(domain_plot)


################FOR MALES ONLY
# ── MALE-ONLY DOMAIN EXPRESSION PLOT ─────────────────────────────────────────

# Filter to male samples only (columns starting with "M", excluding "ensembl")
male_sample_names <- domain_sample_names[grepl("^M", domain_sample_names)]

male_counts_mat <- domain_counts_mat[, male_sample_names]

male_col_data <- data.frame(
  sample = male_sample_names,
  group = case_when(
    grepl("AT$", male_sample_names) ~ "AT",
    grepl("M$",  male_sample_names) ~ "M",
    grepl("PT$", male_sample_names) ~ "PT")
) %>% filter(!is.na(group))
rownames(male_col_data) <- male_col_data$sample
male_counts_mat <- male_counts_mat[, male_col_data$sample]
male_col_data$group <- factor(male_col_data$group, levels = c("AT", "M", "PT"))

# DESeq2 normalization on male samples only
dds_male <- DESeqDataSetFromMatrix(countData = male_counts_mat,
                                   colData   = male_col_data,
                                   design    = ~ group)
dds_male <- DESeq(dds_male)
normalized_counts_male <- counts(dds_male, normalized = TRUE)

# Subset to genes of interest AFTER normalization
male_norm_data <- normalized_counts_male %>%
  as.data.frame() %>%
  rownames_to_column("ensembl") %>%
  filter(ensembl %in% gene_mapping$ensembl) %>%
  pivot_longer(-ensembl, names_to = "sample", values_to = "normalized_count") %>%
  left_join(male_col_data, by = "sample") %>%
  left_join(gene_mapping, by = "ensembl") %>%
  mutate(log_norm_counts = log2(normalized_count + 1))

# Mean expression per gene per domain
male_mean_data <- male_norm_data %>%
  group_by(gene_name, group) %>%
  summarise(mean_log_counts = mean(log_norm_counts, na.rm = TRUE), .groups = "drop") %>%
  mutate(group_num = as.numeric(group))

# Labels at PT end point
male_domain_labels <- male_mean_data %>%
  group_by(gene_name) %>%
  filter(group == "PT") %>%
  ungroup()

# Plot
male_domain_plot <- ggplot(male_mean_data,
                           aes(x = group, y = mean_log_counts,
                               group = gene_name, color = gene_name)) +
  geom_jitter(data    = male_norm_data,
              aes(x = group, y = log_norm_counts,
                  group = gene_name, color = gene_name),
              width = 0.08, size = 3, alpha = 0.25, inherit.aes = FALSE) +
  geom_line(size  = 1.8) +
  geom_point(size = 5.5) +
  geom_text_repel(data         = male_domain_labels,
                  aes(label    = gene_name),
                  size         = 5.5,
                  nudge_x      = 0.3,
                  hjust        = 0,
                  direction    = "y",
                  segment.size = 0.4,
                  segment.color = "grey50") +
  scale_color_manual(values = gene_colors) +
  labs(x     = "Domain",
       y     = "Mean (Log Normalized Counts)",
       title = "Gene Expression Across Domains (Males Only)") +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none",
    axis.text       = element_text(color = "black"),
    plot.title      = element_text(face = "bold", hjust = 0.5),
    axis.title      = element_text(face = "bold")
  ) +
  coord_cartesian(clip = "off") +
  theme(plot.margin = margin(5, 50, 5, 5))

# Save and display
ggsave("male_domain_mean_expression.pdf", plot = male_domain_plot,
       width = 10, height = 8, units = "in")
print(male_domain_plot)

###########################################################################
########################### SEX DE GENES ###################################
###########################################################################
# Uses dds_sex, which was already fit above with design = ~ group

dir.create("line_plots", showWarnings = FALSE)

res_sex    <- results(dds_sex, contrast = c("group", "M", "F"), alpha = 0.05)
res_sex_df <- as.data.frame(res_sex) %>% rownames_to_column("ensembl")

DESeq_results_sex <- subset(res_sex_df, padj < 0.05 & abs(log2FoldChange) > 1)

write_csv(DESeq_results_sex, "line_plots/DESeq_sig_genes_sex.csv")

