# Load required libraries
library(edgeR)
library(ggplot2)
library(ggpubr)

# File paths
metaFile <- "meta.txt"
exprFile <- "aaronAllTissues_count_matrix.txt.gz"
gwasFile <- "Egg_spot_Theis_Top_cand.txt"
gtFile <- "Egg_spot_Theis_Top_cand_GT_v2.txt"
plotDir <- "temp_eggspot_STAT2/"

# Tissues to analyze
tissues <- c("analfin", "eye", "gills", "liver")

# Create output directory if needed
if (!dir.exists(plotDir)) {
  dir.create(plotDir)
}

# Read in data
meta <- read.delim(metaFile, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
gt <- read.delim(gtFile, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
gwas <- read.delim(gwasFile, sep = "\t", header = TRUE, stringsAsFactors = FALSE)

counts <- read.delim(exprFile, sep = "\t", stringsAsFactors = FALSE, row.names = 1)
counts <- as.matrix(counts)
dge <- DGEList(counts = counts)
dge <- calcNormFactors(dge)
cpm <- cpm.DGEList(dge)

# Define ecotype color scheme (matching template)
ecotype_colors <- c("Benthic" = "steelblue2", "Intermediate" = "grey", "Littoral" = "goldenrod2")

# Main loop
for (tissue in tissues) {
  tissColName <- paste0(tissue, "ID_RNAseq")
  tmeta <- meta[!is.na(meta[[tissColName]]), ]
  tcpm <- cpm[, tmeta[[tissColName]]]
  
  for (j in 1:nrow(gwas)) {
    tgwas <- gwas[j, ]
    snp <- paste0(tgwas$chr, ":", tgwas$ps)
    gene <- tgwas$geneID
    tgt <- gt[gt$CHROM.POS == snp, ]
    
    # Ensure sample order matches
    if (!all(tmeta$sampleID %in% colnames(tgt))) next
    tgt <- tgt[, tmeta$sampleID]
    
    snpcpm <- tcpm[gene, ]
    names(snpcpm) <- tmeta$sampleID
    
    # Prepare dataframe
    df <- data.frame(
      cpm = snpcpm,
      gt = factor(as.integer(tgt), levels = c(0, 1, 2)),
      col = tmeta$pca_col,
      sampleID = tmeta$sampleID,
      stringsAsFactors = FALSE
    )
    
    # Assign ecotype labels based on pca_col (adjust mapping as needed)
    df$Group_Label <- case_when(
      df$col == "steelblue2" ~ "Benthic",
      df$col == "grey" ~ "Intermediate", 
      df$col == "goldenrod2" ~ "Littoral",
      TRUE ~ "Other"
    )
    
    # Rename genotype levels to allele labels
    levels(df$gt) <- c(
      paste0(tgwas$allele0, tgwas$allele0),
      paste0(tgwas$allele0, tgwas$allele1),
      paste0(tgwas$allele1, tgwas$allele1)
    )
    
    # Plot title and output name
    pName <- gene
    if (tgwas$geneName != "") {
      pName <- paste0(pName, " ", tgwas$geneName)
    }
    pName <- paste0(pName, " ", tgwas$chr, ":", tgwas$ps)
    plotName <- paste0(plotDir, tissue, "_", gsub(" |:", "_", pName), ".pdf")
    
    # Generate comparison list for statistics
    comparison_list <- combn(levels(df$gt), 2, simplify = FALSE)
    
    # Generate plot following template style
    cat("Writing", plotName, "\n")
    p <- ggplot(df, aes(x = gt, y = cpm)) +
      
      # Violin plot
      geom_violin(
        width = 0.8,
        alpha = 0.3,
        color = "black",
        linewidth = 0.4,
        trim = FALSE
      ) +
      
      # Boxplot overlay (narrow)
      geom_boxplot(
        width = 0.1,
        alpha = 0.6,
        outlier.shape = NA,
        color = "black",
        linewidth = 0.4
      ) +
      
      # Jittered data points colored by ecotype
      geom_jitter(
        aes(color = Group_Label),
        width = 0.07,
        height = 0,
        size = 2,
        alpha = 0.2,
        stroke = 0.7
      ) +
      
      # Statistical comparisons
      stat_compare_means(
        method = "t.test",
        comparisons = comparison_list,
        label = "p.signif",
        hide.ns = TRUE
      ) +
      
      # Custom color scales for ecotype
      scale_color_manual(
        values = ecotype_colors,
        name = "Ecotype"
      ) +
      
      # Labels and theme
      labs(
        title = pName,
        y = "CPM",
        x = "Genotype"
      ) +
      theme_classic(base_size = 16) +
      theme(
        legend.position = "right",
        plot.title = element_text(hjust = 0.5, size = 12),
        axis.text = element_text(color = "black", size = 10),
        axis.title = element_text(size = 11),
        axis.line = element_line(color = "black", linewidth = 0.8),
        axis.ticks = element_line(color = "black")
      )
    
    ggsave(plotName, p, width = 6, height = 5, dpi = 300)
  }
}
