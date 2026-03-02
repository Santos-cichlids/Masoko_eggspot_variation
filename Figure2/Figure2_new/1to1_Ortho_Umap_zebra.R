# Setup -------------------------------------------------------------------
#setwd("~/Library/CloudStorage/OneDrive-UniversityofCambridge/Single_cell_Lauren/drive-download-20260211T101947Z-1-001")

# Clear environment
rm(list = ls())

library(Seurat)
library(ggplot2)
library(dplyr)
library(viridis)

# Load data ---------------------------------------------------------------
pigment_cds <- readRDS("zscape_ref_pigment_15k/cds_object.rds")

# Extract data directly
counts_mat <- pigment_cds@assays@data@listData$counts
metadata_df <- as.data.frame(pigment_cds@colData)

# Create Seurat object
seurat_obj <- CreateSeuratObject(counts = counts_mat, meta.data = metadata_df)
seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)

# Add UMAP coordinates
umap_mat <- pigment_cds@int_colData@listData$reducedDims$UMAP
seurat_obj[["umap"]] <- CreateDimReducObject(embeddings = umap_mat, key = "UMAP_")

# Clean up original object
rm(pigment_cds)

# Gene name mapping (based on 1-1 orthology from ENSEMBL) ----------------
gene_name_map <- c(
  "ENSDARG00000061303" = "oca2",
  "ENSDARG00000081044" = "miR125A",
  "ENSDARG00000025953" = "RHOQ",
  "ENSDARG00000069382" = "mxtx1",
  "ENSDARG00000109902" = "TMEM254",
  "ENSDARG00000078417" = "PLXDC1",
  "ENSDARG00000069099" = "TMCO1",
  "ENSDARG00000069100" = "aldh9a1a.1",
  "ENSDARG00000117353" = "lncRNA",
  "ENSDARG00000042018" ="flh2b")

# Genes of interest (only the mapped genes) -------------------------------
gene_list <- names(gene_name_map)

# Find which genes exist in dataset
genes_present <- gene_list[gene_list %in% rownames(seurat_obj)]
cat("Found these genes:", paste(genes_present, collapse = ", "), "\n")
cat("Missing genes:", paste(gene_list[!gene_list %in% rownames(seurat_obj)], collapse = ", "), "\n")

# Create plots directory
dir.create("plots", showWarnings = FALSE)

# Define cell type colors
cell_colors <- c(
  "iridophore" = "grey",
  "melanophore" = "grey10",
  "xanthophore" = "orange"
)

# Plot cell types ---------------------------------------------------------
p1 <- DimPlot(
  seurat_obj, 
  reduction = "umap", 
  group.by = "cell_type_broad",
  pt.size = 0.8, 
  label = TRUE, 
  label.size = 5,
  cols = cell_colors
)
ggsave("plots/Suhaib_pigment_celltypes.pdf", p1, height = 5, width = 6)

# Plot individual genes with GREY ZEROS -----------------------------------
for (gene in genes_present) {
  # Get display name from mapping
  display_name <- gene_name_map[gene]
  
  # Get expression data
  expr_values <- FetchData(seurat_obj, vars = gene)[,1]
  
  # Create a custom color vector: grey for zeros, viridis for non-zeros
  color_data <- data.frame(
    expression = expr_values,
    is_zero = expr_values == 0
  )
  
  p <- FeaturePlot(seurat_obj, features = gene, reduction = "umap", 
                   pt.size = 0.8, order = TRUE) +
    scale_color_gradientn(
      colors = c("grey80", viridis(100, option = "cividis")),
      values = c(0, seq(0.001, 1, length.out = 100)),
      na.value = "grey80"
    ) +
    theme_void() +
    labs(title = display_name)
  
  ggsave(paste0("plots/Suhaib_pigment_", display_name, "_expression.pdf"),
         p, height = 5, width = 6)
  cat("Saved plot for", display_name, "\n")
}

# Plot all genes together with GREY ZEROS ---------------------------------
if (length(genes_present) > 0) {
  p_all <- FeaturePlot(seurat_obj, features = genes_present, reduction = "umap", 
                       pt.size = 0.5, ncol = 3, order = TRUE) & 
    scale_color_gradientn(
      colors = c("grey80", viridis(100, option = "cividis")),
      values = c(0, seq(0.001, 1, length.out = 100)),
      na.value = "grey80"
    )
  
  # Update facet labels with gene names
  for (i in seq_along(genes_present)) {
    p_all[[i]] <- p_all[[i]] + 
      labs(title = gene_name_map[genes_present[i]])
  }
  
  ggsave("plots/Suhaib_pigment_all_genes.pdf", p_all, height = 10, width = 11)
}

# Correlation plots WITHOUT regression lines ------------------------------

# Define gene pairs to compare
gene1_id <- "ENSDARG00000042018"  # lncRNA
gene2_id <- "ENSDARG00000025953"  # oca2

gene1_name <- gene_name_map[gene1_id]
gene2_name <- gene_name_map[gene2_id]

# Check if both genes are present
if (gene1_id %in% genes_present & gene2_id %in% genes_present) {
  
  # Extract expression data
  expr_data <- data.frame(
    lncRNA = FetchData(seurat_obj, vars = gene1_id)[,1],
    oca2 = FetchData(seurat_obj, vars = gene2_id)[,1],
    cell_type = seurat_obj$cell_type_broad
  )
  
  # Remove cells where either gene has zero expression
  expr_data_nonzero <- expr_data %>%
    filter(lncRNA > 0 & oca2 > 0)
  

  # Count by cell type
  print(table(expr_data_nonzero$cell_type))
  
  # Overall correlation (non-zero cells only)
  if (nrow(expr_data_nonzero) > 2) {
    cor_all <- cor.test(expr_data_nonzero$lncRNA, expr_data_nonzero$oca2, method = "pearson")
    print(cor_all)
    
    # Scatter plot - ALL POINTS (with zeros as grey)
    p_scatter_all_with_zeros <- ggplot(expr_data, aes(x = oca2, y = lncRNA)) +
      geom_point(data = expr_data %>% filter(lncRNA == 0 | oca2 == 0), 
                 alpha = 0.3, size = 1.5, color = "grey70") +
      geom_point(data = expr_data_nonzero, aes(color = cell_type),
                 alpha = 0.6, size = 2) +
      scale_color_manual(values = cell_colors) +
      labs(
        title = paste0(gene1_name, " vs ", gene2_name, " - All Points"),
        x = paste0(gene2_name, " Expression"),
        y = paste0(gene1_name, " Expression"),
        color = "Cell Type"
      ) +
      theme_minimal() +
      theme(legend.position = "right")
    
    ggsave("plots/Suhaib_correlation_FHL2B_vs_RHOQ_all_with_zeros.pdf", 
           p_scatter_all_with_zeros, height = 5, width = 7)
    
    # Scatter plot - ZEROS DROPPED (no line)
    p_scatter_all_nonzero <- ggplot(expr_data_nonzero, aes(x = oca2, y = lncRNA, color = cell_type)) +
      geom_point(alpha = 0.6, size = 2) +
      scale_color_manual(values = cell_colors) +
      labs(
        title = paste0(gene1_name, " vs ", gene2_name, " - Non-zero Only"),
        x = paste0(gene2_name, " Expression"),
        y = paste0(gene1_name, " Expression"),
        color = "Cell Type"
      ) +
      theme_minimal() +
      theme(legend.position = "right") +
      annotate("text", x = Inf, y = Inf, 
               label = sprintf("r = %.3f\np = %.2e\nn = %d", 
                               cor_all$estimate, cor_all$p.value, nrow(expr_data_nonzero)),
               hjust = 1.1, vjust = 1.5, size = 4)
    
    ggsave("plots/Suhaib_correlation_FHL2B_vs_RHOQ_nonzero_no_line.pdf", 
           p_scatter_all_nonzero, height = 5, width = 7)
  } else {
    cat("\nNot enough cells for overall correlation\n")
  }
  
  # Correlation by cell type (non-zero cells only)
  cat("\n=== Correlation by Cell Type (non-zero cells) ===\n")
  cor_by_type <- expr_data_nonzero %>%
    group_by(cell_type) %>%
    summarise(
      n_cells = n(),
      correlation = ifelse(n() > 2, cor(lncRNA, oca2, method = "pearson"), NA),
      p_value = ifelse(n() > 2, cor.test(lncRNA, oca2)$p.value, NA),
      .groups = "drop"
    )
  print(cor_by_type)
  
  # Separate scatter plots for EACH cell type (NO REGRESSION LINE)
  for (ct in unique(expr_data_nonzero$cell_type)) {
    subset_data <- expr_data_nonzero %>% filter(cell_type == ct)
    
    if (nrow(subset_data) > 2) {
      
      cor_test <- cor.test(subset_data$lncRNA, subset_data$oca2, method = "pearson")
      
      p_scatter_ct <- ggplot(subset_data, aes(x = oca2, y = lncRNA)) +
        geom_point(alpha = 0.6, size = 2.5, color = cell_colors[ct]) +
        labs(
          title = paste0(gene1_name, " vs ", gene2_name, " - ", ct),
          x = paste0(gene2_name, " Expression"),
          y = paste0(gene1_name, " Expression")
        ) +
        theme_minimal() +
        annotate("text", x = Inf, y = Inf, 
                 label = sprintf("r = %.3f\np = %.2e\nn = %d", 
                                 cor_test$estimate, cor_test$p.value, nrow(subset_data)),
                 hjust = 1.1, vjust = 1.5, size = 5)
      
      ggsave(paste0("plots/Suhaib_correlation_FHL2B_vs_RHOQ_", ct, "_nonzero_no_line.pdf"),
             p_scatter_ct, height = 5, width = 6)
      
      cat("Saved correlation plot for", ct, "\n")
    } else {
      cat("Skipping", ct, "- not enough cells (n =", nrow(subset_data), ") for correlation/plot\n")
    }
  }
  
  cat("\nCorrelation analysis complete!\n")
  
} else {
  cat("\nOne or both genes not found in dataset\n")
}






