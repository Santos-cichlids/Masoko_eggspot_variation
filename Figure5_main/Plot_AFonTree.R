# =============================================================================
# Circular phylogeny — two concentric rings + black branches & tip labels
# Ring 1 (inner): AF colour scale
# Ring 2 (outer): avg_eggspot_theis colour scale
# Branches and tip labels: all black
# Output saved to ./two_circles/
# =============================================================================

library(ape)
library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(scales)
suppressPackageStartupMessages({ library(ggtree) })
suppressPackageStartupMessages({ library(ggtreeExtra) })

# =============================================================================
# STEP 1: Load AF + eggspot data
# =============================================================================

af <- read_table("compiled_AF_SNPs.txt", show_col_types = FALSE)
af$snp_label <- paste0(af$CHROM, ":", af$POS)
af$AF[is.na(af$AF)] <- NA

spots <- read_csv("AF_Malawi_oca2_12210500_12213000.csv", show_col_types = FALSE)

eggspot_df <- spots %>%
  group_by(species_name) %>%
  summarise(eggspot = mean(avg_eggspot_theis, na.rm = TRUE), .groups = "drop")

# =============================================================================
# STEP 2: Remap species names
# =============================================================================

remap <- c(
  "Aulonocara_sp_blue_orange_chilumba_2016"   = "Aulonocara_stuartgranti-Chilumba",
  "Aulonocara_sp_malembo_orange_2016"         = "Aulonocara_malembo-orange",
  "Aulonocara_sp_minutus"                     = "Aulonocara_minutus",
  "Chilotilapia_rhoadesii"                    = "Chilotilapia_rhoadesi",
  "Diplotaxodon_sp_big_eye_black_dorsal_2016" = "Diplotaxodon_bigeye",
  "Fossorochromis_rostratus"                  = "Fossochromis_rostratus",
  "Rhamphochromis_sp_grey_back_Chiofu_2016"   = "Rhamphochromis_greyback",
  "Tropheops_sp_chilumba_2023"                = "Tropheops_chilumba",
  "Tropheops_sp_chiofu_yellow_2016"           = "Tropheops_chiofu-yellow",
  "Tropheops_sp_yellow_gular"                 = "Tropheops_yellow-gular"
)
af$species <- ifelse(af$species %in% names(remap), remap[af$species], af$species)

eggspot_df$species_name <- ifelse(
  eggspot_df$species_name %in% names(remap),
  remap[eggspot_df$species_name],
  eggspot_df$species_name
)
eggspot_df$species_name <- gsub("_", " ", gsub("[-]", " ", eggspot_df$species_name))

# =============================================================================
# STEP 3: Load tree and build tip -> species mapping
# =============================================================================

tree <- ape::read.tree("uninverted_genome.ngsdist.nwk")

extract_species <- function(tip) {
  parts <- strsplit(tip, "__")[[1]]
  if (length(parts) >= 3) {
    sp    <- parts[3]
    words <- strsplit(sp, "_")[[1]]
    paste(words[1], paste(words[-1], collapse = "_"), sep = "_")
  } else { tip }
}

tip_species_map <- data.frame(
  tip_label = tree$tip.label,
  species   = sapply(tree$tip.label, extract_species),
  stringsAsFactors = FALSE
) %>%
  mutate(display_name = gsub("[-_]", " ", species))

matched       <- intersect(unique(af$species), unique(tip_species_map$species))
still_missing <- setdiff(unique(af$species), unique(tip_species_map$species))
cat("Matched:", length(matched), "/", length(unique(af$species)), "AF species\n")
if (length(still_missing) > 0) { cat("Still unmatched:\n"); print(still_missing) }

# =============================================================================
# STEP 4: Output directory
# =============================================================================

out_dir <- "two_circles"
if (!dir.exists(out_dir)) dir.create(out_dir)

# =============================================================================
# STEP 5: Colour scales (rings only — branches and labels stay black)
# =============================================================================

# AF ring: blue -> green gradient
af_ring_colours <- c("#2166AC", "#4393C3", "#92C5DE", "#B8E186", "#66BD63", "darkgreen")
af_ring_values  <- rescale(c(0, 0.20, 0.50, 0.70, 0.85, 1.0))

make_af_fill_scale <- function(bar_width = 8) {
  scale_fill_gradientn(
    name     = "ALT allele frequency",
    colors   = af_ring_colours,
    values   = af_ring_values,
    limits   = c(0, 1),
    na.value = "grey78",
    breaks   = c(0, 0.25, 0.5, 0.75, 1),
    labels   = c("0", "0.25", "0.5", "0.75", "1"),
    guide    = guide_colorbar(
      title.position = "top",
      title.hjust    = 0.5,
      barwidth       = unit(bar_width, "cm"),
      barheight      = unit(0.5, "cm"),
      ticks.colour   = "grey30",
      frame.colour   = "grey30"
    )
  )
}

# Eggspot ring: yellow -> red
make_eggspot_fill_scale <- function(bar_width = 8) {
  scale_fill_gradientn(
    name     = "Avg. eggspot number",
    colors   = c("yellow3", "orange2", "orangered3", "red4"),
    na.value = "grey90",
    guide    = guide_colorbar(
      title.position = "top",
      title.hjust    = 0.5,
      barwidth       = unit(bar_width, "cm"),
      barheight      = unit(0.5, "cm"),
      ticks.colour   = "grey30",
      frame.colour   = "grey30"
    )
  )
}

# Because ggtreeExtra uses the active fill scale, we need two separate
# fill aesthetics. We achieve this with ggnewscale::new_scale_fill()
# between the two geom_fruit() calls.
if (!requireNamespace("ggnewscale", quietly = TRUE)) {
  install.packages("ggnewscale")
}
library(ggnewscale)

# =============================================================================
# STEP 6: Build per-tip annotation for one SNP
# =============================================================================

build_tip_df <- function(snp) {
  d <- af %>% filter(snp_label == snp)
  tip_species_map %>%
    left_join(d %>% select(species, AF), by = "species") %>%
    mutate(
      has_af            = !is.na(AF),
      tip_label_display = ifelse(has_af,
                                 paste0(display_name, "  \u2605"),
                                 display_name)
    ) %>%
    left_join(eggspot_df, by = c("display_name" = "species_name")) %>%
    select(label = tip_label, AF, has_af, tip_label_display, display_name, eggspot)
}

# =============================================================================
# STEP 7A: FULL circular tree — all tips, two rings
# =============================================================================

make_full_tree <- function(snp) {
  
  tip_df <- build_tip_df(snp)
  
  p <- ggtree(tree,
              layout        = "fan",
              open.angle    = 15,
              branch.length = "branch.length",
              size          = 0.35,
              color         = "black") %<+% tip_df +   # <-- branches: black
    
    # ── Ring 1 (inner): AF ──────────────────────────────────────────────────
    geom_fruit(
      geom    = geom_tile,
      mapping = aes(y = label, fill = AF),
      width   = 0.001,
      offset  = 1.2,        # pushes ring well outside the tree
      color   = NA
    ) +
    make_af_fill_scale(bar_width = 9) +
    
    # Switch to a fresh fill scale for ring 2
    new_scale_fill() +
    
    # ── Ring 2 (outer): eggspot ─────────────────────────────────────────────
    geom_fruit(
      geom    = geom_tile,
      mapping = aes(y = label, fill = eggspot),
      width   = 0.001,
      offset  = 0.35,       # gap between ring 1 and ring 2
      color   = NA
    ) +
    make_eggspot_fill_scale(bar_width = 9) +
    
    # ── Tip labels: black, italic ────────────────────────────────────────────
    geom_tiplab(
      aes(label    = tip_label_display,
          size     = has_af,
          fontface = ifelse(has_af, "bold.italic", "italic")),
      color  = "black",
      offset = 0.000001,
      hjust  = 0,
      align  = TRUE
    ) +
    
    scale_size_manual(
      values = c("TRUE" = 3.5, "FALSE" = 2.5),
      guide  = "none"
    ) +
    
    labs(
      title    = snp,
      subtitle = "\u2605 = species with AF data   |   inner ring = ALT allele freq.   |   outer ring = avg. eggspot number"
    ) +
    
    theme_tree2() +
    theme(
      plot.title       = element_text(size = 12, face = "bold",    hjust = 0.5,
                                      margin = margin(b = 3)),
      plot.subtitle    = element_text(size =  8, color = "grey40", hjust = 0.5),
      legend.position  = "bottom",
      legend.direction = "horizontal",
      legend.box       = "horizontal",
      legend.title     = element_text(size = 9, face = "bold"),
      legend.text      = element_text(size = 8),
      plot.margin      = margin(20, 20, 20, 20),
      axis.text        = element_blank(),
      axis.ticks       = element_blank(),
      panel.grid       = element_blank()
    )
  
  p
}

# =============================================================================
# STEP 7B: AF-ONLY pruned circular tree — two rings
# =============================================================================

make_af_only_tree <- function(snp) {
  
  tip_df  <- build_tip_df(snp)
  af_tips <- tip_df %>% filter(has_af) %>% pull(label)
  
  if (length(af_tips) < 4) {
    message("Too few AF tips for pruned tree: ", snp); return(NULL)
  }
  
  pruned   <- ape::drop.tip(tree, setdiff(tree$tip.label, af_tips))
  tip_df_p <- tip_df %>% filter(label %in% pruned$tip.label)
  
  p <- ggtree(pruned,
              layout        = "fan",
              open.angle    = 15,
              branch.length = "branch.length",
              size          = 1.0,
              color         = "black") %<+% tip_df_p +   # <-- branches: black
    
    # ── Ring 1 (inner): AF ──────────────────────────────────────────────────
    geom_fruit(
      geom    = geom_tile,
      mapping = aes(y = label, fill = AF),
      width   = 0.001,
      offset  = 1.5,        # pushes ring well outside the tree
      color   = NA
    ) +
    make_af_fill_scale(bar_width = 10) +
    
    new_scale_fill() +
    
    # ── Ring 2 (outer): eggspot ─────────────────────────────────────────────
    geom_fruit(
      geom    = geom_tile,
      mapping = aes(y = label, fill = eggspot),
      width   = 0.001,
      offset  = 0.40,       # gap between ring 1 and ring 2
      color   = NA
    ) +
    make_eggspot_fill_scale(bar_width = 10) +
    
    # ── Tip labels: black, bold italic ──────────────────────────────────────
    geom_tiplab(
      aes(label = tip_label_display),
      color    = "black",
      size     = 3.5,
      offset   = 0.000001,
      hjust    = 0,
      align    = TRUE,
      fontface = "bold.italic"
    ) +
    
    labs(
      title    = paste0(snp, "   [AF species only]"),
      subtitle = "Pruned to species with AF data   |   inner ring = ALT allele freq.   |   outer ring = avg. eggspot number"
    ) +
    
    theme_tree2() +
    theme(
      plot.title       = element_text(size = 12, face = "bold",    hjust = 0.5),
      plot.subtitle    = element_text(size =  9, color = "grey40", hjust = 0.5),
      legend.position  = "bottom",
      legend.direction = "horizontal",
      legend.box       = "horizontal",
      legend.title     = element_text(size = 10, face = "bold"),
      legend.text      = element_text(size =  9),
      plot.margin      = margin(20, 20, 20, 20),
      axis.text        = element_blank(),
      axis.ticks       = element_blank(),
      panel.grid       = element_blank()
    )
  
  p
}

# =============================================================================
# STEP 8: Generate and save
# =============================================================================

snp_labels <- unique(af$snp_label)
full_plots <- list()
af_plots   <- list()

cat("\nGenerating plots for", length(snp_labels), "SNPs...\n")

for (snp in snp_labels) {
  
  safe <- gsub(":", "_", snp)
  
  p_full <- make_full_tree(snp)
  full_plots[[snp]] <- p_full
  f1 <- file.path(out_dir, paste0("full_tree_", safe, ".pdf"))
  ggsave(f1, plot = p_full, width = 24, height = 24,
         device = "pdf", limitsize = FALSE)
  cat("Saved:", f1, "\n")
  
  p_af <- make_af_only_tree(snp)
  if (!is.null(p_af)) {
    af_plots[[snp]] <- p_af
    f2 <- file.path(out_dir, paste0("AF_only_tree_", safe, ".pdf"))
    ggsave(f2, plot = p_af, width = 16, height = 16,
           device = "pdf", limitsize = FALSE)
    cat("Saved:", f2, "\n")
  }
}

# =============================================================================
# STEP 9: Combined overview PDFs
# =============================================================================

if (requireNamespace("patchwork", quietly = TRUE)) {
  library(patchwork)
  
  make_overview <- function(plot_list, fname, pw, ph, title, subtitle) {
    ncols <- ifelse(length(plot_list) <= 2, 1, 2)
    ov <- patchwork::wrap_plots(plot_list, ncol = ncols, guides = "collect") +
      patchwork::plot_annotation(
        title = title, subtitle = subtitle,
        theme = theme(
          plot.title    = element_text(size = 14, face = "bold", hjust = 0.5),
          plot.subtitle = element_text(size =  9, color = "grey40", hjust = 0.5),
          legend.position = "bottom"
        )
      )
    ggsave(file.path(out_dir, fname), plot = ov,
           width  = ncols * pw,
           height = ceiling(length(plot_list) / ncols) * ph,
           device = "pdf", limitsize = FALSE)
    cat("Saved overview:", fname, "\n")
  }
  
  make_overview(
    full_plots, "overview_full_trees.pdf", 25, 25,
    "Malawi cichlids — two concentric rings (AF + eggspot), black branches",
    "Inner ring: ALT allele frequency   |   Outer ring: avg. eggspot number"
  )
  
  if (length(af_plots) > 0) {
    make_overview(
      af_plots, "overview_AF_only_trees.pdf", 17, 17,
      "Malawi cichlids — AF species only, two concentric rings",
      "Inner ring: ALT allele frequency   |   Outer ring: avg. eggspot number"
    )
  }
}

