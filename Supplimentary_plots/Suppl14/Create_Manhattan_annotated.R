## Combined Manhattan plot of Tajima's D across all 4 populations
## (OVERALL, LITTORAL, MIDDLE, BENTHIC), stacked as facets in one figure.
##
## Thresholds: single empirical 2.5% / 97.5% "candidate" line, computed
## SEPARATELY per population, drawn as one dotted black line on each facet.
##
## Gene annotations: three genes of interest (+5kb flank on each end),
## pulled from the reference GTF by Ensembl gene ID, shaded and labeled
## on every facet.

library(data.table)
library(ggplot2)

## chromosomes present in your files (note: chr21 is absent throughout)
chr_levels <- paste0("chr", c(1:20, 22, 23))

## one entry per population: label used in facet titles + filename glob pattern
## order here also sets the top-to-bottom facet order
pops <- list(
  OVERALL  = "10KB_WINDOW/Overall_10KB_HighSD.*\\.Tajima\\.D$",
  LITTORAL = "10KB_WINDOW/Littoral_10KB_HighSD.*\\.Tajima\\.D$",   # lowercase 'kb'
  MIDDLE   = "10KB_WINDOW/Middle_10KB_HighSD.*\\.Tajima\\.D$",
  BENTHIC  = "10KB_WINDOW/Benth_10KB_HighSD.*\\.Tajima\\.D$"
)

## GTF, used only to look up the three genes of interest below
gtf_file <- "fAstCal1.2.111.gtf.gz"

genes_of_interest <- data.table(
  gene_id   = c("ENSACLG00000024190", "ENSACLG00000029630", "ENSACLG00000017026"),
  gene_name = c("oca2",               "lncRNA",              "TMEM254")
)
flank_bp <- 5000   # +/- 5kb on each end of the gene

## 1. read + combine everything, computing BPcum and thresholds per population --
read_one_pop <- function(pattern, pop_name) {

  files <- list.files(".", pattern = pattern, full.names = TRUE)
  if (length(files) == 0) stop("No files found for pattern: ", pattern)

  d <- rbindlist(lapply(files, fread), use.names = TRUE)
  d <- d[!is.na(TajimaD)]
  d[, CHROM := factor(CHROM, levels = chr_levels)]
  d <- d[!is.na(CHROM)]
  d <- d[order(CHROM, BIN_START)]

  chr_len <- d[, .(len = max(BIN_START)), by = CHROM]
  chr_len[, offset := cumsum(shift(len, fill = 0))]
  d <- merge(d, chr_len[, .(CHROM, offset)], by = "CHROM")
  d[, BPcum := BIN_START + offset]
  d[, Population := pop_name]

  d
}

all_data <- rbindlist(Map(read_one_pop, pops, names(pops)))
all_data[, Population := factor(Population, levels = names(pops))]

## x-axis chromosome label positions (same across populations, since BIN_START/offset match)
axis_df <- all_data[Population == names(pops)[1],
                     .(center = (min(BPcum) + max(BPcum)) / 2), by = CHROM]

## unified per-chromosome offset lookup (BPcum - BIN_START), taken from the
## first population, consistent with the assumption already made for axis_df
offset_lookup <- unique(all_data[Population == names(pops)[1],
                                  .(CHROM, offset = BPcum - BIN_START)])

## 2. population-specific thresholds: single 2.5% / 97.5% dotted black line ------
thresh_df <- all_data[, {
  q <- quantile(TajimaD, probs = c(0.025, 0.975))
  .(type  = c("lower_cand", "upper_cand"),
    value = c(q[["2.5%"]], q[["97.5%"]]),
    label = sprintf(c("2.5%%: %.2f", "97.5%%: %.2f"), c(q[["2.5%"]], q[["97.5%"]])))
}, by = Population]

## 3. look up genes of interest in the GTF, expand by flank_bp on each end -------
gtf <- fread(cmd = sprintf("zcat %s", gtf_file), header = FALSE, sep = "\t", quote = "")
setnames(gtf, c("seqname", "source", "feature", "start", "end",
                "score", "strand", "frame", "attribute"))

gtf_genes <- gtf[feature == "gene"]
gtf_genes[, gene_id := sub('.*gene_id "([^"]+)".*', "\\1", attribute)]

genes_of_interest <- merge(
  genes_of_interest,
  gtf_genes[, .(gene_id, seqname, start, end)],
  by = "gene_id", all.x = TRUE
)

if (anyNA(genes_of_interest$seqname)) {
  warning("Could not find these gene IDs in the GTF: ",
          paste(genes_of_interest$gene_id[is.na(genes_of_interest$seqname)], collapse = ", "))
}

## match the "chr1", "chr2", ... naming style used in the Tajima's D files;
## adjust this line if your GTF already uses "chr"-prefixed names
genes_of_interest[, CHROM := ifelse(grepl("^chr", seqname), seqname, paste0("chr", seqname))]
genes_of_interest[, CHROM := factor(CHROM, levels = chr_levels)]

genes_of_interest[, `:=`(
  flank_start = pmax(1, start - flank_bp),
  flank_end   = end + flank_bp
)]

genes_of_interest <- merge(genes_of_interest, offset_lookup, by = "CHROM", all.x = TRUE)
genes_of_interest[, `:=`(
  xmin = flank_start + offset,
  xmax = flank_end   + offset
)]

## flag points that fall inside any gene region of interest, so they can be
## colored orange regardless of their chromosome's grey shade
all_data[, in_gene_region := FALSE]
for (i in seq_len(nrow(genes_of_interest))) {
  all_data[BPcum >= genes_of_interest$xmin[i] & BPcum <= genes_of_interest$xmax[i],
           in_gene_region := TRUE]
}

## 4. plot -----------------------------------------------------------------------
p <- ggplot() +
  ## background points, alternating grey by chromosome (drawn first)
  geom_point(data = all_data[in_gene_region == FALSE],
             aes(x = BPcum, y = TajimaD, color = CHROM), size = 0.5, show.legend = FALSE) +
  scale_color_manual(values = rep(c("grey30", "grey70"), length.out = length(chr_levels))) +
  ## points inside a gene region of interest, drawn on top in orange
  ## (fixed color, not mapped via aes(), so it doesn't touch the scale above)
  geom_point(data = all_data[in_gene_region == TRUE],
             aes(x = BPcum, y = TajimaD), color = "orange", size = 1.2, show.legend = FALSE) +
  scale_x_continuous(labels = axis_df$CHROM, breaks = axis_df$center) +

  ## gene highlights (recycled across every facet since Population isn't in this data)
  geom_rect(data = genes_of_interest,
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
            inherit.aes = FALSE, fill = "gold", alpha = 0.25) +
  geom_text(data = genes_of_interest,
            aes(x = (xmin + xmax) / 2, y = Inf, label = gene_name),
            inherit.aes = FALSE, vjust = 1.4, size = 3, fontface = "italic") +

  ## single dotted black threshold line (2.5% / 97.5%), per population
  geom_hline(data = thresh_df, aes(yintercept = value),
             color = "black", linetype = "dotted") +
  geom_text(data = thresh_df, aes(x = max(all_data$BPcum), y = value, label = label),
            color = "black", inherit.aes = FALSE, hjust = 1, size = 2.8,
            vjust = ifelse(thresh_df$type == "upper_cand", -0.4, 1.4)) +

  facet_wrap(~ Population, ncol = 1, scales = "free_y") +
  labs(x = "Chromosome", y = "Tajima's D") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 7),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", hjust = 0))

ggsave("tajima_manhattan_all_populations_354HighSD.png", p, width = 14, height = 16, dpi = 300)
ggsave("tajima_manhattan_all_populations_354HighSD.pdf", p, width = 14, height = 16)