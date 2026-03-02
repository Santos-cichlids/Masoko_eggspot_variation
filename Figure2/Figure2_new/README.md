This is the revised version of figure 2. Previous Figure 2 will be included within figure 1

1. QC
2. Trimming
3. Annotation and counts
4.  Plotting:
  _Adult dataset_:
    data -
      lncRNA dataset from adult anal fins, gonads and muscles = "TPM_counts_lncRNA.txt"
      small RNA dataset from anal fins, gonads and muscles = "Masoko_NormCounts_smallRNA.txt"
    script -
      use _Plot_violin_corr_line.R_ to plot Violin plots, Line plots and correlations for the Genes of interest.

_Domain & Sex dataset_:
    data -
      ACAL_gene_counts_Extracted.csv
      Domain_expression.csv
    scripts - 
      use _line_plot_domain_divergence.R_ to plot the line plots for genes of interest
      use _heatmap_domain_sex.R_ to plot heatmap for domain and sex at the time of trait initiation and divergence.

  _Single cell dataset_
     data -
       cds_object.rds
     scripts -
       use _1to1_ortho_Umap_zebra.R_ to plot the pigment cell expression for each gene orthologous to our gene of interest.

