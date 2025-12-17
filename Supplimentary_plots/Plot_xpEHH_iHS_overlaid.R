#===========================================================#
# Three-panel selection scan with overlaid iHH
#===========================================================#

library(dplyr)
library(rtracklayer)

#-----------------------------------------------------------#
# Input files
#-----------------------------------------------------------#
ehh_benthic_file <- "EHH_Benthic_chr23.tsv"
ehh_littoral_file <- "EHH_Littoral_chr23.tsv"
xpihs_file       <- "HighSEQ_BENTHvsLITT_chr23_IHH_xpEHH_combined.tsv"
gwas_file        <- "/rds/project/rds-8b3VcZwY7rY/projects/cichlid/suhaib/masoko/GWAS/masoko_GWAS_eggSpot_theis/Only_SL_as_Covar/output/Masoko_GWAS_with_covariates_SL_covar.assoc.txt"

#---------------------------#
# Plot region
#---------------------------#
focal_chrom <- "chr23"
win_start <- 12100000
win_end   <- 12300000

###############################TO PLOT TMEM_RHOQ ############################

#-----------------------------------------------------------#
# Input files
#-----------------------------------------------------------#
#ehh_benthic_file <- "EHH_Benthic_chr13.tsv"
#ehh_littoral_file <- "EHH_Littoral_chr13.tsv"
#xpihs_file       <- "HighSEQ_BENTHvsLITT_chr13_IHH_xpEHH_combined.tsv"
#gwas_file  <- "/rds/project/rds-8b3VcZwY7rY/projects/cichlid/suhaib/masoko/GWAS/masoko_HUE/output/Masoko_SUBSET_FILTERED_160_HUE_SL.assoc.txt"

#---------------------------#
# Plot region
#---------------------------#
#focal_chrom <- "chr13"
#win_start <- 4350000
#win_end   <- 4850000
##############################################################################

#-----------------------------------------------------------#
# Helper: read and trim EHH data
#-----------------------------------------------------------#
read_ehh <- function(file, pop) {
  dat <- read.table(file, header = TRUE, sep = "\t")
  # Standardize column names
  names(dat) <- toupper(names(dat))
  
  dat %>%
    filter(CHR == focal_chrom, POSITION >= win_start, POSITION <= win_end) %>%
    mutate(Pop = pop)
}

benthic <- read_ehh(ehh_benthic_file, "Benthic")
littoral <- read_ehh(ehh_littoral_file, "Littoral")

#-----------------------------------------------------------#
# xpEHH data
#-----------------------------------------------------------#
xpi <- read.table(xpihs_file, header = TRUE, sep = "\t")
# Standardize column names
names(xpi) <- tolower(names(xpi))

xpi <- xpi %>%
  transmute(CHR = chr,
            POSITION = position,
            xpEHH_logp = logpvalue) %>%
  filter(CHR == focal_chrom, POSITION >= win_start, POSITION <= win_end)

#-----------------------------------------------------------#
# GWAS data
#-----------------------------------------------------------#
gwas <- read.table(gwas_file, header = TRUE)
# Standardize column names
names(gwas) <- tolower(names(gwas))

gwas <- gwas %>%
  transmute(CHR = chr,
            POSITION = ps,
            Pvalue = p_wald) %>%
  filter(CHR == as.numeric(gsub("chr", "", focal_chrom)), 
         POSITION >= win_start, 
         POSITION <= win_end)

#-----------------------------------------------------------#
# Plot setup: 3 panels
#-----------------------------------------------------------#
pdf("chr23_oca2_OVERLAID_IHH_3panel.pdf", width = 10, height = 9)
par(mfrow = c(3, 1), mar = c(2, 4, 2, 1), oma = c(3, 0, 2, 0))

#-----------------------------------------------------------#
# Panel 1: Overlaid iHH (Benthic and Littoral)
#-----------------------------------------------------------#
# Calculate combined y-axis range
y_range <- range(c(benthic$IHH_MAJ, benthic$IHH_MIN, 
                   littoral$IHH_MAJ, littoral$IHH_MIN), na.rm = TRUE)

plot(benthic$POSITION, benthic$IHH_MAJ, type = "l", col = "steelblue2",
     ylim = y_range, ylab = "iHH", xlab = "", 
     main = "Benthic vs Littoral iHH")

# Add benthic minor
lines(benthic$POSITION, benthic$IHH_MIN, col = "lightblue", lty = 2)

# Add littoral major and minor
lines(littoral$POSITION, littoral$IHH_MAJ, col = "goldenrod2")
lines(littoral$POSITION, littoral$IHH_MIN, col = "lightgoldenrod", lty = 2)

legend("topright", 
       legend = c("Benthic IHH_MAJ", "Benthic IHH_MIN", 
                  "Littoral IHH_MAJ", "Littoral IHH_MIN"),
       col = c("steelblue2", "lightblue", "goldenrod2", "lightgoldenrod"), 
       lty = c(1, 2, 1, 2), 
       bty = "n", cex = 0.7)

#-----------------------------------------------------------#
# Panel 2: xpEHH
#-----------------------------------------------------------#
par(mar = c(2, 4, 2, 1))
plot(xpi$POSITION, xpi$xpEHH_logp, type = "l", col = "purple",
     ylab = "xpEHH -logP", xlab = "", 
     main = "xpEHH (Benthic vs Littoral)")

#-----------------------------------------------------------#
# Panel 3: GWAS
#-----------------------------------------------------------#
par(mar = c(4, 4, 2, 1))
plot(gwas$POSITION, -log10(gwas$Pvalue), pch = 20, col = "black",
     ylab = "-log10(P)", xlab = paste0(focal_chrom, " position (bp)"),
     main = "GWAS")

dev.off()