#USE THIS  SCRIPT TO GET THE GENOME WIDE THRESHOLDS
library(dplyr)
library(tidyr)
library(readr)

# Read and combine files
file_list <- list.files(path = "./", 
                        pattern = "HighSEQ_BENTHvsLITT_df_chr.*_xpEHH.tsv", 
                        full.names = TRUE)

xpEHH_data <- lapply(file_list, read_tsv) %>%
  bind_rows() %>%
  as.data.frame()

# ---- FIX COLUMN NAMES ----
colnames(xpEHH_data)[1:8] <- c("CHR", "POSITION", "freq_maj.benth", "freq_min.benth",
                               "freq_maj.litt", "freq_min.litt", "XP_EHH", "LOGPVALUE")

# Filter missing values
xpEHH_data <- xpEHH_data %>%
  filter(!is.na(XP_EHH))

# ---------- Thresholds ----------
upper_975 <- quantile(xpEHH_data$XP_EHH, 0.975, na.rm = TRUE)
lower_025 <- quantile(xpEHH_data$XP_EHH, 0.025, na.rm = TRUE)

upper_995 <- quantile(xpEHH_data$XP_EHH, 0.995, na.rm = TRUE)
lower_005 <- quantile(xpEHH_data$XP_EHH, 0.005, na.rm = TRUE)

# ---------- Genome-wide means ----------
mean_xpEHH <- mean(xpEHH_data$XP_EHH, na.rm = TRUE)

mean_pos <- xpEHH_data %>% 
  filter(XP_EHH > 0) %>% 
  summarise(mean_pos = mean(XP_EHH, na.rm = TRUE)) %>% 
  pull()

mean_neg <- xpEHH_data %>% 
  filter(XP_EHH < 0) %>% 
  summarise(mean_neg = mean(XP_EHH, na.rm = TRUE)) %>% 
  pull()

# ---------- Print results ----------
cat("=== XP-EHH Thresholds ===\n")
cat("Upper 97.5% =", upper_975, "\n")
cat("Lower 2.5%  =", lower_025, "\n\n")
cat("Upper 99.5% =", upper_995, "\n")
cat("Lower 0.5%  =", lower_005, "\n\n")

cat("=== Genome-wide Means ===\n")
cat("Overall mean xpEHH =", mean_xpEHH, "\n")
cat("Mean xpEHH (>0)    =", mean_pos, "\n")
cat("Mean xpEHH (<0)    =", mean_neg, "\n")

#=== XP-EHH Thresholds ===
#Upper 97.5% = 1.965812 
#Lower 2.5%  = -2.072641 

#Upper 99.5% = 2.743514 
#Lower 0.5%  = -3.064567 

#=== Genome-wide Means ===
#Overall mean xpEHH = -1.075955e-17 
#Mean xpEHH (>0)    = 0.7538642 
#Mean xpEHH (<0)    = -0.7689185