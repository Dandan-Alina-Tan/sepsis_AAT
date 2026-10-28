library(TwoSampleMR)
library(ggplot2)
library(tidyverse)
library(dplyr)


files <- list.files(pattern = "\\.tsv$")

# Exposure file
exposure_file <- "UKB-PPP/AAT/summarystats"
exposure_data <- read_exposure_data(
  filename = exposure_file,
  sep = "\t",
  snp_col = "SNP",
  beta_col = "BETA",
  se_col = "SE",
  effect_allele_col = "ALLELE1",
  other_allele_col = "ALLELE0",
  eaf_col = "A1FREQ",
  pval_col = "P",
  samplesize_col = "N"
)

# Initialize results as a data frame
results <- data.frame(
  outcome_file = character(),
  beta = numeric(),
  se = numeric(),
  odd_ratio = numeric(),
  odd_lower = numeric(),
  odd_upper = numeric(),
  p_val = numeric(),
  stringsAsFactors = FALSE
)

# Loop over outcome files
for (file in files) {
  if (file != exposure_file) {
    outcome_data <- read_outcome_data(
      filename = file,
      sep = "\t",
      snp_col = "SNP",
      beta_col = "BETA",
      se_col = "SE",
      effect_allele_col = "ALLELE1",
      other_allele_col = "ALLELE0",
      eaf_col = "A1FREQ",
      pval_col = "P",
      samplesize_col = "N"
    )
    
    harmonized_data <- harmonise_data(
      exposure_dat = exposure_data,
      outcome_dat = outcome_data,
      action = 2
    )
    mr_res <- mr(harmonized_data)
    if (nrow(mr_res) > 0) {
      mr_res <- mr(harmonized_data)
      OR <- generate_odds_ratios(mr_res)
      protein_name <- sub("^modified_([^_]+)_.*$", "\\1", basename(file))
      # Append to results
      results <- rbind(
        results,
        data.frame(
          outcome_file = protein_name,
          beta = OR$b,
          se = mr_res$se,
          odd_ratio = OR$or,
          odd_lower = OR$or_lci95,
          odd_upper = OR$or_uci95,
          p_val = OR$pval,
          stringsAsFactors = FALSE
        )
      )
    }
  }
}
