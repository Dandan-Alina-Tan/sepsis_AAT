library(TwoSampleMR)
library(ggplot2)
library(tidyverse)
library(dplyr)
library(vroom)

# Outcome data
outcome_data <- read_outcome_data(
  filename = "UKB-PPP/AAT/summarystats,
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

# Exposure files
exposure_files <- list.files(path = ".", pattern = "\\.tsv$", full.names = TRUE)

# Initialize results
results <- data.frame(
  exposure_genes = character(),
  beta = numeric(),
  se = numeric(),
  odds_ratio = numeric(),
  odds_lower = numeric(),
  odds_upper = numeric(),
  p_val = numeric(),
  stringsAsFactors = FALSE
)

for (exposure_file in exposure_files) {
  if(exposure_file!="AAT"){
    # Read exposure data
    exposure_data <- read_exposure_data(
      filename = exposure_file,
      sep = "\t",
      snp_col = "variant_id",
      beta_col = "beta",
      se_col = "standard_error",
      effect_allele_col = "effect_allele",
      other_allele_col = "other_allele",
      eaf_col = "effect_allele_frequency",
      pval_col = "p_value",
      samplesize_col = "n"
    )
    
    # Extract gene name: from file or first row if column exists
    exposure_gene <- tryCatch({
      read.table(exposure_file, header = TRUE, nrows = 1)$gene[1]
    }, error = function(e) basename(exposure_file))
    
    # Harmonise
    harmonized_data <- harmonise_data(
      exposure_dat = exposure_data,
      outcome_dat = outcome_data,
      action = 2
    )
    
    mr_res <- mr(harmonized_data)
    # Run MR if harmonization produced SNPs
    if (nrow(mr_res) > 0) {
      mr_res <- mr(harmonized_data)
      OR <- generate_odds_ratios(mr_res)
      
      # Append results
      results <- rbind(
        results,
        data.frame(
          exposure_genes = exposure_gene,
          methods= OR$method,
          beta = OR$b,
          se = OR$se,
          odds_ratio = OR$or,
          odds_lower = OR$or_lci95,
          odds_upper = OR$or_uci95,
          p_val = OR$pval,
          stringsAsFactors = FALSE
        )
      )
    }
  }

}


)
