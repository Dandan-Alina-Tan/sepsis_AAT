library(TwoSampleMR)
library(ggplot2)
library(tidyverse)
library(dplyr)


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


outcome_data_bronchiectasis <- read_outcome_data(
  filename = outcome_path_bronchiectasis,
  sep = "\t",
  snp_col = "MarkerName",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "non_effect_allele",
  eaf_col = "effect_allele_frequency",
  pval_col = "P"
)

outcome_data_pneumonia <- read_outcome_data(
  filename = outcome_path_pneumonia,
  sep = "\t",
  snp_col = "MarkerName",
  beta_col = "Effect",
  se_col = "StdErr",
  effect_allele_col = "Allele1",
  other_allele_col = "Allele2",
  eaf_col = "Freq1",
  pval_col = "Pvalue"
)
 
outcome_data_copd <- read_outcome_data(
  filename = outcome_path_copd,
  sep = "\t",
  snp_col = "MarkerName",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "non_effect_allele",
  eaf_col = "effect_allele_frequency",
  pval_col = "P"
)

outcome_data_UTI <- read_outcome_data(
  filename = outcome_path_UTI,
  sep = "\t",
  snp_col = "Name",
  beta_col = "regenie_beta",
  se_col = "regenie_se",
  effect_allele_col = "effect_allele",
  other_allele_col = "non_effect_allele",
  eaf_col = "effect_allele_frequency",
  pval_col = "Pval"
)

outcome_data_skin_infection <- read_outcome_data(
  filename = outcome_path_skin_infection,
  sep = "\t",
  snp_col = "Name",
  beta_col = "regenie_beta",
  se_col = "regenie_se",
  effect_allele_col = "effect_allele",
  other_allele_col = "non_effect_allele",
  eaf_col = "effect_allele_frequency",
  pval_col = "Pval"
)

outcome_data_T2D <- read_outcome_data(
  filename = outcome_path_T2D,
  sep = "\t",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  eaf_col = "Freq1",
  pval_col = "P"
)

outcome_data_hypertension <- read_outcome_data(
  filename = outcome_path_hypertension,
  sep = "\t",
  snp_col = "Name",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "alt",
  other_allele_col = "ref",
  eaf_col = "eaf",
  pval_col = "Pval"
)

outcome_data_sepsis <- read_outcome_data(
  filename = outcome_path_sepsis,
  sep = "\t",
  snp_col = "MarkerName",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "non_effect_allele",
  eaf_col = "effect_allele_frequency",
  pval_col = "P"
)

outcome_data_cirrhosis <- read_outcome_data(
  filename = outcome_path_cirrhosis,
  sep = "\t",
  snp_col = "MarkerName",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "non_effect_allele",
  eaf_col = "effect_allele_frequency",
  pval_col = "P"
)

traits <- c("sepsis", "bronchiectasis", "pneumonia", "COPD", "urinary_tract_infection",
            "skin_infection_bacterial", "cirrhosis", "T2D", "hypertension")

outcome_list <- list(
  sepsis = outcome_data_sepsis,
  bronchiectasis = outcome_data_bronchiectasis,
  pneumonia = outcome_data_pneumonia,
  COPD = outcome_data_copd,
  urinary_tract_infection = outcome_data_UTI,
  skin_infection_bacterial = outcome_data_skin_infection,
  cirrhosis = outcome_data_cirrhosis,
  T2D = outcome_data_T2D,
  hypertension = outcome_data_hypertension
)

prevalence_dict <- c(
  sepsis = 0.007, # 747 per 100,000
  bronchiectasis = 0.01, 
  pneumonia = 0.012, 
  COPD = 0.10,
  urinary_tract_infection = 0.2,
  skin_infection_bacterial = 0.1, 
  cirrhosis = 0.014,
  T2D = 0.1,   
  hypertension = 0.32
)

# Create an empty list to store results
results_list <- list()

# Loop through each trait
for (trait in traits) {
  message("Running MR for: ", trait)
  
  outcome_data <- outcome_list[[trait]]
  current_prevalence <- prevalence_dict[[trait]]
  
  # Harmonize exposure and outcome
  harmonized_data <- harmonise_data(
    exposure_dat = exposure_data,
    outcome_dat = outcome_data,
    action = 2
  )

  # Run MR
  mr_res <- mr(harmonized_data, method_list = c("mr_ivw", "mr_ivw_fe"))
  
  # Convert to odds ratios
  OR <- generate_odds_ratios(mr_res)
  
  # Calculate r for Steiger test
  harmonized_data$r.exposure <- get_r_from_pn(
    p = harmonized_data$pval.exposure,
    n = harmonized_data$samplesize.exposure
  )
  
  harmonized_data$r.outcome <- get_r_from_lor(
    lor = harmonized_data$beta.outcome,
    af = harmonized_data$eaf.outcome,
    ncase = harmonized_data$ncase.outcome,       
    ncontrol = harmonized_data$ncontrol.outcome,
    prevalence = current_prevalence
  )
  
  # Run MR-Steiger directionality test
  steiger_res <- directionality_test(harmonized_data)
  
  # Merge Steiger results with the OR dataframe
  OR <- merge(OR, steiger_res, by = c("id.exposure", "id.outcome"), all.x = TRUE)

  # Store results
  results_list[[trait]] <- OR
}
# Calculate F-stat for each individual SNP
harmonized_data$F_stat <- (harmonized_data$beta.exposure / harmonized_data$se.exposure)^2

# Calculate the Mean F-statistic across all SNPs used for this trait
mean_F <- mean(harmonized_data$F_stat, na.rm = TRUE)

# Combine all results into one data frame
final_results <- do.call(rbind, lapply(names(results_list), function(x) {
  df <- results_list[[x]]
  df$trait <- x
  df
}))


###########################     PLOT      #########################################
####################################################################################

trait_order <- c("sepsis","bronchiectasis", "pneumonia", "COPD","urinary_tract_infection",
                 "skin_infection_bacterial", "cirrhosis", "T2D", "asthma", "hypertension")

# Convert trait column to factor with this order
final_results$trait <- factor(final_results$trait, levels = rev(trait_order))

final_results_fe <- final_results %>%
  filter(method == "Inverse variance weighted (fixed effects)") %>%
  filter(trait != "asthma") %>%
  mutate(
    trait = recode(trait,
                   "sepsis" = "Sepsis",
                   "bronchiectasis" = "Bronchiectasis",
                   "pneumonia" = "Pneumonia",
                   "COPD" = "COPD",
                   "urinary_tract_infection" = "Urinary tract infection",
                   "skin_infection_bacterial" = "Skin infection (bacterial)",
                   "cirrhosis" = "Cirrhosis",
                   "T2D" = "T2D",
                   "hypertension" = "Hypertension"
    )
  )

final_results_fe <- final_results_fe %>%
  mutate(effect_sign = ifelse(or < 1, "negative", "positive"))

mr_results_plot <- ggplot(final_results_fe, aes(y = trait, x = or, color = effect_sign)) +
  geom_point(position = position_dodge(width = 0.6), size = 3) +
  geom_errorbarh(aes(xmin = or_lci95, xmax = or_uci95),
                 position = position_dodge(width = 0.6), height = 0.2) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray40") +
  scale_x_log10() +
  scale_color_manual(values = c("negative" = "#67a9cf", "positive" = "#ef8a62")) +
  labs(
    x = "Odds Ratio (95% CI)",
    y = "Traits",
    color = "Effect direction"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.y = element_text(size = 14, color = "black"),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11)
  )
