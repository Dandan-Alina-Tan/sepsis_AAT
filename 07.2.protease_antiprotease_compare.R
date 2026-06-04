library(ggplot2)
library(tidyverse)
library(dplyr)
library(vroom)
library(SomaDataIO)
library(ggpubr)
library(ggbeeswarm)
library(ggsignif)
library(rstatix)
library(stringr)

serpina1_genotype_Z<- read.table("bqc19_sepsis_genotype_Z_allele.raw", header = TRUE)
serpina1_genotype_S<-read.table("bqc19_sepsis_genotype_S_allele.raw", header = TRUE)

drug_usage<-read.csv("drug_administrations.csv")
drug_affect_serpina1 <- drug_usage %>%
  filter(
    str_detect(
      drug_name, 
      regex("DEXAMETHASONE|methylPREDNISolone|Betamethasone|Hydrocortisone|Prednisone|TOCILIZUMAB|Sarilumab", 
            ignore_case = TRUE)
    )
  )

drug_affect_serpina1 <- drug_affect_serpina1 %>% 
  mutate(
    drug_administration_date = as.Date(drug_administration_date)
  ) %>% 
  mutate(
    DSO_corticosteroids_use = as.numeric(drug_administration_date - symptom_begin)
  ) %>% filter(DSO_corticosteroids_use>0 & DSO_corticosteroids_use<120)

###########################
#AAT
#ELANE
non_infection <- merged_data_prot_gene %>% 
  filter(DSO > 60) %>% 
  group_by(BQCID,rs28929474_T,rs17580_A,Case,study_id) %>% 
  summarize(
    DSO_range = paste(DSO, collapse = ","),
    # Calculate the average for DSO
    DSO = mean(DSO, na.rm = TRUE),
    # Calculate the average for your specific protein column
    serpina1 = mean(seq.3580.25, na.rm = TRUE),
    elane = mean(seq.13671.40, na.rm = TRUE),
    # Drops the grouping so it doesn't mess up future code
    .groups = "drop",
  )


infection_not_sepsis<-merged_data_prot_gene %>% 
  filter(DSO <= 60 & Case == 0) %>% 
  group_by(BQCID,rs28929474_T,rs17580_A,Case,study_id) %>% 
  summarize(
    DSO_range = paste(DSO, collapse = ","),
    # Calculate the average for DSO
    DSO = mean(DSO, na.rm = TRUE),
    # Calculate the average for protein 
    serpina1 = mean(seq.3580.25, na.rm = TRUE),
    elane = mean(seq.13671.40, na.rm = TRUE),
    # Drops the grouping so it doesn't mess up future code
    .groups = "drop",
  )

infection_sepsis<-merged_data_prot_gene %>% 
  filter(DSO <= 60 & Case == 1) %>% 
  group_by(BQCID,rs28929474_T,rs17580_A,Case,study_id) %>% 
  summarize(
    DSO_range = paste(DSO, collapse = ","),
    # Calculate the average for DSO
    DSO = mean(DSO, na.rm = TRUE),
    # Calculate the average for protein 
    serpina1 = mean(seq.3580.25, na.rm = TRUE),
    elane = mean(seq.13671.40, na.rm = TRUE),
    # Drops the grouping so it doesn't mess up future code
    .groups = "drop",
  )
###################plot no drug adjustment############
non_infection_plot<- non_infection %>% mutate(
  Group = "Baseline\n(DSO > 60)",
  Adjustment = "No adjustment"
)
infection_not_sepsis_plot<- infection_not_sepsis %>% mutate(
  Group = "Not Septic\nCOVID",
  Adjustment = "No adjustment"
)
infection_sepsis_plot<- infection_sepsis %>% mutate(
  Group = "Septic COVID",
  Adjustment = "No adjustment"
)

combined_data<-rbind(non_infection_plot,infection_not_sepsis_plot,infection_sepsis_plot)
combined_data<- combined_data %>%
  mutate(
    serpina1_log2 = log2(serpina1),
    elane_log2 = log2(elane),
    ratio_of_log2 =serpina1_log2-elane_log2
  )

wilcox_result_serpina1 <- pairwise.wilcox.test(combined_data$serpina1_log2, 
                                      combined_data$Group,
                                      p.adjust.method = "bonferroni",
                                      paired = FALSE)
wilcox_result_serpina1$p.value

wilcox_result_elane <- pairwise.wilcox.test(combined_data$elane_log2, 
                                               combined_data$Group,
                                               p.adjust.method = "bonferroni",
                                               paired = FALSE)
wilcox_result_elane$p.value

# 1. Reshape the data
long_data <- combined_data %>%
  select(Group, serpina1_log2, elane_log2) %>%
  pivot_longer(
    cols = c(serpina1_log2, elane_log2),
    names_to = "Gene",
    values_to = "Log2_Abundance"
  ) %>%
  mutate(
    Gene = ifelse(Gene == "serpina1_log2", "AAT", "ELANE"),
    Gene = factor(Gene, levels = c("AAT", "ELANE")),
    Group = factor(Group, levels = c("Baseline\n(DSO > 60)", "Not Septic\nCOVID", "Septic COVID"))
  )

combined_plot <- ggplot(long_data, aes(x = Group, y = Log2_Abundance, fill = Group)) +
  geom_quasirandom(shape = 21, size = 1.5, width = 0.2, alpha = 0.3) +
  geom_boxplot(width = 0.6, color = "black", outlier.shape = NA) +
  facet_wrap(~ Gene, strip.position = "bottom") +
  scale_fill_manual(values = c(
    "Baseline\n(DSO > 60)" = "#31a354", 
    "Not Septic\nCOVID"    = "#fec44f", 
    "Septic COVID"        = "#e6550d"  
  )) +
  theme_classic() +
  labs(
    x = NULL, # Completely removes the generic x-axis title
    y = expression(bold(Log[2]~"(Protein Abundance)"))
    
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    axis.text.x = element_text(size = 14, color = "black"),
    axis.text.y = element_text(size = 14, color = "black"),
    axis.title = element_text(size = 16, face = "bold"),
        legend.position = "none",
        strip.text = element_text(size = 14, face = "bold"),
    strip.background = element_blank(),
    strip.placement = "outside" # Ensures the Gene names sit below the x-axis text
  )


long_data_ratio <- combined_data %>%
  # Keep only the columns needed for this specific plot
  select(Group,ratio_of_log2) %>%
  pivot_longer(
    cols = c(ratio_of_log2),
    names_to = "Gene",
    values_to = "Log2_Ratio"
  ) %>%
  mutate(
    Group = factor(Group, levels = c("Baseline\n(DSO > 60)", "Not Septic\nCOVID", "Septic COVID"))
  )

combined_plot_ratio <- ggplot(long_data_ratio, aes(x = Group, y = Log2_Ratio, fill = Group)) +
  
  # 1. The Equilibrium Line
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40", linewidth = 1) +
  
  # 2. Individual Data Points
  geom_quasirandom(shape = 21, size = 1.5, width = 0.2, alpha = 0.3) +
  
  # 3. The Boxplot
  geom_boxplot(width = 0.5, color = "black", outlier.shape = NA, alpha = 0.8) +
 
  scale_fill_manual(values = c(
    "Baseline\n(DSO > 60)" = "#31a354", 
    "Not Septic\nCOVID"    = "#fec44f", 
    "Septic COVID"         = "#e6550d"  
  )) +
  
  theme_classic() +
  labs(
    x = NULL,
    y = expression(bold(Log[2]~"(AAT / ELANE Ratio)"))
  ) +
  theme(
    axis.text.x = element_text(size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 14, face = "bold"),
    legend.position = "none" 
  )

mean_summary <- combined_data %>%
  group_by(Group) %>%
  summarize(
    Mean_Log2_Ratio = mean(ratio_of_log2, na.rm = TRUE),
    SD = sd(ratio_of_log2, na.rm = TRUE),
    n = n()
  )

#######drug adjustment#########
remove_steroid_overlap <- function(prot_data, steroid_data, window_days = 5) {
  
  # 1. Identify patients who had steroids within the specific time window
  overlapping_patients <- prot_data %>%
    separate_longer_delim(cols = DSO_range, delim = ",") %>%
    mutate(true_DSO = as.numeric(DSO_range)) %>%
    inner_join(steroid_data, by = "study_id") %>%
    mutate(days_between = true_DSO - DSO_corticosteroids_use) %>%
    filter(days_between >= 0 & days_between <= window_days) %>%
    distinct(study_id)
  
  # 2. Remove those patients
  clean_data <- prot_data %>%
    anti_join(overlapping_patients, by = "study_id")
  
  num_removed <- nrow(prot_data) - nrow(clean_data)
  message(sprintf("Removed %d patients due to corticosteroid use within %d days of a TRUE blood draw.",
                  num_removed, window_days))
  
  return(clean_data)
}

non_infection_clean <- remove_steroid_overlap(
  prot_data = non_infection,
  steroid_data = drug_affect_serpina1
)

infection_not_sepsis_clean <- remove_steroid_overlap(
  prot_data = infection_not_sepsis,
  steroid_data = drug_affect_serpina1
)

infection_sepsis_clean <- remove_steroid_overlap(
  prot_data = infection_sepsis,
  steroid_data = drug_affect_serpina1
)

########plot with drug adjustment ###########
non_infection_plot_clean<- non_infection_clean %>% mutate(
  Group = "Baseline\n(DSO > 60)",
  Adjustment = "Drug adjustment"

)
infection_not_sepsis_plot_clean<- infection_not_sepsis_clean %>% mutate(
  Group = "Not Septic\nCOVID",
  Adjustment = "Drug adjustment"
)
infection_sepsis_plot_clean<- infection_sepsis_clean %>% mutate(
  Group = "Septic COVID",
  Adjustment = "Drug adjustment"
)

combined_data_clean<-rbind(non_infection_plot_clean,infection_not_sepsis_plot_clean,infection_sepsis_plot_clean)
combined_data_clean<- combined_data_clean %>%
  mutate(
    serpina1_log2 = log2(serpina1),
    elane_log2 = log2(elane),
    ratio_of_log2 =serpina1_log2-elane_log2
  )


long_data_ratio_clean <- combined_data_clean %>%
  # Keep only the columns needed for this specific plot
  select(Group,ratio_of_log2) %>%
  pivot_longer(
    cols = c(ratio_of_log2),
    names_to = "Gene",
    values_to = "Log2_Ratio"
  ) %>%
  mutate(
    Group = factor(Group, levels = c("Baseline\n(DSO > 60)", "Not Septic\nCOVID", "Septic COVID"))
  )

combined_plot_ratio_clean <- ggplot(long_data_ratio_clean, aes(x = Group, y = Log2_Ratio, fill = Group)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40", linewidth = 1) +
  
  geom_quasirandom(shape = 21, size = 1.5, width = 0.2, alpha = 0.3) +
  
  geom_boxplot(width = 0.5, color = "black", outlier.shape = NA, alpha = 0.8) +
  
  scale_fill_manual(values = c(
    "Baseline\n(DSO > 60)" = "#31a354", 
    "Not Septic\nCOVID"    = "#fec44f", 
    "Septic COVID"         = "#e6550d"  
  )) +
  
  theme_classic() +
  labs(
    x = NULL,
    y = expression(bold(Log[2]~"(AAT / ELANE Ratio)"))
  ) +
  theme(
    axis.text.x = element_text(size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 14, face = "bold"),
    legend.position = "none" # X-axis does the job of the legend
  )

mean_summary_clean <- combined_data_clean %>%
  group_by(Group) %>%
  summarize(
    Mean_Log2_Ratio = mean(ratio_of_log2, na.rm = TRUE),
    SD = sd(ratio_of_log2, na.rm = TRUE),
    n = n()
  )


#####################################################################
#################plot trajectory with bins###########################
merged_data_bin <- merged_data_prot_gene %>%
  mutate(
    severe_COVID = ifelse(Case == 1, 1, 0),
    not_severe_COVID = ifelse(Case == 0, 1, 0))

# Bin DSO
merged_data_bin <- merged_data_bin %>%
  mutate(DSO_bin = cut(
    DSO,
    breaks = c(-1, 30, 60, 90, 120, Inf),      
    labels = c("0-30", "31-60", "61-90", "91-120", "120+"),
    right = TRUE
  ))

merged_data_bin <- merged_data_bin %>%
  mutate(
    COVID_status = case_when(
      severe_COVID == 1     ~ "Septic COVID",
      not_severe_COVID == 1 ~ "Not Septic COVID"
    )
  )

serpina1_trajectory <- ggplot(merged_data_bin, 
                              aes(x = DSO_bin, y = log2(seq.3580.25), 
                                  color = COVID_status, group = COVID_status)) +
  
  geom_jitter(width = 0.2, alpha = 0.2, size = 1.5, stroke = 0) +
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar",
               width = 0.2, linewidth = 0.8, alpha = 0.8) +
  stat_summary(fun = mean, geom = "line", linewidth = 1.5) +
  stat_summary(fun = mean, geom = "point", size = 1.5, shape = 21, fill = "white", stroke = 1.5) +
  scale_color_manual(values = c(
    "Not Septic COVID" = "#fec44f", # Yellow
    "Septic COVID"     = "#e6550d"  # Orange
  )) +
  labs(
    x = "Days Since Onset (DSO)", 
    y = "AAT level (BQC19)", 
    color = "Cohort"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    
    # Angles the x-axis text so the week/day bins don't overlap
    axis.text.x = element_text(size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 14, face = "bold"),
    
    # Moves the legend cleanly to the top
    legend.position = "top",
    legend.title = element_text(face = "bold", size = 12),
    legend.text = element_text(size = 11)
  )

bin_statistics <- merged_data_bin %>%
  group_by(DSO_bin) %>%
  wilcox_test(seq.3580.25 ~ COVID_status) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")

###with ratio
merged_data_bin <- merged_data_prot_gene_clinical %>%
  mutate(
    severe_COVID = ifelse(Case == 1, 1, 0),
    not_severe_COVID = ifelse(Case == 0, 1, 0))

# Bin DSO
merged_data_bin <- merged_data_bin %>%
  mutate(DSO_bin = cut(
    DSO,
    breaks = c(-1, 30, 60, 90, 120, Inf),     
    labels = c("0-30", "31-60", "61-90", "91-120", "120+"),
    right = TRUE
  ))

merged_data_bin <- merged_data_bin %>%
  mutate(
    COVID_status = case_when(
      severe_COVID == 1     ~ "Septic COVID",
      not_severe_COVID == 1 ~ "Not Septic COVID"
    ) 
  )

merged_data_bin <-merged_data_bin %>%
  mutate(
    serpina1_log2=log2(seq.3580.25),
    elane_log2=log2(seq.13671.40),
    ratio_of_log2 = log2(seq.3580.25) - log2(seq.13671.40)
  )


ratio_trajectory <- ggplot(merged_data_bin, 
                           aes(x = DSO_bin, y = ratio_of_log2, 
                               color = COVID_status, group = COVID_status)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40", linewidth = 1) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, linewidth = 1) +
  stat_summary(fun = mean, geom = "line", linewidth = 1.2) +
  stat_summary(fun = mean, geom = "point", size = 1.5, shape = 21, fill = "white", stroke = 1.5) +
  scale_color_manual(values = c(
    "Not Septic COVID" = "#fec44f", # Yellow
    "Septic COVID"     = "#e6550d" # Orange
    #"Death"            = "#a50f15"  # Dark Red
  )) +
  theme_classic() +
  labs(
    x = "Days Since Onset (DSO)",
    y = expression(bold(Log[2]~"(AAT / ELANE Ratio)")),
    color = "Cohort"
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 14, face = "bold"),
    
    # Place legend cleanly at the top
    legend.position = "top",
    legend.title = element_text(face = "bold", size = 12),
    legend.text = element_text(size = 12)
  )


merged_data_bin_clean <- remove_steroid_overlap(
  prot_data = merged_data_bin,
  steroid_data = drug_affect_serpina1
)

ratio_trajectory_clean <- ggplot(merged_data_bin_clean, 
                           aes(x = DSO_bin, y = ratio_of_log2, 
                               color = COVID_status, group = COVID_status)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40", linewidth = 1) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, linewidth = 1) +
  stat_summary(fun = mean, geom = "line", linewidth = 1.2) +
  stat_summary(fun = mean, geom = "point", size = 1.5, shape = 21, fill = "white", stroke = 1.5) +
  scale_color_manual(values = c(
  "Not Septic COVID" = "#fec44f", # Yellow
  "Septic COVID"     = "#e6550d" # Orange
  )) +
    theme_classic() +
  labs(
    x = "Days Since Onset (DSO)",
    y = expression(bold(Log[2]~"(AAT / ELANE Ratio)")),
    color = "Cohort"
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 14, face = "bold"),
    
    # Place legend cleanly at the top
    legend.position = "top",
    legend.title = element_text(face = "bold", size = 12),
    legend.text = element_text(size = 12)
  )

bin_statistics_clean <- merged_data_bin_clean %>%
  group_by(DSO_bin) %>%
  wilcox_test(ratio_of_log2 ~ COVID_status) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")

