library(ggplot2)
library(tidyverse)
library(dplyr)
library(vroom)
library(SomaDataIO)
library(ggpubr)
library(ggbeeswarm)
library(ggsignif) 

serpina1_genotype_Z<- read.table("../bqc19_genotype_Z_allele.raw", header = TRUE)
serpina1_genotype_S<-read.table("../bqc19_genotype_S_allele.raw", header = TRUE)

#phenotype (septic COVID or not)
merged_data <- phenotype %>%
  left_join(serpina1_genotype_Z, by = c("individual_id"= "FID")) %>% 
  left_join(serpina1_genotype_S,by = c("individual_id"= "FID"))

#remove row with NA in genotype
merged_data <- merged_data %>%
  filter(!is.na(rs28929474_T)) %>%
  filter(!is.na(rs17580_A))

############################################################
#BQC19 AAT
prot_comb <- prot_comb %>% select(SubjectID,seq.3580.25)
merged_data_prot <- merged_data %>%
  left_join(prot_AAT, by = c("SubjectID"))

#case: whether a indivudal is septic COVID or not
merged_data_prot <- merged_data_prot %>% 
  mutate(
    PI_MM_septic_covid = ifelse(rs28929474_T == 0 & rs17580_A == 0 & Case == 1, 1, 0),
    PI_MZ_septic_covid = ifelse(rs28929474_T == 1 & Case == 1, 1, 0),
    PI_MS_septic_covid = ifelse(rs17580_A == 1 & Case == 1, 1, 0),
    #PI_SZ_septic_covid = ifelse(rs17580_A == 1 & rs28929474_T == 1 & Case == 1, 1, 0), no case
    not_severe_covid   = ifelse(Case == 0, 1, 0)
  )


# Define a reference group (usually the healthy/non-severe group)
bqc19_clean_data <- merged_data_prot %>%
  mutate(Analysis_Group = case_when(
    # Septic COVID groups
    PI_MM_septic_covid == 1 ~ "PI*MM\nSeptic COVID",
    PI_MZ_septic_covid == 1 ~ "PI*MZ\nSeptic COVID",
    PI_MS_septic_covid == 1 ~ "PI*MS\nSeptic COVID",

    # Not Septic COVID groups
    not_severe_covid == 1 & rs28929474_T == 0 & rs17580_A == 0 ~ "PI*MM\nNot Septic COVID",
    not_severe_covid == 1 & rs28929474_T == 1 & rs17580_A == 0 ~ "PI*MZ\nNot Septic COVID",
    not_severe_covid == 1 & rs28929474_T == 0 & rs17580_A == 1 ~ "PI*MS\nNot Septic COVID",

    # Catch-all for anyone who doesn't fit the above (e.g., PI*SZ, PI*ZZ, or missing data)
    TRUE ~ "Other"
  ))


# 1. Prepare the data
plot_data <- bqc19_clean_data %>%
  separate(col = Analysis_Group, 
           into = c("Genotype", "Status"), 
           sep = "\n", 
           remove = FALSE) %>%
  
  mutate(
    Status = factor(Status, levels = c("Not Septic COVID", "Septic COVID")),
    Genotype = factor(Genotype, levels = c("PI*MM", "PI*MS", "PI*MZ")),
    Analysis_Group = factor(Analysis_Group, levels = c(
      "PI*MM\nNot Septic COVID", "PI*MS\nNot Septic COVID", "PI*MZ\nNot Septic COVID",
      "PI*MM\nSeptic COVID", "PI*MS\nSeptic COVID", "PI*MZ\nSeptic COVID"
    ))
  )

AAT_plot<-ggplot(plot_data, aes(x = Genotype, y = `Mean SERPINA1`, fill = Analysis_Group)) +
  geom_quasirandom(shape = 21, size = 1.5, width = 0.2, alpha = 0.3) +
  geom_boxplot(width = 0.6, outlier.shape = NA) +
  
  # LAYER 3: The Magic Facet
  # This splits the plot by Status, and 'switch = "x"' moves the label to the bottom
  facet_grid(. ~ Status, switch = "x") +
  
  # LAYER 4: Styling
  theme_classic(base_size = 14) +

  scale_fill_manual(values = c('#a6cee3', '#b2df8a', '#fb9a99',   # Not Septic Colors
                               '#1f78b4', '#33a02c', '#e31a1c')) + # Septic Colors
  
  labs(
    x = NULL,
    y = "Mean AAT level (BQC19)"
  ) +
  
  theme(
    legend.position = "none", # Hides the legend entirely
    
    # These three lines format the bottom "Status" labels to look like a master x-axis
    strip.placement = "outside", 
    strip.background = element_blank(), 
    strip.text = element_text(size = 15,face = "bold"), 
    axis.title.y = element_text(size = 15, color = "black", face = "bold"),
    axis.text.x = element_text(size = 13, color = "black")
  )



#Wilcoxon Rank Sum tests independent and not paired
wilcox_result <- pairwise.wilcox.test(bqc19_clean_data$`Mean SERPINA1`,
                                      bqc19_clean_data$Analysis_Group,
                                      p.adjust.method = "bonferroni",
                                      paired = FALSE)

wilcox_result_adjusted<-wilcox_result$p.value

############################################################
###################UK GAinS cohort############################
############################################################
oxford_clean <- oxford_serpina1 %>%
  rename(Status = status,
         `Min SERPINA1` = min_value,
         `Mean SERPINA1`= mean_value,
         `Max SERPINA1` = max_value) %>%
  pivot_longer(
    cols = c(`Min SERPINA1`, `Mean SERPINA1`, `Max SERPINA1`),
    names_to = "norm_type",
    values_to = "SERPINA1"
  ) %>%
  # CRITICAL STEP: Filter to only one metric (Mean) to avoid triplicate points
  filter(norm_type == "Mean SERPINA1") %>% 
  mutate(
    Status = case_when(
      Status == "sepsis_non_carrier" ~ "PI*MM Sepsis",
      Status == "sepsis_carrier"     ~ "PI*MZ Sepsis",
      Status == "not_sepsis"         ~ "Not Sepsis",
      TRUE ~ Status
    )
  ) %>%
  # Set the order: Control -> MZ -> MM
  mutate(Status = factor(Status, levels = c("Not Sepsis", "PI*MZ Sepsis", "PI*MM Sepsis")))

# 2. GENERATE PLOT
oxford_plot <- ggplot(oxford_clean, aes(x = Status, y = SERPINA1)) +
  geom_quasirandom(aes(fill = Status),shape = 21, size = 1.5, width = 0.2, alpha = 0.3) +
  geom_boxplot(aes(fill = Status),width = 0.6, outlier.shape = NA) +
  theme_classic(base_size = 12) +
  scale_fill_manual(values = c('#66c2a5', '#fc8d62', '#8da0cb')) + # Using first 3 colors
  
  labs(x = "Status", 
       y = "Mean AAT level (UK GAinS Cohort)") +
  
  theme(legend.position = "none",
        axis.text.x = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 13, face = "bold"))
  

oxford_plot2<- oxford_plot + 
  geom_signif(
    comparisons = list(
      c("PI*MM Sepsis", "Not Sepsis"),
      c("PI*MM Sepsis", "PI*MZ Sepsis"),
      c("PI*MZ Sepsis", "Not Sepsis")
    ),
    map_signif_level = TRUE,
    annotation = c("", "", ""), # Force empty strings for all 3 brackets
    step_increase = 0.1,
    tip_length = 0.02
  )


kruskal_res<- kruskal.test(SERPINA1 ~ Status, data = oxford_clean)

#Wilcoxon Rank Sum tests independent and not paired
wilcox_result <- pairwise.wilcox.test(oxford_clean$SERPINA1, 
                                      oxford_clean$Status,
                                      #p.adjust.method = "bonferroni",
                                      paired = FALSE)

wilcox_result$p.value

