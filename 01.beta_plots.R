library(ggplot2)
library(dplyr)
library(gridExtra)
library(grid)
library(gtable)

ukb_eur_gwas<-read.table("extracted_snps_eur/ukb_eur_matched_snps.tsv",header = TRUE)
mvp_eur_gwas<-read.table("extracted_snps_eur/mvp_eur_matched_snps.tsv",header = TRUE)
finngen_gwas<-read.table("extracted_snps_eur/finngen_matched_snps.tsv",header = TRUE)
aou_eur_gwas<-read.table("extracted_snps_eur/aou_eur_matched_snps.tsv",header = TRUE)

eur_meta<-read.table("extracted_snps_eur/eur_meta_matched_snps.tsv",header = TRUE)
eur_meta<-eur_meta %>% rename(Name=MarkerName,
                              alt=effect_allele,
                              ref=non_effect_allele,
                              eaf=effect_allele_frequency)

ukb_eur_gwas$Cohort <- "UKB_EUR"
mvp_eur_gwas$Cohort <- "MVP_EUR"
finngen_gwas$Cohort <- "FinnGen_EUR"
aou_eur_gwas$Cohort <- "AoU_EUR"
eur_meta$Cohort <- "Meta_EUR"


snplist <- unique(eur_meta$Name)
snp_df <- data.frame(SNP = snplist)

standardize_gwas <- function(df, cohort_name) {
  df <- df %>%
    rename_with(tolower) %>%
    rename(SNP = name, BETA = beta, SE = se) %>%
    select(SNP, BETA, SE, ref, alt, eaf) %>%
    mutate(Cohort = cohort_name)
  
  full_join(snp_df, df, by = "SNP")
}

# Apply to each GWAS dataset
gwas_list <- list(
  standardize_gwas(ukb_eur_gwas, "UKB_EUR"),
  standardize_gwas(mvp_eur_gwas, "MVP_EUR"),
  standardize_gwas(finngen_gwas, "FinnGen_EUR"),
  standardize_gwas(aou_eur_gwas, "AoU_EUR"),
  standardize_gwas(eur_meta, "Meta_EUR")
)
# Combine all cohorts into one dataframe
all_data <- bind_rows(gwas_list)

# Calculate 95% confidence intervals
all_data <- all_data %>%
  mutate(
    CI_lower = BETA - 1.96 * SE,
    CI_upper = BETA + 1.96 * SE
  )


all_data$SNP <- factor(all_data$SNP, levels = snplist)
cohort_order <- c("UKB_EUR", "MVP_EUR", "FinnGen_EUR","AoU_EUR",
                  "Meta_EUR")

custom_colors<-c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99')
####### Plot: Faceted by SNP
all_data_filtered <- all_data %>%
  filter(!is.na(BETA) & !is.na(SE)) %>%
  mutate(Cohort = factor(Cohort, levels = rev(cohort_order)))

custom_colors<-c('#e41a1c','#377eb8','#4daf4a','#984ea3','#a6cee3')
order<-c("chr14:94378610:C:T","chr2:632609:A:G","chr19:33096952:C:T","chr6:32224045:G:A")

# Reorder SNP as a factor
all_data_filtered$SNP <- factor(all_data_filtered$SNP,
                                levels = c("chr6:32224045:G:A",
                                           "chr14:94378610:C:T",
                                           "chr19:33096952:C:T",
                                           "chr2:632609:A:G"))

all_data_filtered <- all_data_filtered %>%
  mutate(
    OR = exp(BETA),
    OR_upper = exp(BETA + 1.96 * SE),
    OR_lower = exp(BETA - 1.96 * SE)
  )

OR_forest_plot <- ggplot(all_data_filtered, aes(x = OR, y = Cohort, color = Cohort)) +
  geom_point(position = position_dodge(width = 0.6), size = 2, na.rm = TRUE) +
  geom_errorbarh(aes(xmin = OR_lower, xmax = OR_upper),
                 height = 0.3, position = position_dodge(width = 0.6), na.rm = TRUE) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray50") +
  
  # ADDED: This draws the line above Meta_EUR (which is at y = 1)
  geom_hline(yintercept = 1.5, linetype = "solid", color = "black", linewidth = 0.5) +
  
  facet_wrap(~ SNP, scales = "free_x", ncol = 4) + 
  scale_color_manual(values = custom_colors) +
  theme_bw(base_size = 10) +
  labs(x = "Odds Ratio (95% CI)", y = "Cohort") +
  theme(strip.text = element_text(size = 11),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 8),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        legend.position = "none",
        panel.spacing.x = unit(1.1, "lines"))

OR_forest_plot <- OR_forest_plot +
  theme(
    panel.spacing.x = unit(1.1, "lines")  # horizontal space
    #panel.spacing.y = unit(1, "lines")   # vertical space
  )
ggsave("../Results/OR_forest_plot_eur.pdf",
       OR_forest_plot,width=8,height=2,dpi=300)
