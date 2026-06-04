library(ggplot2)
library(dplyr)
library(tidyverse)
library(purrr)
library(forcats)
library(ggnewscale)
library(extrafont)
library(readxl)
library(ggpubr)

#####################select genes to plot #################################
results<-read.table("../SERPINA1_on_others_decode.tsv",sep = "\t",header = TRUE)

MR_decode<- results %>% mutate(p_adj_bon = p.adjust(p_val, method = "bonferroni"))
MR_decode <- MR_decode %>% rename(
  outcome = outcome_file,
  beta_decode = beta,
  se_decode = se,
  odds_ratio_decode = odd_ratio,
  odd_lower_decode = odd_lower,
  odd_upper_decode = odd_upper,
  p_val_decode = p_val,
  p_adj_bon_decode = p_adj_bon
)
significant_MR_decode <- MR_decode %>% filter(p_adj_bon_decode <0.05)

#filter with UKB PPP sig proteins 
significant_MR_UKBPPP<-read.csv("../SERPINA1_interaction_other_proteins.csv",header = TRUE)

MR_protein_comb <- significant_MR_UKBPPP %>% left_join(significant_MR_decode,by="outcome")
MR_protein_comb <- MR_protein_comb %>%
  mutate(same_dir = ifelse(sign(beta) == sign(beta_decode), 1, 0))

protein_function <- read_excel("../sepsis_protein_func.xlsx")

protein_function<- protein_function %>% select(Proteins,`Subcategory (for sepsis)`) %>%
  rename(Subcategory = `Subcategory (for sepsis)`,
         outcome = Proteins)

MR_protein_comb <- MR_protein_comb %>% left_join(protein_function, by ="outcome")

MR_protein_comb$Subcategory[MR_protein_comb$exposure == "SERPINA6"] <- "Metabolism / Stress"

#############################
MR_protein_comb_plot<- MR_protein_comb %>% 
  mutate(method = "Mendelian Randomization",outcome.name = "AAT",mr.p_cat =1) %>%
  rename(exposure.name = outcome,
         domain.exposure = Subcategory,
         decode_sig_validation = same_dir) %>% 
  mutate(mr.decode_val = ifelse(decode_sig_validation==1, "*"," ")) %>%
  select(method, outcome.name, exposure.name, domain.exposure, beta, se, p_val, p_adj_bon,mr.decode_val,mr.p_cat)
MR_protein_comb_plot$exposure.name[MR_protein_comb_plot$exposure.name == "SERPINA1"] <- "SERPINA6"
MR_protein_comb_plot$mr.decode_val[MR_protein_comb_plot$exposure.name == "SERPINA6"] <- "#"
MR_protein_comb_plot <- MR_protein_comb_plot %>%
  mutate(domain.exposure = factor(domain.exposure, levels = c("Coagulation / Hemostasis", "Endothelial / Vascular Integrity", "Metabolism / Stress", 
                                                              "Immune Response / Inflammation", "Tissue Damage / Repair")))

MR_protein_comb_plot <- MR_protein_comb_plot %>%
  mutate(beta = pmin(pmax(beta, -0.5), 0.5))
MR_protein_comb_plot_1<- MR_protein_comb_plot %>% filter (domain.exposure =="Coagulation / Hemostasis" | domain.exposure =="Endothelial / Vascular Integrity"| domain.exposure=="Metabolism / Stress")
MR_protein_comb_plot_2<- MR_protein_comb_plot %>% filter (domain.exposure =="Immune Response / Inflammation" | domain.exposure =="Tissue Damage / Repair")


loadfonts()
text_size = 7
theme.size = 1.25
astrix.size = 2
line.size = 0.25

data.long<-MR_protein_comb_plot_1
plot1 <- ggplot(data.long) +
  facet_grid(method ~ domain.exposure, scales = "free", space = "free", switch = "y") +
  
  # MR tiles
  geom_tile(data = filter(data.long, method == "Mendelian Randomization" & !is.na(beta)),
            aes(x = exposure.name, y = outcome.name, fill = beta,
                height = mr.p_cat, width = mr.p_cat)) +
  geom_text(data = filter(data.long, method == "Mendelian Randomization" & !is.na(beta)),
            aes(label = mr.decode_val, x = exposure.name, y = outcome.name),
            vjust = 0.75, size = astrix.size) +
  
  # Color scale
  scale_fill_distiller(
    palette = "RdBu",
    direction = -1,
    limits = c(-0.5, 0.5),
    na.value = "white"
  ) +
  
  # Grid lines
  geom_vline(xintercept = seq(0.5, length(unique(data.long$exposure.name)) + 0.5, 1),
             color = "grey90", size = line.size) +
  geom_hline(yintercept = seq(0.5, length(unique(data.long$outcome.name)) + 0.5, 1),
             color = "grey90", size = line.size) +
  
  # Theme
  theme_classic() +
  theme(
    legend.position = 'none',
    axis.text.x = element_text(angle = 45, hjust = 0),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    plot.margin = margin(10, 10, 10, 10, 'pt'),
    text = element_text(size = text_size, family = "sans"),  # safer than Arial
    axis.ticks = element_line(size = line.size),
    axis.line = element_line(size = line.size),
    strip.background = element_blank(),
    
    # Hide method label (y-axis strip)
    strip.text.y.left = element_blank(),
    
    # Bold domain.exposure label (top strip)
    strip.text.x = element_text(face = "bold", size = 6)
  ) +
  
  scale_x_discrete(position = "top")

data.long<-MR_protein_comb_plot_2
plot2 <- ggplot(data.long) +
  facet_grid(method ~ domain.exposure, scales = "free", space = "free", switch = "y") +
  
  # MR tiles
  geom_tile(data = filter(data.long, method == "Mendelian Randomization" & !is.na(beta)),
            aes(x = exposure.name, y = outcome.name, fill = beta,
                height = mr.p_cat, width = mr.p_cat)) +
  geom_text(data = filter(data.long, method == "Mendelian Randomization" & !is.na(beta)),
            aes(label = mr.decode_val, x = exposure.name, y = outcome.name),
            vjust = 0.75, size = astrix.size) +
  
  # Color scale
  scale_fill_distiller(
    palette = "RdBu",
    direction = -1,
    limits = c(-0.5, 0.5),
    na.value = "white"
  ) +
  
  # Grid lines
  geom_vline(xintercept = seq(0.5, length(unique(data.long$exposure.name)) + 0.5, 1),
             color = "grey90", size = line.size) +
  geom_hline(yintercept = seq(0.5, length(unique(data.long$outcome.name)) + 0.5, 1),
             color = "grey90", size = line.size) +
  
  # Theme
  theme_classic() +
  theme(
    legend.position = 'none',
    axis.text.x = element_text(angle = 45, hjust = 0),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    plot.margin = margin(10, 10, 10, 10, 'pt'),
    text = element_text(size = text_size, family = "sans"),  # safer than Arial
    axis.ticks = element_line(size = line.size),
    axis.line = element_line(size = line.size),
    strip.background = element_blank(),
    
    # Hide method label (y-axis strip)
    strip.text.y.left = element_blank(),
    
    # Bold domain.exposure label (top strip)
    strip.text.x = element_text(face = "bold", size = 6)
  ) +
  
  scale_x_discrete(position = "top")


# Stack vertically Mendelian Randomization (MR) results for Bonferroni-adjusted p < 0.05, using SERPINA1 as the exposure and other proteins as outcomes. 
# Red indicates positive effect sizes; blue indicates negative effect sizes. * denotes validation in the DeCODE cohort (Bonferroni-corrected). '#' indicates SERPINA1 used as the outcome.
combined_plot <- ggarrange(
  plot1,
  plot2,
  ncol = 1,        
  nrow = 2,        
  align = "v",    
  heights = c(1, 1)
)


