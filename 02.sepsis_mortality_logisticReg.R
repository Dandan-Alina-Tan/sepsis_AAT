library(dplyr)
library(data.table)
library(vroom)
library(purrr)
library(ggplot2)
library(broom)
library(gridExtra)
library(cowplot)

#make sure the higher dose(2) increase the risk of sepsis: positive beta
geno_chr14<-read.table("/path/to/genotypefile", header = TRUE)
geno_chr19<-read.table("/path/to/genotypefile", header = TRUE)
geno_chr2<-read.table("/path/to/genotypefile", header = TRUE)
geno_chr6<-read.table("/path/to/genotypefile", header = TRUE)

# sepsis_data have sepsis_diagnose (1/0) and the date of diagnose (2000-01-01), date of mortality (2000-01-01), as well as death_within_35 (1/0)
sepsis_data <- vroom("/path/to/file")
sepsis_cases <- sepsis_data %>%
  filter(sepsis_diagnose == 1)

sepsis_cases$sex <- factor(sepsis_cases$sex)

model_snps_combined_case_control <- glm(
  death_within_35 ~ chr14.94378610.C.T_T + chr19.33096952.C.T_T + chr2.632609.A.G_A + chr6.32224045.G.A_A + age_enrolment + sex + age2,
  data = sepsis_data,
  family = "binomial"
)

model_snps_combined_case <- glm(
  death_within_35 ~ chr14.94378610.C.T_T + chr19.33096952.C.T_T + chr2.632609.A.G_A + chr6.32224045.G.A_A + age_enrolment + sex + age2,
  data = sepsis_cases,
  family = "binomial"
)

res_case <- tidy(model_snps_combined_case, conf.int = TRUE)
res_all  <- tidy(model_snps_combined_case_control, conf.int = TRUE)

# add model labels
res_case$model <- "Sepsis mortality in individuals with sepsis"
res_all$model  <- "Sepsis mortality in all participants"
# combine
res <- rbind(res_case, res_all)
res_snps <- res[grepl("^chr", res$term), ]

rename_vars <- c(
  "chr14.94378610.C.T_T" = "chr14:94378610:C:T",
  "chr19.33096952.C.T_T" = "chr19:33096952:C:T",
  "chr2.632609.A.G_A"    = "chr2:632609:A:G",
  "chr6.32224045.G.A_A"  = "chr6:32224045:G:A"
)

res_snps$term <- dplyr::recode(res_snps$term, !!!rename_vars)
res_snps <- res_snps %>% 
  dplyr::mutate(
    OR = exp(estimate),
    OR_low = exp(conf.low),
    OR_high = exp(conf.high)
  )

order <- c("chr14:94378610:C:T", "chr2:632609:A:G", "chr19:33096952:C:T", "chr6:32224045:G:A")
res_snps$term <- factor(res_snps$term, levels = rev(order))  # for coord_flip

res_snps$model <- factor(
  res_snps$model,
  levels = c("Sepsis mortality in individuals with sepsis", 
             "Sepsis mortality in all participants")  # legend order
)
pd <- position_dodge2(width = 0.5, preserve = "single")

mortality_plot <- ggplot(
  res_snps,
  aes(
    x = OR,
    y = term,
    color = model
  )
) +
  geom_point(position = pd, size = 3) +
  geom_errorbarh(
    aes(xmin = OR_low, xmax = OR_high),
    position = pd,
    height = 0.5,
    linewidth = 0.7
  ) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey40") +
  scale_color_manual(values = c("#ef8a62", "#67a9cf")) +
  guides(color = guide_legend(reverse = TRUE)) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "right",
    legend.title = element_blank(),
    legend.text = element_text(size = 8),
    legend.key.size = unit(0.6, "cm"),  # smaller keys
    axis.text.y = element_text(size = 14, color = "black"),
    axis.text.x = element_text(size = 14, color = "black"),
    plot.title = element_text(face = "bold", size = 15),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.ticks = element_line(color = "black", linewidth = 0.5)
  ) +
  labs(
    x = "Odds Ratio (95% CI)",
    y = "SNP"
  )

