library(dplyr)
library(data.table)
library(vroom)
library(purrr)
library(ggplot2)
library(broom)
library(gridExtra)
library(cowplot)

geno_chr14<-read.table("genotypefile", header = TRUE)

covariates<-covariates %>% select(age_enrolment,sex,age_sex,age2,smoking_status,eid) 
mortaility<-mortaility %>% rename(eid=Participant.ID)

sepsis_data_all <- disease_information %>%
  mutate(
    sepsis_diagnose = ifelse(rowSums(!is.na(cbind(p130068, p130070, p132286))) > 0, 1, 0),
    copd_diagnose = ifelse(!is.na(p131486) | !is.na(p131488) | !is.na(p131490) | !is.na(p131492), 1, 0),
    # takes the earliest date from the 4 columns
    copd_diagnose_date = as.Date(
      pmin(
        as.numeric(p131486), 
        as.numeric(p131488), 
        as.numeric(p131490), 
        as.numeric(p131492), 
        na.rm = TRUE
      ),
      origin = "1970-01-01"
    ),
    
    cirrhosis_diagnose = as.integer(!is.na(p131666)),
    cirrhosis_diagnose_date = p131666,
    
    bronchiectasis_diagnose = as.integer(!is.na(p131498)),
    bronchiectasis_diagnose_date = p131498
  ) %>%
  select(
    copd_diagnose, copd_diagnose_date, 
    cirrhosis_diagnose, cirrhosis_diagnose_date, 
    bronchiectasis_diagnose, bronchiectasis_diagnose_date,
    p130068, p130070, p132286, sepsis_diagnose, eid
  ) %>%
  right_join(geno_chr14, by = "eid") %>%
  mutate(
    # Ensure binary indicators are 0 instead of NA after joins
    copd_diagnose = ifelse(is.na(copd_diagnose), 0, copd_diagnose),
    cirrhosis_diagnose = ifelse(is.na(cirrhosis_diagnose), 0, cirrhosis_diagnose),
    bronchiectasis_diagnose = ifelse(is.na(bronchiectasis_diagnose), 0, bronchiectasis_diagnose),
    sepsis_diagnose = ifelse(is.na(sepsis_diagnose), 0, sepsis_diagnose)
  ) %>%
  left_join(distinct(mortaility, eid, .keep_all = TRUE), by = "eid") %>%
  left_join(covariates, by = "eid")

sepsis_data_all <- sepsis_data_all %>%
  mutate(
    sepsis_diagnose_date = as.Date(
      pmin(
        as.numeric(p130068), 
        as.numeric(p130070), 
        as.numeric(p132286), 
        na.rm = TRUE
      ),
      origin = "1970-01-01"
    ),
    
    # 2. Complete the ifelse statements
    copd_diagnose_before_sepsis = ifelse(
      !is.na(copd_diagnose_date) & !is.na(sepsis_diagnose_date) & copd_diagnose_date < sepsis_diagnose_date, 
      1, 0
    ),
    
    cirrhosis_diagnose_before_sepsis = ifelse(
      !is.na(cirrhosis_diagnose_date) & !is.na(sepsis_diagnose_date) & cirrhosis_diagnose_date < sepsis_diagnose_date, 
      1, 0
    ),
    
    bronchiectasis_diagnose_before_sepsis = ifelse(
      !is.na(bronchiectasis_diagnose_date) & !is.na(sepsis_diagnose_date) & bronchiectasis_diagnose_date < sepsis_diagnose_date, 
      1, 0
    )
  )

sepsis_data <- sepsis_data_all %>%
  mutate(across(c(p130068, p130070, p132286, Date.of.death), as.Date))

sepsis_data <- sepsis_data %>%
  mutate(diagnose_date = pmap_chr(
    list(p130068, p130070, p132286),
    ~ suppressWarnings(
      as.character(max(as.Date(c(...)), na.rm = TRUE))
      )
    ), 
    diagnose_date = as.Date(diagnose_date))

# Compute binary death outcome
sepsis_data <- sepsis_data %>%
  mutate(
    death_within_35 = if_else(
      !is.na(Date.of.death) & !is.na(diagnose_date) &
        as.numeric(Date.of.death - diagnose_date) >= -5 &
        as.numeric(Date.of.death - diagnose_date) <= 30,
      1L,
      0L
    )
  )

sepsis_cases<-sepsis_data %>%
  filter(sepsis_diagnose == 1)

sepsis_cases$sex <- factor(sepsis_cases$sex)
sepsis_cases$smoking_status <- factor(sepsis_cases$smoking_status)

sepsis_data$sex <-factor(sepsis_data$sex)
sepsis_data$smoking_status <- factor(sepsis_data$smoking_status)

sum(sepsis_data$sepsis_diagnose == 1, na.rm = TRUE)
sum(sepsis_data$death_within_35 ==1, na.rm = TRUE)


model_smoking_risk_interaction <- glm(
  sepsis_diagnose ~ chr14.94378610.C.T_T * smoking_status+ age_enrolment + sex + age2,
  data = sepsis_data,
  family = "binomial"
)
summary_df_smoking <- as.data.frame(summary(model_smoking_risk_interaction)$coefficients)
colnames(summary_df_smoking) <- c("Estimate", "Std_Error", "z_value", "P_value")

model_other_disease_risk <- glm(
  sepsis_diagnose ~ chr14.94378610.C.T_T+ copd_diagnose_before_sepsis + bronchiectasis_diagnose_before_sepsis + cirrhosis_diagnose_before_sepsis + age_enrolment + sex + age2,
  data = sepsis_data,
  family = "binomial" 
)
summary_df_other_diseases <- as.data.frame(summary(model_other_disease_risk)$coefficients)
colnames(summary_df_other_diseases) <- c("Estimate", "Std_Error", "z_value", "P_value")

