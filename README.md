# Human Genetic Analysis Reveals Circulating Alpha-1 Antitrypsin Level as a Protective Factor in Sepsis

Sepsis represents a dysregulated host response to infection and remains a leading cause of global mortality. Despite its prevalence, effective targeted therapies are lacking. This repository contains the analytical pipeline and code for a proteogenomic framework that integrates large-scale human genetics with circulating proteomics to identify and validate novel therapeutic targets for sepsis.

### Project Overview
This study shifts the focus toward genetic causal inference by employing a meta-analysis of genome-wide association studies (GWAS) encompassing over 1.5 million individuals. By leveraging Mendelian randomization (MR) and colocalization analyses, we provide human genetic evidence for the causal role of specific proteins in sepsis pathogenesis, prioritizing drug repurposing candidates.

#### Note on Code Reproducibility: The scripts in this repository are provided to assist researchers in replicating our results or to serve as a framework for similar analyses; however, they are not designed to function "out-of-the-box" and will require modification to align with your specific local computing environment, file paths, and data formats.

### Methodology

This repository includes the code for the following core analytical steps:

#### 1. GWAS & Primary Associations
- `01.beta_plots.R`: Generates visualizations of genetic effect sizes (betas) for the identified loci.
- `02.sepsis_mortality_logisticReg.R`: Conducts logistic regression models evaluating the association between key variants (e.g., SERPINA1 missense variant) and 30-day sepsis mortality.
- `03.LZ_plot_sepsis.R`: Creates LocusZoom plots to visualize the regional association landscape around genome-wide significant loci.

#### 2. Mendelian Randomization (MR) Analyses
- `04.1.AAT_on_other_proteins_MR.R`: Two-sample MR estimating the causal effect of circulating AAT levels on the broader proteome.
- `04.2.other_proteins_on_AAT_MR.R`: Reverse MR analyses evaluating the effects of other circulating proteins on AAT levels.
- `04.3.plot_protein_protein_MR_results.R`: Generates visualizations (e.g., heatmap) summarizing the protein-to-protein MR findings.
- `05.MR_AAT_on_other_traits.R`: Evaluates the causal effect of AAT on various acute infectious and non-infectious clinical phenotypes to assess pleiotropy and specificity.

#### 3. Covariate & Cohort Validations
- `06.regression_smoking_diseases.R`: Regression models adjusting for key clinical covariates and comorbidities (e.g., smoking status, underlying diseases).
- `07.1.septic_cohort_validations.R`: Scripts for analyzing AAT levels during acute illness in independent validation cohorts (UK GAinS and BQC19).
- `07.2.protease_antiprotease_comp.R`: Analyzes the dose-dependent attenuation of AAT in missense variant carriers to evaluate the protease-antiprotease balance.

### Datasets
The analyses within this repository utilize data from several large-scale biobanks and clinical cohorts:
- UK Biobank
- UK Genomic Advances in Sepsis (GAinS)
- Biobanque Québécois sur la COVID-19 (BQC19)

*(Note: Individual-level data is subject to data use agreements and is not hosted in this repository. Summary statistics will be made available upon publication).*

### Citation
If you use this code or build upon our findings in your own research, please cite our preprint:

> **Human Genetic Analysis Reveals Circulating Alpha-1 Antitrypsin Level as a Protective Factor in Sepsis** *medRxiv* (2026). Available at: https://doi.org/10.64898/2026.03.25.26349312

### License
Distributed under the MIT License. See `LICENSE` for more information.
