# ----------------------------------------------------------------------------------------------------------------------------------
# Data analysis and results
# ----------------------------------------------------------------------------------------------------------------------------------
# 
# The following code summarizes the results from both assays.
# 
# - Traits assay: for each prey trait and for the prey growth rate, a linear
#   mixed model (LMM) is fitted, with selection and assay stressors as fixed
#   effects, and with selection line as random effect.
#   Post-hoc tests for investigating hypothesis 3 ("evolution in the multi-
#   stressor environment is hindered compared to that in single-stressor 
#   environments") are fitted to growth rates and trait data within single-
#   stressor assay environments.
# 
# - Functional response assay: the summary table of the model fitted in "R scripts/02_FR_search_grid.R" 
#   is given again, so that all results can be checked within the same R script.
# 
# NOTE: this script was run on R version 4.5.2 (provided in folder
# "R versions/R version for data analysis"). It will not run on R version 4.3.2 due to package
# compatibility issues.
# 
# ----------------------------------------------------------------------------------------------------------------------------------

rm(list = ls())

# Libraries
library(tidyverse)
library(lmerTest)
library(broom)
library(emmeans)

# Data
load("Data/Traits_assay.RData")

# --------------------------------------------------------------------------------------------------
# Hypothesis 1: prey growth rates and traits exhibit plastic responses to assay stressors
# --------------------------------------------------------------------------------------------------

# Growth rates
GR_4_way_LMM <- lmer(GR ~ Selection_T * Selection_S * Assay_T * Assay_S + (1|Selection_line),
                     traits_assay)
summary(GR_4_way_LMM)

# Cell size
size_4_way_LMM <- lmer(Mean_cell_area ~ Selection_T * Selection_S * Assay_T * Assay_S + (1|Selection_line),
                       traits_assay)
summary(size_4_way_LMM)

# Movement speed
speed_4_way_LMM <- lmer(Mean_movement_speed ~ Selection_T * Selection_S * Assay_T * Assay_S + (1|Selection_line),
                        traits_assay)
summary(speed_4_way_LMM)

# --------------------------------------------------------------------------------------------------
# Hypothesis 2: evolutionary responses of prey traits and growth rates to single stressors weaken
# their plastic responses to the same stressor
# 
# Hypothesis 3: a second stressor in the selection environment will decrease the evolutionary
# response to the first stressor
# --------------------------------------------------------------------------------------------------

# NOTE: to facilitate the interpretation of the emmeans results, we re-fit the linear mixed models
# with selection environment and assay environment encoded as factors with levels C, S, T, and ST,
# instead of using selection temperature, selection salinity, assay temperature, and assay salinity.
# This does not affect the results in any way, other than changing the names of the variables.

# Relevant contrasts
cont <- list(
  deltaS = c(-1, 1, 0, 0),    # deltaS = S - C
  deltaT = c(-1, 0, 1, 0),    # deltaT = T - C
  deltaST = c(-1, 0, 0, 1),   # deltaST = ST - C
  `deltaST - deltaS` = c(0, -1, 0, 1),   # deltaST - deltaS = (ST - C) - (S - C) = ST - S
  `deltaST - deltaT` = c(0, 0, -1, 1)    # deltaST - deltaT = (ST - C) - (T - T) = ST - T
)

# --------------- Growth rate
GR_2_way_LMM <- lmer(GR ~ Selection_environment * Assay_environment + (1|Selection_line),
                     traits_assay)
emm_GR <- emmeans(GR_2_way_LMM, ~ Selection_environment | Assay_environment)

# Checking evolutionary responses within single-stressor assay environments (Table 1 in main text)
contrast(emm_GR, cont, by = "Assay_environment") %>%
  as.data.frame() %>%
  filter((Assay_environment == "S" & contrast %in% c("deltaS", "deltaST", "deltaST - deltaS")) |
           (Assay_environment == "T" & contrast %in% c("deltaT", "deltaST", "deltaST - deltaT")))

# --------------- Cell size
size_2_way_LMM <- lmer(Mean_cell_area ~ Selection_environment * Assay_environment + (1|Selection_line),
                       traits_assay)
emm_size <- emmeans(size_2_way_LMM, ~ Selection_environment | Assay_environment)

# Checking evolutionary responses within single-stressor assay environments (Table 1 in main text)
contrast(emm_size, cont, by = "Assay_environment") %>%
  as.data.frame() %>%
  filter((Assay_environment == "S" & contrast %in% c("deltaS", "deltaST", "deltaST - deltaS")) |
           (Assay_environment == "T" & contrast %in% c("deltaT", "deltaST", "deltaST - deltaT")))

# --------------- Movement speed
speed_2_way_LMM <- lmer(Mean_movement_speed ~ Selection_environment * Assay_environment + (1|Selection_line),
                       traits_assay)
emm_speed <- emmeans(speed_2_way_LMM, ~ Selection_environment | Assay_environment)

# Checking evolutionary responses within single-stressor assay environments (Table 1 in main text)
contrast(emm_speed, cont, by = "Assay_environment") %>%
  as.data.frame() %>%
  filter((Assay_environment == "S" & contrast %in% c("deltaS", "deltaST", "deltaST - deltaS")) |
           (Assay_environment == "T" & contrast %in% c("deltaT", "deltaST", "deltaST - deltaT")))

# save(emm_GR, emm_size, emm_speed,
#      file = "Data/Marginal_means_estimates.RData")

# --------------------------------------------------------------------------------------------------
# Functional response assay
# --------------------------------------------------------------------------------------------------

# Model fitting to FR assay data was done in "R scripts/02_FR_search_grid.R"
load("Data/FR assay/03 - FR search grid/FR_ANOVA_like_best_fit.RData")

(FR_coefs <- generics::tidy(best_FR_fit))

# Getting the coefficients for each prey selection environment (used in 05_Manuscript_figures.R)
vcov_matrix <- bbmle::vcov(best_FR_fit)

FR_coefs <- generics::tidy(best_FR_fit)

FRC_coefs <- FR_coefs[grep("C", FR_coefs$term),]
FRS_coefs <- FR_coefs[grep("C", FR_coefs$term) + 1,]
FRT_coefs <- FR_coefs[grep("C", FR_coefs$term) + 2,]
FRST_coefs <- FR_coefs[grep("ST", FR_coefs$term),]

# Computing the sums for each parameter in each environment
FRST_coefs$estimate <- FRST_coefs$estimate + FRS_coefs$estimate + FRT_coefs$estimate + FRC_coefs$estimate
FRST_coefs$std.error <- c(sqrt(sum(vcov_matrix[c("bC.log", "bS.log", "bT.log", "bST.log"), c("bC.log", "bS.log", "bT.log", "bST.log")])),
                          sqrt(sum(vcov_matrix[c("hC", "hS", "hT", "hST"), c("hC", "hS", "hT", "hST")])),
                          sqrt(sum(vcov_matrix[c("mC", "mS", "mT", "mST"), c("mC", "mS", "mT", "mST")])))


FRS_coefs$estimate <- FRS_coefs$estimate + FRC_coefs$estimate
FRS_coefs$std.error <- c(sqrt(sum(vcov_matrix[c("bC.log", "bS.log"), c("bC.log", "bS.log")])),
                         sqrt(sum(vcov_matrix[c("hC", "hS"), c("hC", "hS")])),
                         sqrt(sum(vcov_matrix[c("mC", "mS"), c("mC", "mS")])))


FRT_coefs$estimate <- FRT_coefs$estimate + FRC_coefs$estimate
FRT_coefs$std.error <- c(sqrt(sum(vcov_matrix[c("bC.log", "bT.log"), c("bC.log", "bT.log")])),
                         sqrt(sum(vcov_matrix[c("hC", "hT"), c("hC", "hT")])),
                         sqrt(sum(vcov_matrix[c("mC", "mT"), c("mC", "mT")])))

# Saving coefficient estimates for plotting (in 05_Manuscript_figures.R)
# save(FRC_coefs, FRS_coefs, FRT_coefs, FRST_coefs,
#      file = "Data/FR assay/04 - FR fits and simulations/Best_fit_estimates_per_environment.RData")
