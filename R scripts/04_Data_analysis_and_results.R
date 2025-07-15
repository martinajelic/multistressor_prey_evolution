# ----------------------------------------------------------------------------------------------------------------------------------
# Data analysis and results
# ----------------------------------------------------------------------------------------------------------------------------------
# 
# The following code summarizes the results from both assays.
# 
# - FR assay: the summary table of the model fitted in "R scripts/02_FR_search_grid.R" 
#   is given again, so that all results can be checked within the same R script.
# 
# - Traits assay: for each prey trait and for the prey growth rate, a linear
#   mixed model (LMM) is fitted, with selection and assay stressors as fixed
#   effects, and with selection line as random effect.
#   Post-hoc tests for investigating hypothesis 3 ("evolution in the multi-
#   stressor environment is hindered compared to that in single-stressor 
#   environments") are fitted to growth rates and trait data within single-
#   stressor assay environments.
# 
# ----------------------------------------------------------------------------------------------------------------------------------

rm(list = ls())

library(tidyverse)
library(lmerTest)

# --------------------------------------------------------------------------------------------------
# FR assay
# --------------------------------------------------------------------------------------------------

# Model fitting to FR assay data was done in "R scripts/02_FR_search_grid.R"
load("Data/FR assay/03 - FR search grid/FR_ANOVA_like_best_fit.RData")

summary(best_FR_fit)

# --------------------------------------------------------------------------------------------------
# Traits assay
# --------------------------------------------------------------------------------------------------

load("Data/Traits_assay.RData")

# ------------------------------------------
# 4-way LMM of growth rate data
# ------------------------------------------
# Visualization of growth rate data
# ggplot() +
#   geom_boxplot(data = traits_assay, aes(x = Assay_environment, y = GR,
#                                                           color = Selection_environment, fill = Selection_environment),
#                alpha = 0.3, linewidth = 0.2) +
#   geom_point(data = traits_assay, aes(x = Assay_environment, y = GR,
#                                                         color = Selection_environment, fill = Selection_environment, shape = Biorep),
#              alpha = 1, size = 1, position = position_jitterdodge(jitter.width = 0.09, dodge.width = 0.75)) +
#   theme_bw() +
#   labs(x = "Assay environment",
#        y = latex2exp::TeX("Growth rate $ln(N_{t}/N_{0})/\\Delta t$"),
#        color = "Selection environment",
#        fill = "Selection environment",
#        shape = "Replicate") +
#   scale_color_manual(values = c("#000000", "#E69F00", "#56B4E9", "#009E73")) +
#   scale_fill_manual(values = c("#000000", "#E69F00", "#56B4E9", "#009E73")) +
#   geom_hline(yintercept = 0, linetype = 2, color = "gray") +
#   theme(strip.background = element_blank())
# ------------------------------------------

GR_4_way_LMM <- lmer(GR ~ Selection_T * Selection_S * Assay_T * Assay_S + (1|Selection_line),
                     traits_assay)
summary(GR_4_way_LMM)

# ------------------------------------------
# 4-way LMM of cell size data
# ------------------------------------------
# Visualization of cell size data
# ggplot() +
#   geom_boxplot(data = traits_assay, aes(x = Assay_environment, y = Mean_cell_area,
#                                                           color = Selection_environment, fill = Selection_environment),
#                alpha = 0.3, linewidth = 0.2) +
#   geom_point(data = traits_assay, aes(x = Assay_environment, y = Mean_cell_area,
#                                                         color = Selection_environment, fill = Selection_environment, shape = Biorep),
#              alpha = 1, size = 1, position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.75)) +
#   theme_bw() +
#   labs(x = "Assay environment",
#        y = latex2exp::TeX("Cell area $(µm^2)$"),
#        color = "Selection environment",
#        fill = "Selection environment",
#        shape = "Replicate") +
#   scale_color_manual(values = c("#000000", "#E69F00", "#56B4E9", "#009E73")) +
#   scale_fill_manual(values = c("#000000", "#E69F00", "#56B4E9", "#009E73")) +
#   theme(strip.background = element_blank())
# ------------------------------------------

size_4_way_LMM <- lmer(Mean_cell_area ~ Selection_T * Selection_S * Assay_T * Assay_S + (1|Selection_line),
                       traits_assay)
summary(size_4_way_LMM)

# ------------------------------------------
# 4-way LMM of movement speed data
# ------------------------------------------
# Visualization of movement speed data
# ggplot() +
#   geom_boxplot(data = traits_assay, aes(x = Assay_environment, y = Mean_movement_speed,
#                                                           color = Selection_environment, fill = Selection_environment),
#                alpha = 0.3, linewidth = 0.2) +
#   geom_point(data = traits_assay, aes(x = Assay_environment, y = Mean_movement_speed,
#                                                         color = Selection_environment, fill = Selection_environment, shape = Biorep),
#              alpha = 1, size = 1, position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.75)) +
#   theme_bw() +
#   labs(x = "Assay environment",
#        y = latex2exp::TeX("Movement speed ($µm$ $sec^{-1}$)"),
#        color = "Selection environment",
#        fill = "Selection environment",
#        shape = "Replicate") +
#   scale_color_manual(values = c("#000000", "#E69F00", "#56B4E9", "#009E73")) +
#   scale_fill_manual(values = c("#000000", "#E69F00", "#56B4E9", "#009E73")) +
#   theme(strip.background = element_blank())
# ------------------------------------------

speed_4_way_LMM <- lmer(Mean_movement_speed ~ Selection_T * Selection_S * Assay_T * Assay_S + (1|Selection_line),
                        traits_assay)
summary(speed_4_way_LMM)

# ------------------------------------------
# Post-hoc analyses within single-stressor assay environments
# ------------------------------------------
# Is evolution in multi-stressor environments hindered compared to evolution in single-stressor
# environments? (Hypothesis 3 in the main text)

# Model fitting
GR_MS_evo_S <- lmer(GR ~ Selection_T * Selection_S + (1|Selection_line),
                    traits_assay %>% filter(Assay_environment == "S"))
size_MS_evo_S <- lmer(Mean_cell_area ~ Selection_T * Selection_S + (1|Selection_line),
                    traits_assay %>% filter(Assay_environment == "S"))
speed_MS_evo_S <- lmer(Mean_movement_speed ~ Selection_T * Selection_S + (1|Selection_line),
                    traits_assay %>% filter(Assay_environment == "S"))

GR_MS_evo_T <- lmer(GR ~ Selection_T * Selection_S + (1|Selection_line),
                    traits_assay %>% filter(Assay_environment == "T"))
size_MS_evo_T <- lmer(Mean_cell_area ~ Selection_T * Selection_S + (1|Selection_line),
                      traits_assay %>% filter(Assay_environment == "T"))
speed_MS_evo_T <- lmer(Mean_movement_speed ~ Selection_T * Selection_S + (1|Selection_line),
                       traits_assay %>% filter(Assay_environment == "T"))

# Summary tables
summary(GR_MS_evo_S)
summary(size_MS_evo_S)
summary(speed_MS_evo_S)

summary(GR_MS_evo_T)
summary(size_MS_evo_T)
summary(speed_MS_evo_T)

# NOTE: multiple testing problem!
# A Bonferroni correction would set the significance level to 0.05/6 = 0.008, but since most of the
# significant results would remain significant, the interpretation also remains similar, only more
# nuanced.
