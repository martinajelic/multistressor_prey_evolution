# ----------------------------------------------------------------------------------------------------------------------------------
# Manuscript figures
# ----------------------------------------------------------------------------------------------------------------------------------

# --------------------------------------------------------------------------------------------------
# Figure 1: FR parameters and fits
# --------------------------------------------------------------------------------------------------

rm(list = ls())

library(tidyverse)
library(generics)
library(broom)

# ---------------------- Loading FR dataset ----------------------
load("Data/FR assay/01 - Raw data/Dataset_all.RData")


# ---------------------- Loading best model fits ----------------------
load("Data/FR assay/03 - FR search grid/FR_ANOVA_like_best_fit.RData")

vcov_matrix <- vcov(best_FR_fit)

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

# ---------------------- Loading predicted values from best model fits ----------------------
# Dataset created in FR_fits_simulations.R
load("Data/FR assay/04 - FR fits and simulations/Best_fit_curves_and_95%CIs.RData") # FR_predictions
FR_predictions$SEL_environment <- factor(FR_predictions$SEL_environment,
                                           levels = c("C", "S", "T", "ST"))

# ---------------------- Estimated parameters and uncertainty ----------------------
param_df <- rbind(FRC_coefs %>% dplyr::select(term, estimate, std.error) %>% mutate(SEL_environment = "C"),
                  FRS_coefs %>% dplyr::select(term, estimate, std.error) %>% mutate(SEL_environment = "S"),
                  FRT_coefs %>% dplyr::select(term, estimate, std.error) %>% mutate(SEL_environment = "T"),
                  FRST_coefs %>% dplyr::select(term, estimate, std.error) %>% mutate(SEL_environment = "ST")) %>%
  dplyr::rename(Parameter = term,
                Estimate = estimate,
                StdError = std.error)

param_df$Parameter[grep("b", param_df$Parameter)] <- "b.log"
param_df$Parameter[grep("h", param_df$Parameter)] <- "h"
param_df$Parameter[grep("m", param_df$Parameter)] <- "m"

param_df$CI_lower <- param_df$Estimate - 1.96 * param_df$StdError
param_df$CI_upper <- param_df$Estimate + 1.96 * param_df$StdError

param_df$StdUpper <- param_df$Estimate + param_df$StdError
param_df$StdLower <- param_df$Estimate - param_df$StdError

# Changing the b.log facet to exp(b.log)
param_df_exp <- param_df

# Exp the b.log values
param_df_exp$Estimate <- ifelse(param_df_exp$Parameter == "b.log", exp(param_df_exp$Estimate), param_df_exp$Estimate)

# Change b.log label to b
param_df_exp$Parameter <- ifelse(param_df_exp$Parameter == "b.log", "b", param_df_exp$Parameter)

# Compute upper and lower values of Std. error bars
param_df_exp$StdUpper <- ifelse(param_df_exp$Parameter == "b",
                                exp(param_df_exp$StdUpper),
                                param_df_exp$StdUpper)
param_df_exp$StdLower <- ifelse(param_df_exp$Parameter == "b",
                                exp(param_df_exp$StdLower),
                                param_df_exp$StdLower)

# Compute upper and lower values of CIs
param_df_exp$CI_lower <- ifelse(param_df_exp$Parameter == "b",
                                exp(param_df_exp$CI_lower),
                                param_df_exp$CI_lower)
param_df_exp$CI_upper <- ifelse(param_df_exp$Parameter == "b",
                                exp(param_df_exp$CI_upper),
                                param_df_exp$CI_upper)

# Relevel SEL_environment
param_df_exp$SEL_environment <- factor(param_df_exp$SEL_environment,
                                       levels = c("C", "S", "T", "ST"))

# ---------------------- Figure 1 ----------------------
Pb <- ggplot(filter(param_df_exp, Parameter == "b"), aes(x = SEL_environment, y = Estimate, color = SEL_environment)) +
  geom_point(size = 3, position = position_dodge(width = 0.7)) +
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper), width = 0.2, size = 0.7,
                position = position_dodge(width = 0.7)) +
  theme_bw() +
  scale_color_manual(values = c("#716868", "#E69F00", "#56B4E9", "#009E73")) +
  labs(x = "Prey selection environment",
       y = latex2exp::TeX("Space clearance rate $\\textit{b}$ (mL $day^{-1}$)"),
       color = "Selection environment") +
  theme(legend.position = "none") +
  geom_hline(yintercept = 0, linetype = 2, color = "gray") +
  theme(strip.background = element_blank()) +
  theme(plot.margin = margin(t = 10, r = 5, b = 0, l = 23, "pt"))

Ph <- ggplot(filter(param_df_exp, Parameter == "h"), aes(x = SEL_environment, y = Estimate, color = SEL_environment)) +
  geom_point(size = 3, position = position_dodge(width = 0.7)) +
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper), width = 0.2, size = 0.7,
                position = position_dodge(width = 0.7)) +
  theme_bw() +
  scale_color_manual(values = c("#716868", "#E69F00", "#56B4E9", "#009E73")) +
  labs(x = "Prey selection environment",
       y = latex2exp::TeX("Handling time $\\textit{h}$ $(day)$"),
       color = "Selection environment") +
  theme(legend.position = "none") +
  geom_hline(yintercept = 0, linetype = 2, color = "gray") +
  theme(strip.background = element_blank()) +
  theme(plot.margin = margin(t = 10, r = 5, b = 0, l = 23, "pt"))

Pm <- ggplot(filter(param_df_exp, Parameter == "m"), aes(x = SEL_environment, y = Estimate, color = SEL_environment)) +
  geom_point(size = 3, position = position_dodge(width = 0.7)) +
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper), width = 0.3, size = 0.7,
                position = position_dodge(width = 0.7)) +
  theme_bw() +
  scale_color_manual(values = c("#716868", "#E69F00", "#56B4E9", "#009E73")) +
  labs(x = "Prey selection environment",
       y = latex2exp::TeX("Predator mortality $\\textit{m}$ $(day^{-1})$"),
       color = "Selection environment") +
  theme(legend.position = "none") +
  geom_hline(yintercept = 0, linetype = 2, color = "gray") +
  theme(strip.background = element_blank()) +
  theme(plot.margin = margin(t = 10, r = 5, b = 0, l = 23, "pt"))

PFR_CS <- ggplot(filter(FR_predictions, Prey_start_density < 1000,
                        SEL_environment %in% c("C", "S")), aes(x = Prey_start_density, y = Prey_consumed, color = SEL_environment, group = SEL_environment)) +
  geom_line(size = 1.1) +
  geom_ribbon(aes(x = Prey_start_density, y = Prey_consumed, ymin = CI_lower, ymax = CI_upper, fill = SEL_environment, color = NULL), alpha = 0.15) +
  theme_bw() +
  scale_color_manual(values = c("#716868", "#E69F00")) +
  scale_fill_manual(values = c("#716868", "#E69F00")) +
  labs(x = latex2exp::TeX("Prey start density ($mL^{-1}$)"),
       y = latex2exp::TeX("Prey eaten ($predator^{-1} day^{-1}$)"),
       color = "Selection environment",
       fill = "Selection environment") +
  theme(strip.background = element_blank(),
        legend.position = "bottom") +
  theme(plot.margin = margin(t = 10, r = 5, b = 0, l = 20, "pt")) +
  ylim(0, 125)


PFR_CT <- ggplot(filter(FR_predictions, Prey_start_density < 1000,
                        SEL_environment %in% c("C", "T")), aes(x = Prey_start_density, y = Prey_consumed, color = SEL_environment, group = SEL_environment)) +
  geom_line(size = 1.1) +
  geom_ribbon(aes(x = Prey_start_density, y = Prey_consumed, ymin = CI_lower, ymax = CI_upper, fill = SEL_environment, color = NULL), alpha = 0.15) +
  theme_bw() +
  scale_color_manual(values = c("#716868", "#56B4E9")) +
  scale_fill_manual(values = c("#716868", "#56B4E9")) +
  labs(x = latex2exp::TeX("Prey start density ($mL^{-1}$)"),
       y = latex2exp::TeX("Prey eaten ($predator^{-1} day^{-1}$)"),
       color = "Selection environment",
       fill = "Selection environment") +
  theme(strip.background = element_blank(),
        legend.position = "bottom") +
  theme(plot.margin = margin(t = 10, r = 5, b = 0, l = 20, "pt")) +
  ylim(0, 125)

PFR_CST <- ggplot(filter(FR_predictions, Prey_start_density < 1000,
                         SEL_environment %in% c("C", "ST")), aes(x = Prey_start_density, y = Prey_consumed, color = SEL_environment, group = SEL_environment)) +
  geom_line(size = 1.1) +
  geom_ribbon(aes(x = Prey_start_density, y = Prey_consumed, ymin = CI_lower, ymax = CI_upper, fill = SEL_environment, color = NULL), alpha = 0.15) +
  theme_bw() +
  scale_color_manual(values = c("#716868", "#009E73")) +
  scale_fill_manual(values = c("#716868", "#009E73")) +
  labs(x = latex2exp::TeX("Prey start density ($mL^{-1}$)"),
       y = latex2exp::TeX("Prey eaten ($predator^{-1} day^{-1}$)"),
       color = "Selection environment",
       fill = "Selection environment") +
  theme(strip.background = element_blank(),
        legend.position = "bottom") +
  theme(plot.margin = margin(t = 10, r = 5, b = 0, l = 20, "pt")) +
  ylim(0, 125)

ggpubr::ggarrange(Pb, Ph, Pm, PFR_CS, PFR_CT, PFR_CST, ncol = 3, nrow = 2, labels = "auto", common.legend = T)
# Save as pdf with figure dimensions 10 inches x 6 inches landscape mode




# --------------------------------------------------------------------------------------------------
# Figure 2: prey traits and growth rate
# --------------------------------------------------------------------------------------------------

rm(list = ls())

load("Data/Traits_assay.RData")

P_GR <- ggplot() +
  geom_boxplot(data = traits_assay, aes(x = Assay_environment, y = GR,
                                                    color = Selection_environment, fill = Selection_environment),
               alpha = 0.3, linewidth = 0.2) +
  geom_point(data = traits_assay, aes(x = Assay_environment, y = GR,
                                                  color = Selection_environment, fill = Selection_environment, shape = Biorep),
             alpha = 1, size = 1, position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.75)) +
  theme_bw() +
  labs(x = "Assay environment",
       y = latex2exp::TeX("Growth rate $ln(N_{t}/N_{0})/\\Delta t$"),
       color = "Selection environment",
       fill = "Selection environment",
       shape = "Replicate") +
  scale_color_manual(values = c("#000000", "#E69F00", "#56B4E9", "#009E73")) +
  scale_fill_manual(values = c("#000000", "#E69F00", "#56B4E9", "#009E73")) +
  theme(legend.position = "none") +
  geom_hline(yintercept = 0, linetype = 2, color = "gray") +
  theme(strip.background = element_blank())

P_size <- ggplot() +
  geom_boxplot(data = traits_assay, aes(x = Assay_environment, y = Mean_cell_area,
                                                    color = Selection_environment, fill = Selection_environment),
               alpha = 0.3, linewidth = 0.2) +
  geom_point(data = traits_assay, aes(x = Assay_environment, y = Mean_cell_area,
                                                  color = Selection_environment, fill = Selection_environment, shape = Biorep),
             alpha = 1, size = 1, position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.75)) +
  theme_bw() +
  labs(x = "Assay environment",
       y = latex2exp::TeX("Cell area $(µm^2)$"),
       color = "Selection environment",
       fill = "Selection environment",
       shape = "Replicate") +
  scale_color_manual(values = c("#000000", "#E69F00", "#56B4E9", "#009E73")) +
  scale_fill_manual(values = c("#000000", "#E69F00", "#56B4E9", "#009E73")) +
  theme(legend.position = "none") +
  theme(strip.background = element_blank())

P_speed <- ggplot() +
  geom_boxplot(data = traits_assay, aes(x = Assay_environment, y = Mean_movement_speed,
                                                    color = Selection_environment, fill = Selection_environment),
               alpha = 0.3, linewidth = 0.2) +
  geom_point(data = traits_assay, aes(x = Assay_environment, y = Mean_movement_speed,
                                                  color = Selection_environment, fill = Selection_environment, shape = Biorep),
             alpha = 1, size = 1, position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.75)) +
  theme_bw() +
  labs(x = "Assay environment",
       y = latex2exp::TeX("Movement speed ($µm$ $sec^{-1}$)"),
       color = "Selection environment",
       fill = "Selection environment",
       shape = "Replicate") +
  scale_color_manual(values = c("#000000", "#E69F00", "#56B4E9", "#009E73")) +
  scale_fill_manual(values = c("#000000", "#E69F00", "#56B4E9", "#009E73")) +
  theme(legend.position = "none") +
  theme(strip.background = element_blank())

ggpubr::ggarrange(P_GR, P_size, P_speed,
                  ncol = 3, nrow = 1, common.legend = T, labels = "auto")
# Save as pdf 3.25 x 10 inches landscape mode



