# ----------------------------------------------------------------------------------------------------------------------------------
# Manuscript figures
# ----------------------------------------------------------------------------------------------------------------------------------

# --------------------------------------------------------------------------------------------------
# Figure 1
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
# Data created in FR_ANOVA_simulations.R
load("FR data/Model predictions with simulated parameters/ANOVA-like fits/FR_predicted_fits_ANOVA_c.RData") # FR_predictions_c
FR_predictions_c$SEL_environment <- factor(FR_predictions_c$SEL_environment,
                                           levels = c("C", "S", "T", "ST"))


# Make a long version for plotting
FR_predictions_c_long <- FR_predictions_c %>%
  tidyr::pivot_longer(names_to = 'Consumer_presence', values_to = 'Resource_end_density', colnames(FR_predictions_c)[2:3])

FR_predictions_c_long$Consumer_presence <- factor(ifelse(FR_predictions_c_long$Consumer_presence == "Resource_end_density_without", 0, 1))

# ---------------------- Simulated parameter values ----------------------
n_sim <- 1000

load("FR data/Model predictions with simulated parameters/ANOVA-like fits/1000 simulations/FR_ANOVA_simulation_parameters_c.RData") # loads param_sim

# ---------------------- Loading simulations with varying parameters ----------------------
# Simulated FR values
load("FR data/Model predictions with simulated parameters/ANOVA-like fits/1000 simulations/FR_ANOVA_multiple_model_predictions_C_c.RData")
load("FR data/Model predictions with simulated parameters/ANOVA-like fits/1000 simulations/FR_ANOVA_multiple_model_predictions_S_c.RData")
load("FR data/Model predictions with simulated parameters/ANOVA-like fits/1000 simulations/FR_ANOVA_multiple_model_predictions_T_c.RData")
load("FR data/Model predictions with simulated parameters/ANOVA-like fits/1000 simulations/FR_ANOVA_multiple_model_predictions_ST_c.RData")

# Simulated resource growth
load("FR data/Model predictions with simulated parameters/ANOVA-like fits/1000 simulations/FR_ANOVA_multiple_LV_predictions_C_c.RData")
load("FR data/Model predictions with simulated parameters/ANOVA-like fits/1000 simulations/FR_ANOVA_multiple_LV_predictions_S_c.RData")
load("FR data/Model predictions with simulated parameters/ANOVA-like fits/1000 simulations/FR_ANOVA_multiple_LV_predictions_T_c.RData")
load("FR data/Model predictions with simulated parameters/ANOVA-like fits/1000 simulations/FR_ANOVA_multiple_LV_predictions_ST_c.RData")

# Plotting datasets
load("FR data/Model predictions with simulated parameters/ANOVA-like fits/1000 simulations/FR_ANOVA_sim_plot_data_C_c.RData")
load("FR data/Model predictions with simulated parameters/ANOVA-like fits/1000 simulations/FR_ANOVA_sim_plot_data_S_c.RData")
load("FR data/Model predictions with simulated parameters/ANOVA-like fits/1000 simulations/FR_ANOVA_sim_plot_data_T_c.RData")
load("FR data/Model predictions with simulated parameters/ANOVA-like fits/1000 simulations/FR_ANOVA_sim_plot_data_ST_c.RData")

# ---------------------- All CIs of consumption ----------------------
load("FR data/Model predictions with simulated parameters/ANOVA-like fits/1000 simulations/FR_CIs_ANOVA_like_c.RData")
# File created in FR_ANOVA_fits_simulations.Rmd
FRall_CIs_c$SEL_environment <- factor(FRall_CIs_c$SEL_environment,
                                      levels = c("C", "S", "T", "ST"))

# Adding info on LV fit uncertainty for each selection environment
FRall_CIs_c$mean_LV <- c(LV_mean_C_c,
                         LV_mean_S_c,
                         LV_mean_T_c,
                         LV_mean_ST_c)
FRall_CIs_c$lower_LV <- c(LV_lower_C_c,
                          LV_lower_S_c,
                          LV_lower_T_c,
                          LV_lower_ST_c)
FRall_CIs_c$upper_LV <- c(LV_upper_C_c,
                          LV_upper_S_c,
                          LV_upper_T_c,
                          LV_upper_ST_c)
FRall_CIs_c$Consumer_presence <- factor(1, levels = c(0, 1)) # Simulations were done with predators

# Adding one row for legend colors to work
FRall_CIs_c <- rbind(FRall_CIs_c, data.frame(mean = NA,
                                             lower = NA,
                                             upper = NA,
                                             mean_change_with = NA,
                                             lower_change_with = NA,
                                             upper_change_with = NA,
                                             N = 100,
                                             SEL_environment = "C",
                                             mean_LV = NA,
                                             lower_LV = NA,
                                             upper_LV = NA,
                                             Consumer_presence = factor(0, levels = c(0, 1))))


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










# ---------------------- Functional response fits with 95% CIs ----------------------
# Adding info on predator presence and one row without predators for plot colors
FR_predictions_c$Consumer_presence <- factor(1, levels = c(0, 1))
FR_predictions_c <- rbind(FR_predictions_c, data.frame(Resource_start_density = NA,
                                                       Resource_end_density_without = NA,
                                                       Resource_end_density_with = NA,
                                                       Resource_consumed = NA,           
                                                       SEL_environment = "C",
                                                       Consumer_presence = factor(0, levels = c(0, 1))))

P1 <- ggplot(FR_predictions_c, aes(x = Resource_start_density, y = Resource_consumed)) +
  geom_line(aes(color = Consumer_presence)) +
  geom_ribbon(aes(x = N, y = mean, ymin = lower, ymax = upper, fill = Consumer_presence), 
              FRall_CIs_c, alpha = 0.3) +
  theme_bw() +
  scale_color_manual(values = c("black", "red")) +
  scale_fill_manual(values = c("black", "red")) +
  facet_grid(cols = vars(SEL_environment)) +
  labs(x = latex2exp::TeX("Prey start density ($mL^{-1}$)"),
       y = latex2exp::TeX("Prey eaten ($predator^{-1} day^{-1}$)"),
       color = "Predator presence",
       fill = "Predator presence") +
  theme(strip.background = element_blank(),
        legend.position = "bottom") +
  theme(plot.margin = margin(t = 10, r = 5, b = 0, l = 20, "pt"))

# ---------------------- Density changes across prey start densities ----------------------
P2 <- ggplot(filter(FR_predictions_c_long, Resource_start_density >= 50), aes(x = (Resource_start_density), 
                                                                              y = Resource_end_density - Resource_start_density, 
                                                                              color = Consumer_presence)) +
  geom_line() +
  geom_point(aes(x = (Resource_start_density), y = Resource_end_density - Resource_start_density, color = Consumer_presence), 
             FR_all, size = 0.6) +
  labs(x = latex2exp::TeX("Prey start density ($mL^{-1}$)"),
       y = latex2exp::TeX("Change in prey density $(mL^{-1} day^{-1})$"),
       color = "Predator presence",
       fill = "Predator presence") +
  facet_grid(cols = vars(SEL_environment)) +
  theme_bw() +
  # scale_x_log10() + # uncomment this to have the log10-scaled plots for prey density change
  theme(strip.background = element_blank(),
        legend.position = "bottom") +
  geom_ribbon(aes(x = (N), y = mean, ymin = lower_LV - N, ymax = upper_LV - N, color = NULL, fill = Consumer_presence), 
              filter(FRall_CIs_c, N >= 50), alpha = 0.3) + 
  scale_color_manual(values = c("black", "red")) + 
  scale_fill_manual(values = c("black", "red")) +
  theme(plot.margin = margin(t = 10, r = 5, b = 0, l = 20, "pt"))

# ---------------------- Figure 1 unlogged ----------------------
ggpubr::ggarrange(P1, P2, nrow = 2, common.legend = T,
                  legend = "bottom", labels = "auto")
# Save with figure dimensions 5 inches x 10 inches landscape mode