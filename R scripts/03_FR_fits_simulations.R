# ----------------------------------------------------------------------------------------------------------------------------------
# Functional response fits and simulations of uncertainty around fit
# ----------------------------------------------------------------------------------------------------------------------------------
# 
# This script is divided in two parts:
# 
#   Part 1: best model fits
#     We use coefficient estimates from the best model fit obtained from the FR 
#     search grid (script 02_FR_search_grid.R) to get a FR curve for each FR 
#     measured (one for each prey selection environment)
# 
#   Part 2: simulations of uncertainty around fits
#     We use variance-covariance matrices to simulate 1000 realizations of model
#     coefficients within the uncertainty given by the best model fit. We use 
#     these to simulate 1000 FR fits for each prey selection environment and 
#     compute 95% CIs around the fits given by the best fit coefficients.
# 
# NOTE: this script was run on R version 4.3.2 (provided in folder
# "R versions/R version for search grids and ODE fitting/"). It will not run on R version 4.5.2
# due to package compatibility issues). The "odeintr" package is available on github for R version
# 4.3.2 and is unsupported for newer versions of R at the moment. A working version of the package
# can be installed with the following command:
# install.packages("R versions/R version for search grids and ODE fitting/odeintr_1.7.1_patched_clean.tar.gz", repos = NULL, type = "source")
# 
# ----------------------------------------------------------------------------------------------------------------------------------

library(tidyverse)
library(MASS)
library(broom)

#---------------------------------------------------------------------------------------------------
# Part 1: best model fits
#---------------------------------------------------------------------------------------------------

rm(list = ls())

# ---------------------- Loading data and functions ----------------------
source("R scripts/Function_FR_model_ANOVA_like.R")
load("Data/FR assay/01 - Raw data/Dataset_all.RData")
load("Data/FR assay/03 - FR search grid/FR_ANOVA_like_best_fit.RData")

FR_coefs <- generics::tidy(best_FR_fit)

# ---------------------- Best FR fit with prey C ----------------------
FRC_dd <- filter(FR_all, SEL_environment == "C")

FRC_predictions <- data.frame(Prey_start_density = 1:max(FRC_dd$Prey_start_density))

# Predicted end prey density without predators
FRC_predictions$Prey_end_density_without <- unlist(lapply(FRC_predictions$Prey_start_density, FUN = function(N) {
  Nt <- eq.odeint.general.pred(start.v = c(N, 0),
                            b = exp(coef(best_FR_fit)[["bC.log"]]),
                            h = coef(best_FR_fit)[["hC"]],
                            r = coef(best_FR_fit)[["rC"]],
                            K = exp(coef(best_FR_fit)[["KC.log"]]),
                            c = coef(best_FR_fit)[["c"]],
                            m = coef(best_FR_fit)[["mC"]],
                            Tt = 1,
                            timesteplength = 1/100)[101,2]
  return(Nt)
}))

# Predicted end prey density with 5 predators
FRC_predictions$Prey_end_density_with <- unlist(lapply(FRC_predictions$Prey_start_density, FUN = function(N) {
  Nt <- eq.odeint.general.pred(start.v = c(N, 5),
                                  b = exp(coef(best_FR_fit)[["bC.log"]]),
                                  h = coef(best_FR_fit)[["hC"]],
                                  r = coef(best_FR_fit)[["rC"]],
                                  K = exp(coef(best_FR_fit)[["KC.log"]]),
                                  c = coef(best_FR_fit)[["c"]],
                                  m = coef(best_FR_fit)[["mC"]],
                                  Tt = 1,
                                  timesteplength = 1/100)[101,2]
  return(Nt)
}))

# Checking if the predicted prey growth (lines) with and without predators matches the observed prey growth (points)
ggplot(FRC_predictions) +
  geom_line(aes(x = Prey_start_density, y = Prey_end_density_without - Prey_start_density), color = 1) +
  geom_line(aes(x = Prey_start_density, y = Prey_end_density_with - Prey_start_density), color = 2) +
  geom_point(aes(x = Prey_start_density, y = Prey_end_density - Prey_start_density,
                 color = Predator_presence), FRC_dd) +
  labs(x = latex2exp::TeX("Prey starting density ($mL^{-1}$)"),
       y = latex2exp::TeX("Change in prey density ($mL^{-1}$)"),
       color = "Predator presence",
       title = "Selection environment C") +
  scale_color_manual(values = 1:2) +
  theme_bw()

# Getting the predicted type 2 functional response curve from the best model fit
FRC_predictions <- FRC_predictions %>%
  mutate(Prey_consumed = exp(coef(best_FR_fit)[["bC.log"]]) * Prey_start_density / (1 + exp(coef(best_FR_fit)[["bC.log"]]) * coef(best_FR_fit)[["hC"]] * Prey_start_density))




# ---------------------- Best FR fit with prey T ----------------------
FRT_dd <- filter(FR_all, SEL_environment == "T")

FRT_predictions <- data.frame(Prey_start_density = 1:max(FRT_dd$Prey_start_density))

# Predicted end prey density without predators
FRT_predictions$Prey_end_density_without <- unlist(lapply(FRT_predictions$Prey_start_density, FUN = function(N) {
  Nt <- eq.odeint.general.pred(start.v = c(N, 0),
                                  b = exp(coef(best_FR_fit)[["bC.log"]] + coef(best_FR_fit)[["bT.log"]]),
                                  h = coef(best_FR_fit)[["hC"]] + coef(best_FR_fit)[["hT"]],
                                  r = coef(best_FR_fit)[["rT"]],
                                  K = exp(coef(best_FR_fit)[["KT.log"]]),
                                  c = coef(best_FR_fit)[["c"]],
                                  m = coef(best_FR_fit)[["mC"]] + coef(best_FR_fit)[["mT"]],
                                  Tt = 1,
                                  timesteplength = 1/100)[101,2]
  return(Nt)
}))

# Predicted end prey density with 5 predators
FRT_predictions$Prey_end_density_with <- unlist(lapply(FRT_predictions$Prey_start_density, FUN = function(N) {
  Nt <- eq.odeint.general.pred(start.v = c(N, 5),
                                     b = exp(coef(best_FR_fit)[["bC.log"]] + coef(best_FR_fit)[["bT.log"]]),
                                     h = coef(best_FR_fit)[["hC"]] + coef(best_FR_fit)[["hT"]],
                                     r = coef(best_FR_fit)[["rT"]],
                                     K = exp(coef(best_FR_fit)[["KT.log"]]),
                                     c = coef(best_FR_fit)[["c"]],
                                     m = coef(best_FR_fit)[["mC"]] + coef(best_FR_fit)[["mT"]],
                                  Tt = 1,
                                  timesteplength = 1/100)[101,2]
  return(Nt)
}))

# Checking if the predicted prey growth (lines) with and without predators matches the observed prey growth (points)
ggplot(FRT_predictions) +
  geom_line(aes(x = Prey_start_density, y = Prey_end_density_without - Prey_start_density), color = 1) +
  geom_line(aes(x = Prey_start_density, y = Prey_end_density_with - Prey_start_density), color = 2) +
  geom_point(aes(x = Prey_start_density, y = Prey_end_density - Prey_start_density,
                 color = Predator_presence), FRT_dd) +
  labs(x = latex2exp::TeX("Prey starting density ($mL^{-1}$)"),
       y = latex2exp::TeX("Change in prey density ($mL^{-1}$)"),
       color = "Predator presence",
       title = "Selection environment T") +
  scale_color_manual(values = 1:2) +
  theme_bw()

# Getting the predicted type 2 functional response curve from the best model fit
FRT_predictions <- FRT_predictions %>%
  mutate(Prey_consumed = exp(coef(best_FR_fit)[["bC.log"]] + coef(best_FR_fit)[["bT.log"]]) * Prey_start_density / (1 + exp(coef(best_FR_fit)[["bC.log"]] + coef(best_FR_fit)[["bT.log"]]) * (coef(best_FR_fit)[["hC"]] + coef(best_FR_fit)[["hT"]]) * Prey_start_density))


# ---------------------- Best FR fit with prey S ----------------------
FRS_dd <- filter(FR_all, SEL_environment == "S")

FRS_predictions <- data.frame(Prey_start_density = 1:max(FRS_dd$Prey_start_density))

# Predicted end prey density without predators
FRS_predictions$Prey_end_density_without <- unlist(lapply(FRS_predictions$Prey_start_density, FUN = function(N) {
  Nt <- eq.odeint.general.pred(start.v = c(N, 0),
                                     b = exp(coef(best_FR_fit)[["bC.log"]] + coef(best_FR_fit)[["bS.log"]]),
                                     h = coef(best_FR_fit)[["hC"]] + coef(best_FR_fit)[["hS"]],
                                     r = coef(best_FR_fit)[["rS"]],
                                     K = exp(coef(best_FR_fit)[["KS.log"]]),
                                     c = coef(best_FR_fit)[["c"]],
                                     m = coef(best_FR_fit)[["mC"]] + coef(best_FR_fit)[["mS"]],
                                  Tt = 1,
                                  timesteplength = 1/100)[101,2]
  return(Nt)
}))

# Predicted end prey density with 5 predators
FRS_predictions$Prey_end_density_with <- unlist(lapply(FRS_predictions$Prey_start_density, FUN = function(N) {
  Nt <- eq.odeint.general.pred(start.v = c(N, 5),
                                     b = exp(coef(best_FR_fit)[["bC.log"]] + coef(best_FR_fit)[["bS.log"]]),
                                     h = coef(best_FR_fit)[["hC"]] + coef(best_FR_fit)[["hS"]],
                                     r = coef(best_FR_fit)[["rS"]],
                                     K = exp(coef(best_FR_fit)[["KS.log"]]),
                                     c = coef(best_FR_fit)[["c"]],
                                     m = coef(best_FR_fit)[["mC"]] + coef(best_FR_fit)[["mS"]],
                                  Tt = 1,
                                  timesteplength = 1/100)[101,2]
  return(Nt)
}))

# Checking if the predicted prey growth (lines) with and without predators matches the observed prey growth (points)
ggplot(FRS_predictions) +
  geom_line(aes(x = Prey_start_density, y = Prey_end_density_without - Prey_start_density), color = 1) +
  geom_line(aes(x = Prey_start_density, y = Prey_end_density_with - Prey_start_density), color = 2) +
  geom_point(aes(x = Prey_start_density, y = Prey_end_density - Prey_start_density,
                 color = Predator_presence), FRS_dd) +
  labs(x = latex2exp::TeX("Prey starting density ($mL^{-1}$)"),
       y = latex2exp::TeX("Change in prey density ($mL^{-1}$)"),
       color = "Predator presence",
       title = "Selection environment S") +
  scale_color_manual(values = 1:2) +
  theme_bw()

# Getting the predicted type 2 functional response curve from the best model fit
FRS_predictions <- FRS_predictions %>%
  mutate(Prey_consumed = exp(coef(best_FR_fit)[["bC.log"]] + coef(best_FR_fit)[["bS.log"]]) * Prey_start_density / (1 + exp(coef(best_FR_fit)[["bC.log"]] + coef(best_FR_fit)[["bS.log"]]) * (coef(best_FR_fit)[["hC"]] + coef(best_FR_fit)[["hS"]]) * Prey_start_density))


# ---------------------- Best FR fit with prey ST ----------------------
FRST_dd <- filter(FR_all, SEL_environment == "ST")

FRST_predictions <- data.frame(Prey_start_density = 1:max(FRST_dd$Prey_start_density))

# Predicted end prey density without predators
FRST_predictions$Prey_end_density_without <- unlist(lapply(FRST_predictions$Prey_start_density, FUN = function(N) {
  Nt <- eq.odeint.general.pred(start.v = c(N, 0),
                                     b = exp(coef(best_FR_fit)[["bC.log"]] + coef(best_FR_fit)[["bS.log"]] + coef(best_FR_fit)[["bT.log"]] + coef(best_FR_fit)[["bST.log"]]),
                                     h = coef(best_FR_fit)[["hC"]] + coef(best_FR_fit)[["hS"]] + coef(best_FR_fit)[["hT"]] + coef(best_FR_fit)[["hST"]],
                                     r = coef(best_FR_fit)[["rST"]],
                                     K = exp(coef(best_FR_fit)[["KST.log"]]),
                                     c = coef(best_FR_fit)[["c"]],
                                     m = coef(best_FR_fit)[["mC"]] + coef(best_FR_fit)[["mS"]] + coef(best_FR_fit)[["mT"]] + coef(best_FR_fit)[["mST"]],
                                  Tt = 1,
                                  timesteplength = 1/100)[101,2]
  return(Nt)
}))

# Predicted end prey density with 5 predators
FRST_predictions$Prey_end_density_with <- unlist(lapply(FRST_predictions$Prey_start_density, FUN = function(N) {
  Nt <- eq.odeint.general.pred(start.v = c(N, 5),
                                     b = exp(coef(best_FR_fit)[["bC.log"]] + coef(best_FR_fit)[["bS.log"]] + coef(best_FR_fit)[["bT.log"]] + coef(best_FR_fit)[["bST.log"]]),
                                     h = coef(best_FR_fit)[["hC"]] + coef(best_FR_fit)[["hS"]] + coef(best_FR_fit)[["hT"]] + coef(best_FR_fit)[["hST"]],
                                     r = coef(best_FR_fit)[["rST"]],
                                     K = exp(coef(best_FR_fit)[["KST.log"]]),
                                     c = coef(best_FR_fit)[["c"]],
                                     m = coef(best_FR_fit)[["mC"]] + coef(best_FR_fit)[["mS"]] + coef(best_FR_fit)[["mT"]] + coef(best_FR_fit)[["mST"]],
                                  Tt = 1,
                                  timesteplength = 1/100)[101,2]
  return(Nt)
}))

# Checking if the predicted prey growth (lines) with and without predators matches the observed prey growth (points)
ggplot(FRST_predictions) +
  geom_line(aes(x = Prey_start_density, y = Prey_end_density_without - Prey_start_density), color = 1) +
  geom_line(aes(x = Prey_start_density, y = Prey_end_density_with - Prey_start_density), color = 2) +
  geom_point(aes(x = Prey_start_density, y = Prey_end_density - Prey_start_density,
                 color = Predator_presence), FRST_dd) +
  labs(x = latex2exp::TeX("Prey starting density ($mL^{-1}$)"),
       y = latex2exp::TeX("Change in prey density ($mL^{-1}$)"),
       color = "Predator presence",
       title = "Selection environment ST") +
  scale_color_manual(values = 1:2) +
  theme_bw()

# Getting the predicted type 2 functional response curve from the best model fit
FRST_predictions <- FRST_predictions %>%
  mutate(Prey_consumed = exp(coef(best_FR_fit)[["bC.log"]] + coef(best_FR_fit)[["bS.log"]] + coef(best_FR_fit)[["bT.log"]] + coef(best_FR_fit)[["bST.log"]]) * Prey_start_density / (1 + exp(coef(best_FR_fit)[["bC.log"]] + coef(best_FR_fit)[["bS.log"]] + coef(best_FR_fit)[["bT.log"]] + coef(best_FR_fit)[["bST.log"]]) * (coef(best_FR_fit)[["hC"]] + coef(best_FR_fit)[["hS"]] + coef(best_FR_fit)[["hT"]] + coef(best_FR_fit)[["hST"]]) * Prey_start_density))

# ----------------- Merging the four datasets -----------------
FRC_predictions$SEL_environment = "C"
FRS_predictions$SEL_environment = "S"
FRT_predictions$SEL_environment = "T"
FRST_predictions$SEL_environment = "ST"

FR_predictions <- rbind(FRC_predictions,
                        FRS_predictions,
                        FRT_predictions,
                        FRST_predictions)

#---------------------------------------------------------------------------------------------------
# Part 2: simulations of uncertainty around fits
#---------------------------------------------------------------------------------------------------

# ---------------------- Cleaning workspace ----------------------
rm(list = c("FRC_dd", "FRS_dd", "FRT_dd", "FRST_dd",
            "FRC_predictions", "FRS_predictions", "FRT_predictions", "FRST_predictions"))

# ---------------------- Simulating 1000 realizations ----------------------
set.seed(42)

# Simulate 1000 realizations of FR parameter combinations
n_sim <- 1000
param_sim <- mvrnorm(n_sim, mu = FR_coefs$estimate, Sigma = vcov(best_FR_fit)) %>%
  data.frame()

# Adding the r and K information for each prey selection environment
param_sim$rC <- coef(best_FR_fit)[["rC"]]
param_sim$rS <- coef(best_FR_fit)[["rS"]]
param_sim$rT <- coef(best_FR_fit)[["rT"]]
param_sim$rST <- coef(best_FR_fit)[["rST"]]
param_sim$KC.log <- coef(best_FR_fit)[["KC.log"]]
param_sim$KS.log <- coef(best_FR_fit)[["KS.log"]]
param_sim$KT.log <- coef(best_FR_fit)[["KT.log"]]
param_sim$KST.log <- coef(best_FR_fit)[["KST.log"]]

param_sim <- param_sim %>%
  as.matrix()

# Getting unique x-axes for each prey selection environment
N_C <- 1:max(FR_predictions$Prey_start_density[FR_predictions$SEL_environment == "C"])
N_S <- 1:max(FR_predictions$Prey_start_density[FR_predictions$SEL_environment == "S"])
N_T <- 1:max(FR_predictions$Prey_start_density[FR_predictions$SEL_environment == "T"])
N_ST <- 1:max(FR_predictions$Prey_start_density[FR_predictions$SEL_environment == "ST"])

# -------------------- Getting FR fits for each realization --------------------
# For each prey selection environment, we compute the number of prey eaten across 
# prey densities and realized parameter combinations.
# Output will be a matrix where each column is the estimated FR curve for
# one of the parameter combinations (1000 columns x N rows).
predictions_C <- apply(param_sim, 1, function(p) {
  Ndead <- exp(p[["bC.log"]]) * N_C / (1 + exp(p[["bC.log"]]) * p[["hC"]] * N_C)
})
predictions_S <- apply(param_sim, 1, function(p) {
  Ndead <- exp(p[["bC.log"]] + p[["bS.log"]]) * N_S / (1 + exp(p[["bC.log"]] + p[["bS.log"]]) * (p[["hC"]] + p[["hS"]]) * N_S)
})
predictions_T <- apply(param_sim, 1, function(p) {
  Ndead <- exp(p[["bC.log"]] + p[["bT.log"]]) * N_T / (1 + exp(p[["bC.log"]] + p[["bT.log"]]) * (p[["hC"]] + p[["hT"]]) * N_T)
})
predictions_ST <- apply(param_sim, 1, function(p) {
  Ndead <- exp(p[["bC.log"]] + p[["bS.log"]] + p[["bT.log"]] + p[["bST.log"]]) * N_ST / (1 + exp(p[["bC.log"]] + p[["bS.log"]] + p[["bT.log"]] + p[["bST.log"]]) * (p[["hC"]] + p[["hS"]] + p[["hT"]] + p[["hST"]]) * N_ST)
})

# ------------------- Computing 95% CIs -------------------
# --- Uncomment the following code to create datasets for plotting 1000 fits ---
# predictions_dd_C <- data.frame(predictions_C) %>%
#   pivot_longer(cols = 1:n_sim, names_to = "Simulation_nr", values_to = "Ndead")
# predictions_dd_C$N <- rep(N_C, each = n_sim)
# predictions_dd_C$b.log.estimate <- rep(param_sim[,"bC.log"], times = length(N_C))
# predictions_dd_C$h.estimate <- rep(param_sim[,"hC"], times = length(N_C))
# predictions_dd_C$m.estimate <- rep(param_sim[,"mC"], times = length(N_C))
# 
# predictions_dd_S <- data.frame(predictions_S) %>%
#   pivot_longer(cols = 1:n_sim, names_to = "Simulation_nr", values_to = "Ndead")
# predictions_dd_S$N <- rep(N_S, each = n_sim)
# predictions_dd_S$b.log.estimate <- rep(param_sim[,"bC.log"] + param_sim[,"bS.log"], times = length(N_S))
# predictions_dd_S$h.estimate <- rep(param_sim[,"hC"] + param_sim[,"hS"], times = length(N_S))
# predictions_dd_S$m.estimate <- rep(param_sim[,"mC"] + param_sim[,"mS"], times = length(N_S))
# 
# predictions_dd_T <- data.frame(predictions_T) %>%
#   pivot_longer(cols = 1:n_sim, names_to = "Simulation_nr", values_to = "Ndead")
# predictions_dd_T$N <- rep(N_T, each = n_sim)
# predictions_dd_T$b.log.estimate <- rep(param_sim[,"bC.log"] + param_sim[,"bT.log"], times = length(N_T))
# predictions_dd_T$h.estimate <- rep(param_sim[,"hC"] + param_sim[,"hT"], times = length(N_T))
# predictions_dd_T$m.estimate <- rep(param_sim[,"mC"] + param_sim[,"mT"], times = length(N_T))
# 
# predictions_dd_ST <- data.frame(predictions_ST) %>%
#   pivot_longer(cols = 1:n_sim, names_to = "Simulation_nr", values_to = "Ndead")
# predictions_dd_ST$N <- rep(N_ST, each = n_sim)
# predictions_dd_ST$b.log.estimate <- rep(param_sim[,"bC.log"] + param_sim[,"bS.log"] + param_sim[,"bT.log"] + param_sim[,"bST.log"], times = length(N_ST))
# predictions_dd_ST$h.estimate <- rep(param_sim[,"hC"] + param_sim[,"hS"] + param_sim[,"hT"] + param_sim[,"hST"], times = length(N_ST))
# predictions_dd_ST$m.estimate <- rep(param_sim[,"mC"] + param_sim[,"mS"] + param_sim[,"mT"] + param_sim[,"mST"], times = length(N_ST))

# ------------------ prey C ------------------
# --- Uncomment the following code to plot 1000 fits ---
# ggplot(predictions_dd_C, aes(x = N, y = Ndead, group = Simulation_nr)) +
#   geom_line(linewidth = .5, alpha = .3)

# Compute confidence intervals
response_lower_C <- apply(predictions_C, 1, quantile, probs = 0.025, na.rm = T)
response_upper_C <- apply(predictions_C, 1, quantile, probs = 0.975, na.rm = T)

# ------------------ prey S ------------------
# --- Uncomment the following code to plot 1000 fits ---
# ggplot(predictions_dd_S, aes(x = N, y = Ndead, group = Simulation_nr)) +
#   geom_line(linewidth = .5, alpha = .3)

# Compute confidence intervals
response_lower_S <- apply(predictions_S, 1, quantile, probs = 0.025, na.rm = T)
response_upper_S <- apply(predictions_S, 1, quantile, probs = 0.975, na.rm = T)

# ------------------ prey T ------------------
# --- Uncomment the following code to plot 1000 fits ---
# ggplot(predictions_dd_T, aes(x = N, y = Ndead, group = Simulation_nr)) +
#   geom_line(linewidth = .5, alpha = .3)

# Compute confidence intervals
response_lower_T <- apply(predictions_T, 1, quantile, probs = 0.025, na.rm = T)
response_upper_T <- apply(predictions_T, 1, quantile, probs = 0.975, na.rm = T)

# ------------------ prey ST ------------------
# --- Uncomment the following code to plot 1000 fits ---
# ggplot(predictions_dd_ST, aes(x = N, y = Ndead, group = Simulation_nr)) +
#   geom_line(linewidth = .5, alpha = .3)

# Compute confidence intervals
response_lower_ST <- apply(predictions_ST, 1, quantile, probs = 0.025, na.rm = T)
response_upper_ST <- apply(predictions_ST, 1, quantile, probs = 0.975, na.rm = T)


# ----------------- Adding 95% CIs to dataset from part 1 -----------------
FRC_CIs <- data.frame(lower = response_lower_C,
                      upper = response_upper_C,
                      Prey_start_density = N_C,
                      SEL_environment = factor("C", levels = c("C", "S", "T", "ST")))
FRS_CIs <- data.frame(lower = response_lower_S,
                      upper = response_upper_S,
                      Prey_start_density = N_S,
                      SEL_environment = factor("S", levels = c("C", "S", "T", "ST")))
FRT_CIs <- data.frame(lower = response_lower_T,
                      upper = response_upper_T,
                      Prey_start_density = N_T,
                      SEL_environment = factor("T", levels = c("C", "S", "T", "ST")))
FRST_CIs <- data.frame(lower = response_lower_ST,
                       upper = response_upper_ST,
                       Prey_start_density = N_ST,
                       SEL_environment = factor("ST", levels = c("C", "S", "T", "ST")))

FRall_CIs <- rbind(FRC_CIs,
                   FRS_CIs,
                   FRT_CIs,
                   FRST_CIs)

FR_predictions <- merge(FR_predictions, FRall_CIs, by = c("SEL_environment", "Prey_start_density")) %>%
  arrange(SEL_environment, Prey_start_density) %>%
  dplyr::select(-Prey_end_density_with, -Prey_end_density_without) %>%
  rename(CI_lower = lower,
         CI_upper = upper)

# save(FR_predictions, file = "Data/FR assay/04 - FR fits and simulations/Best_fit_curves_and_95%CIs.RData")




