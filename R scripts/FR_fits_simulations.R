# ----------------------------------------------------------------------------------------------------------------------------------
# Functional response fits and simulations of uncertainty around fit
# ----------------------------------------------------------------------------------------------------------------------------------

rm(list = ls())

library(tidyverse)
library(MASS)
library(broom)

source("R scripts/Function_FR_model_ANOVA_like.R")

load("FR data/Dataset_all.RData")

# ---------------------- Loading best model fits ----------------------
load("FR data/FR search grid/ANOVA-like fits/Full_c_best_typeII_fit.RData")

FR_coefs <- generics::tidy(best.mod.fit.full.c)

# FRC_coefs <- FR_coefs[grep("C", FR_coefs$term),]
# FRS_coefs <- FR_coefs[grep("C", FR_coefs$term) + 1,]
# FRT_coefs <- FR_coefs[grep("C", FR_coefs$term) + 2,]
# FRST_coefs <- FR_coefs[grep("ST", FR_coefs$term),]
# 
# FRST_coefs$estimate <- FRST_coefs$estimate + FRS_coefs$estimate + FRT_coefs$estimate + FRC_coefs$estimate
# FRST_coefs$std.error <- sqrt(FRST_coefs$std.error^2 + FRS_coefs$std.error^2 + FRT_coefs$std.error^2 + FRC_coefs$std.error^2)
# 
# FRS_coefs$estimate <- FRS_coefs$estimate + FRC_coefs$estimate
# FRS_coefs$std.error <- sqrt(FRS_coefs$std.error^2 + FRC_coefs$std.error^2)
# 
# 
# FRT_coefs$estimate <- FRT_coefs$estimate + FRC_coefs$estimate
# FRT_coefs$std.error <- sqrt(FRT_coefs$std.error^2 + FRC_coefs$std.error^2)





# ---------------------- Loading predicted fits from best models ----------------------
load("FR data/Model predictions with simulated parameters/ANOVA-like fits/FR_predicted_fits_ANOVA_c.RData")

# This dataset was created with the commented code below

# ####### Simulating FR fits using best model estimates for parameters
# # ---------------------- FRC ----------------------
# FRC_dd <- filter(FR_all, SEL_environment == "C")
# 
# ggplot(FRC_dd, aes(x = Resource_start_density, y = Resource_end_density - Resource_start_density, color = Consumer_presence)) +
#   geom_point()
# 
# FRC_predictions <- data.frame(Resource_start_density = 1:max(FRC_dd$Resource_start_density))
# 
# 
# # Predicted end density without consumers
# FRC_predictions$Resource_end_density_without <- unlist(lapply(FRC_predictions$Resource_start_density, FUN = function(N) {
#   Nt <- eq.odeint.general.pred.t2.mP(start.v = c(N, 0),
#                             b = exp(coef(best.mod.fit.full.c)[["bC.log"]]),
#                             h = coef(best.mod.fit.full.c)[["hC"]],
#                             r = coef(best.mod.fit.full.c)[["rC"]],
#                             K = exp(coef(best.mod.fit.full.c)[["KC.log"]]),
#                             c = coef(best.mod.fit.full.c)[["c"]],
#                             m = coef(best.mod.fit.full.c)[["mC"]],
#                             Tt = 1,
#                             timesteplength = 1/100)[101,2]
#   return(Nt)
# }))
# 
# # Predicted end density with 5 consumers
# FRC_predictions$Resource_end_density_with <- unlist(lapply(FRC_predictions$Resource_start_density, FUN = function(N) {
#   Nt <- eq.odeint.general.pred.t2.mP(start.v = c(N, 5),
#                                   b = exp(coef(best.mod.fit.full.c)[["bC.log"]]),
#                                   h = coef(best.mod.fit.full.c)[["hC"]],
#                                   r = coef(best.mod.fit.full.c)[["rC"]],
#                                   K = exp(coef(best.mod.fit.full.c)[["KC.log"]]),
#                                   c = coef(best.mod.fit.full.c)[["c"]],
#                                   m = coef(best.mod.fit.full.c)[["mC"]],
#                                   Tt = 1,
#                                   timesteplength = 1/100)[101,2]
#   return(Nt)
# }))
# 
# 
# GRC_without <- ggplot(FRC_predictions) +
#   geom_line(aes(x = Resource_start_density, y = Resource_end_density_without - Resource_start_density)) +
#   # geom_line(aes(x = Resource_start_density, y = Resource_end_density_with - Resource_start_density), color = 2) +
#   geom_point(aes(x = Resource_start_density, y = Resource_end_density - Resource_start_density), filter(FRC_dd, Consumer_presence == 0)) +
#   labs(x = "Resource starting density (mL-1)",
#        y = "End - start resource density",
#        title = "Selection environment C - Growth without consumers") +
#   theme_bw()
# 
# GRC_with <- ggplot(FRC_predictions) +
#   geom_line(aes(x = Resource_start_density, y = Resource_end_density_with - Resource_start_density)) +
#   geom_point(aes(x = Resource_start_density, y = Resource_end_density - Resource_start_density), filter(FRC_dd, Consumer_presence == 1)) +
#   labs(x = "Resource starting density (mL-1)",
#        y = "End - start resource density",
#        title = "Selection environment C - Growth with consumers") +
#   theme_bw()
# 
# ggpubr::ggarrange(GRC_without, GRC_with, ncol = 2)
# 
# 
# ggplot(FRC_predictions) +
#   geom_line(aes(x = Resource_start_density, y = Resource_end_density_without - Resource_start_density), color = 1) +
#   geom_line(aes(x = Resource_start_density, y = Resource_end_density_with - Resource_start_density), color = 2) +
#   geom_point(aes(x = Resource_start_density, y = Resource_end_density - Resource_start_density,
#                  color = Consumer_presence), FRC_dd) +
#   labs(x = "Resource starting density (mL-1)",
#        y = "End - start resource density",
#        color = "Consumer presence",
#        title = "Selection environment C") +
#   scale_color_manual(values = 1:2) +
#   theme_bw()
# 
# # What about predicted consumption?
# FRC_predictions <- FRC_predictions %>%
#   mutate(Resource_consumed = exp(coef(best.mod.fit.full.c)[["bC.log"]]) * Resource_start_density / (1 + exp(coef(best.mod.fit.full.c)[["bC.log"]]) * coef(best.mod.fit.full.c)[["hC"]] * Resource_start_density))
# 
# ggplot(FRC_predictions, aes(x = Resource_start_density, y = Resource_consumed)) +
#   geom_line() +
#   theme_bw()
# 
# 
# 
# 
# # ---------------------- FRT ----------------------
# # Let's start by seeing how well resource growth is predicted by the model for environment C
# FRT_dd <- filter(FR_all, SEL_environment == "T")
# 
# ggplot(FRT_dd, aes(x = Resource_start_density, y = Resource_end_density - Resource_start_density, color = Consumer_presence)) +
#   geom_point()
# 
# FRT_predictions <- data.frame(Resource_start_density = 1:max(FRT_dd$Resource_start_density))
# 
# # Predicted end density without consumers
# FRT_predictions$Resource_end_density_without <- unlist(lapply(FRT_predictions$Resource_start_density, FUN = function(N) {
#   Nt <- eq.odeint.general.pred.t2.mP(start.v = c(N, 0),
#                                   b = exp(coef(best.mod.fit.full.c)[["bC.log"]] + coef(best.mod.fit.full.c)[["bT.log"]]),
#                                   h = coef(best.mod.fit.full.c)[["hC"]] + coef(best.mod.fit.full.c)[["hT"]],
#                                   r = coef(best.mod.fit.full.c)[["rT"]],
#                                   K = exp(coef(best.mod.fit.full.c)[["KT.log"]]),
#                                   c = coef(best.mod.fit.full.c)[["c"]],
#                                   m = coef(best.mod.fit.full.c)[["mC"]] + coef(best.mod.fit.full.c)[["mT"]],
#                                   Tt = 1,
#                                   timesteplength = 1/100)[101,2]
#   return(Nt)
# }))
# 
# # Predicted end density with 5 consumers
# FRT_predictions$Resource_end_density_with <- unlist(lapply(FRT_predictions$Resource_start_density, FUN = function(N) {
#   Nt <- eq.odeint.general.pred.t2.mP(start.v = c(N, 5),
#                                      b = exp(coef(best.mod.fit.full.c)[["bC.log"]] + coef(best.mod.fit.full.c)[["bT.log"]]),
#                                      h = coef(best.mod.fit.full.c)[["hC"]] + coef(best.mod.fit.full.c)[["hT"]],
#                                      r = coef(best.mod.fit.full.c)[["rT"]],
#                                      K = exp(coef(best.mod.fit.full.c)[["KT.log"]]),
#                                      c = coef(best.mod.fit.full.c)[["c"]],
#                                      m = coef(best.mod.fit.full.c)[["mC"]] + coef(best.mod.fit.full.c)[["mT"]],
#                                   Tt = 1,
#                                   timesteplength = 1/100)[101,2]
#   return(Nt)
# }))
# 
# 
# GRT_without <- ggplot(FRT_predictions) +
#   geom_line(aes(x = Resource_start_density, y = Resource_end_density_without - Resource_start_density)) +
#   # geom_line(aes(x = Resource_start_density, y = Resource_end_density_with - Resource_start_density), color = 2) +
#   geom_point(aes(x = Resource_start_density, y = Resource_end_density - Resource_start_density), filter(FRT_dd, Consumer_presence == 0)) +
#   labs(x = "Resource starting density (mL-1)",
#        y = "End - start resource density",
#        title = "Selection environment T - Growth without consumers") +
#   theme_bw()
# 
# GRT_with <- ggplot(FRT_predictions) +
#   geom_line(aes(x = Resource_start_density, y = Resource_end_density_with - Resource_start_density)) +
#   geom_point(aes(x = Resource_start_density, y = Resource_end_density - Resource_start_density), filter(FRT_dd, Consumer_presence == 1)) +
#   labs(x = "Resource starting density (mL-1)",
#        y = "End - start resource density",
#        title = "Selection environment T - Growth with consumers") +
#   theme_bw()
# 
# ggpubr::ggarrange(GRT_without, GRT_with, ncol = 2)
# 
# 
# ggplot(FRT_predictions) +
#   geom_line(aes(x = Resource_start_density, y = Resource_end_density_without - Resource_start_density), color = 1) +
#   geom_line(aes(x = Resource_start_density, y = Resource_end_density_with - Resource_start_density), color = 2) +
#   geom_point(aes(x = Resource_start_density, y = Resource_end_density - Resource_start_density,
#                  color = Consumer_presence), FRT_dd) +
#   labs(x = "Resource starting density (mL-1)",
#        y = "End - start resource density",
#        color = "Consumer presence",
#        title = "Selection environment T") +
#   scale_color_manual(values = 1:2) +
#   theme_bw()
# 
# # What about predicted consumption?
# FRT_predictions <- FRT_predictions %>%
#   mutate(Resource_consumed = exp(coef(best.mod.fit.full.c)[["bC.log"]] + coef(best.mod.fit.full.c)[["bT.log"]]) * Resource_start_density / (1 + exp(coef(best.mod.fit.full.c)[["bC.log"]] + coef(best.mod.fit.full.c)[["bT.log"]]) * (coef(best.mod.fit.full.c)[["hC"]] + coef(best.mod.fit.full.c)[["hT"]]) * Resource_start_density))
# 
# ggplot(FRT_predictions, aes(x = Resource_start_density, y = Resource_consumed)) +
#   geom_line() +
#   geom_line(aes(x = Resource_start_density, y = Resource_consumed), FRC_predictions, color = 2) +
#   theme_bw()
# 
# 
# # ---------------------- FRS ----------------------
# # Let's start by seeing how well resource growth is predicted by the model for environment C
# FRS_dd <- filter(FR_all, SEL_environment == "S")
# 
# ggplot(FRS_dd, aes(x = Resource_start_density, y = Resource_end_density - Resource_start_density, color = Consumer_presence)) +
#   geom_point()
# 
# FRS_predictions <- data.frame(Resource_start_density = 1:max(FRS_dd$Resource_start_density))
# 
# # Predicted end density without consumers
# FRS_predictions$Resource_end_density_without <- unlist(lapply(FRS_predictions$Resource_start_density, FUN = function(N) {
#   Nt <- eq.odeint.general.pred.t2.mP(start.v = c(N, 0),
#                                      b = exp(coef(best.mod.fit.full.c)[["bC.log"]] + coef(best.mod.fit.full.c)[["bS.log"]]),
#                                      h = coef(best.mod.fit.full.c)[["hC"]] + coef(best.mod.fit.full.c)[["hS"]],
#                                      r = coef(best.mod.fit.full.c)[["rS"]],
#                                      K = exp(coef(best.mod.fit.full.c)[["KS.log"]]),
#                                      c = coef(best.mod.fit.full.c)[["c"]],
#                                      m = coef(best.mod.fit.full.c)[["mC"]] + coef(best.mod.fit.full.c)[["mS"]],
#                                   Tt = 1,
#                                   timesteplength = 1/100)[101,2]
#   return(Nt)
# }))
# 
# # Predicted end density with 5 consumers
# FRS_predictions$Resource_end_density_with <- unlist(lapply(FRS_predictions$Resource_start_density, FUN = function(N) {
#   Nt <- eq.odeint.general.pred.t2.mP(start.v = c(N, 5),
#                                      b = exp(coef(best.mod.fit.full.c)[["bC.log"]] + coef(best.mod.fit.full.c)[["bS.log"]]),
#                                      h = coef(best.mod.fit.full.c)[["hC"]] + coef(best.mod.fit.full.c)[["hS"]],
#                                      r = coef(best.mod.fit.full.c)[["rS"]],
#                                      K = exp(coef(best.mod.fit.full.c)[["KS.log"]]),
#                                      c = coef(best.mod.fit.full.c)[["c"]],
#                                      m = coef(best.mod.fit.full.c)[["mC"]] + coef(best.mod.fit.full.c)[["mS"]],
#                                   Tt = 1,
#                                   timesteplength = 1/100)[101,2]
#   return(Nt)
# }))
# 
# 
# GRS_without <- ggplot(FRS_predictions) +
#   geom_line(aes(x = Resource_start_density, y = Resource_end_density_without - Resource_start_density)) +
#   # geom_line(aes(x = Resource_start_density, y = Resource_end_density_with - Resource_start_density), color = 2) +
#   geom_point(aes(x = Resource_start_density, y = Resource_end_density - Resource_start_density), filter(FRS_dd, Consumer_presence == 0)) +
#   labs(x = "Resource starting density (mL-1)",
#        y = "End - start resource density",
#        title = "Selection environment S - Growth without consumers") +
#   theme_bw()
# 
# GRS_with <- ggplot(FRS_predictions) +
#   geom_line(aes(x = Resource_start_density, y = Resource_end_density_with - Resource_start_density)) +
#   geom_point(aes(x = Resource_start_density, y = Resource_end_density - Resource_start_density), filter(FRS_dd, Consumer_presence == 1)) +
#   labs(x = "Resource starting density (mL-1)",
#        y = "End - start resource density",
#        title = "Selection environment S - Growth with consumers") +
#   theme_bw()
# 
# ggpubr::ggarrange(GRS_without, GRS_with, ncol = 2)
# 
# 
# ggplot(FRS_predictions) +
#   geom_line(aes(x = Resource_start_density, y = Resource_end_density_without - Resource_start_density), color = 1) +
#   geom_line(aes(x = Resource_start_density, y = Resource_end_density_with - Resource_start_density), color = 2) +
#   geom_point(aes(x = Resource_start_density, y = Resource_end_density - Resource_start_density,
#                  color = Consumer_presence), FRS_dd) +
#   labs(x = "Resource starting density (mL-1)",
#        y = "End - start resource density",
#        color = "Consumer presence",
#        title = "Selection environment S") +
#   scale_color_manual(values = 1:2) +
#   theme_bw()
# 
# # What about predicted consumption?
# FRS_predictions <- FRS_predictions %>%
#   mutate(Resource_consumed = exp(coef(best.mod.fit.full.c)[["bC.log"]] + coef(best.mod.fit.full.c)[["bS.log"]]) * Resource_start_density / (1 + exp(coef(best.mod.fit.full.c)[["bC.log"]] + coef(best.mod.fit.full.c)[["bS.log"]]) * (coef(best.mod.fit.full.c)[["hC"]] + coef(best.mod.fit.full.c)[["hS"]]) * Resource_start_density))
# 
# ggplot(FRS_predictions, aes(x = Resource_start_density, y = Resource_consumed)) +
#   geom_line() +
#   geom_line(aes(x = Resource_start_density, y = Resource_consumed), FRC_predictions, color = 2) +
#   geom_line(aes(x = Resource_start_density, y = Resource_consumed), FRT_predictions, color = 3) +
#   theme_bw()
# 
# 
# # ---------------------- FRST ----------------------
# # Let's start by seeing how well resource growth is predicted by the model for environment C
# FRST_dd <- filter(FR_all, SEL_environment == "ST")
# 
# ggplot(FRST_dd, aes(x = Resource_start_density, y = Resource_end_density - Resource_start_density, color = Consumer_presence)) +
#   geom_point()
# 
# FRST_predictions <- data.frame(Resource_start_density = 1:max(FRST_dd$Resource_start_density))
# 
# # Predicted end density without consumers
# FRST_predictions$Resource_end_density_without <- unlist(lapply(FRST_predictions$Resource_start_density, FUN = function(N) {
#   Nt <- eq.odeint.general.pred.t2.mP(start.v = c(N, 0),
#                                      b = exp(coef(best.mod.fit.full.c)[["bC.log"]] + coef(best.mod.fit.full.c)[["bS.log"]] + coef(best.mod.fit.full.c)[["bT.log"]] + coef(best.mod.fit.full.c)[["bST.log"]]),
#                                      h = coef(best.mod.fit.full.c)[["hC"]] + coef(best.mod.fit.full.c)[["hS"]] + coef(best.mod.fit.full.c)[["hT"]] + coef(best.mod.fit.full.c)[["hST"]],
#                                      r = coef(best.mod.fit.full.c)[["rST"]],
#                                      K = exp(coef(best.mod.fit.full.c)[["KST.log"]]),
#                                      c = coef(best.mod.fit.full.c)[["c"]],
#                                      m = coef(best.mod.fit.full.c)[["mC"]] + coef(best.mod.fit.full.c)[["mS"]] + coef(best.mod.fit.full.c)[["mT"]] + coef(best.mod.fit.full.c)[["mST"]],
#                                   Tt = 1,
#                                   timesteplength = 1/100)[101,2]
#   return(Nt)
# }))
# 
# # Predicted end density with 5 consumers
# FRST_predictions$Resource_end_density_with <- unlist(lapply(FRST_predictions$Resource_start_density, FUN = function(N) {
#   Nt <- eq.odeint.general.pred.t2.mP(start.v = c(N, 5),
#                                      b = exp(coef(best.mod.fit.full.c)[["bC.log"]] + coef(best.mod.fit.full.c)[["bS.log"]] + coef(best.mod.fit.full.c)[["bT.log"]] + coef(best.mod.fit.full.c)[["bST.log"]]),
#                                      h = coef(best.mod.fit.full.c)[["hC"]] + coef(best.mod.fit.full.c)[["hS"]] + coef(best.mod.fit.full.c)[["hT"]] + coef(best.mod.fit.full.c)[["hST"]],
#                                      r = coef(best.mod.fit.full.c)[["rST"]],
#                                      K = exp(coef(best.mod.fit.full.c)[["KST.log"]]),
#                                      c = coef(best.mod.fit.full.c)[["c"]],
#                                      m = coef(best.mod.fit.full.c)[["mC"]] + coef(best.mod.fit.full.c)[["mS"]] + coef(best.mod.fit.full.c)[["mT"]] + coef(best.mod.fit.full.c)[["mST"]],
#                                   Tt = 1,
#                                   timesteplength = 1/100)[101,2]
#   return(Nt)
# }))
# 
# 
# GRST_without <- ggplot(FRST_predictions) +
#   geom_line(aes(x = Resource_start_density, y = Resource_end_density_without - Resource_start_density)) +
#   # geom_line(aes(x = Resource_start_density, y = Resource_end_density_with - Resource_start_density), color = 2) +
#   geom_point(aes(x = Resource_start_density, y = Resource_end_density - Resource_start_density), filter(FRST_dd, Consumer_presence == 0)) +
#   labs(x = "Resource starting density (mL-1)",
#        y = "End - start resource density",
#        title = "Selection environment ST - Growth without consumers") +
#   theme_bw()
# 
# GRST_with <- ggplot(FRST_predictions) +
#   geom_line(aes(x = Resource_start_density, y = Resource_end_density_with - Resource_start_density)) +
#   geom_point(aes(x = Resource_start_density, y = Resource_end_density - Resource_start_density), filter(FRST_dd, Consumer_presence == 1)) +
#   labs(x = "Resource starting density (mL-1)",
#        y = "End - start resource density",
#        title = "Selection environment ST - Growth with consumers") +
#   theme_bw()
# 
# ggpubr::ggarrange(GRST_without, GRST_with, ncol = 2)
# 
# 
# ggplot(FRST_predictions) +
#   geom_line(aes(x = Resource_start_density, y = Resource_end_density_without - Resource_start_density), color = 1) +
#   geom_line(aes(x = Resource_start_density, y = Resource_end_density_with - Resource_start_density), color = 2) +
#   geom_point(aes(x = Resource_start_density, y = Resource_end_density - Resource_start_density,
#                  color = Consumer_presence), FRST_dd) +
#   labs(x = "Resource starting density (mL-1)",
#        y = "End - start resource density",
#        color = "Consumer presence",
#        title = "Selection environment ST") +
#   scale_color_manual(values = 1:2) +
#   theme_bw()
# 
# # What about predicted consumption?
# FRST_predictions <- FRST_predictions %>%
#   mutate(Resource_consumed = exp(coef(best.mod.fit.full.c)[["bC.log"]] + coef(best.mod.fit.full.c)[["bS.log"]] + coef(best.mod.fit.full.c)[["bT.log"]] + coef(best.mod.fit.full.c)[["bST.log"]]) * Resource_start_density / (1 + exp(coef(best.mod.fit.full.c)[["bC.log"]] + coef(best.mod.fit.full.c)[["bS.log"]] + coef(best.mod.fit.full.c)[["bT.log"]] + coef(best.mod.fit.full.c)[["bST.log"]]) * (coef(best.mod.fit.full.c)[["hC"]] + coef(best.mod.fit.full.c)[["hS"]] + coef(best.mod.fit.full.c)[["hT"]] + coef(best.mod.fit.full.c)[["hST"]]) * Resource_start_density))
# 
# ggplot(FRST_predictions, aes(x = Resource_start_density, y = Resource_consumed)) +
#   geom_line() +
#   geom_line(aes(x = Resource_start_density, y = Resource_consumed), FRC_predictions, color = 2) +
#   geom_line(aes(x = Resource_start_density, y = Resource_consumed), FRT_predictions, color = 3) +
#   geom_line(aes(x = Resource_start_density, y = Resource_consumed), FRS_predictions, color = 4) +
#   theme_bw()
# 
# 
# 
# # ----------------- All together -----------------
# FRC_predictions$SEL_environment = "C"
# FRS_predictions$SEL_environment = "S"
# FRT_predictions$SEL_environment = "T"
# FRST_predictions$SEL_environment = "ST"
# 
# FR_predictions <- rbind(FRC_predictions,
#                         FRS_predictions,
#                         FRT_predictions,
#                         FRST_predictions)
# 
# ggplot(FR_predictions, aes(x = Resource_start_density, y = Resource_consumed, color = SEL_environment)) +
#   geom_line() +
#   theme_bw()
#   
# FR_predictions_c <- FR_predictions
# 
# # save(FR_predictions_c, file = "FR data/Model predictions with simulated parameters/ANOVA-like fits/FR_predicted_fits_ANOVA_c.RData")






###### Simulations with varying parameters for 95% CIs ####

# ---------------------- CIs ----------------------
set.seed(42)

# Simulate 4000 realizations of FR parameter combinations
n_sim <- 1000
param_sim <- mvrnorm(n_sim, mu = FR_coefs$estimate, Sigma = vcov(best.mod.fit.full.c)) %>%
  data.frame()

# Adding the r and K information for each SEL environment
param_sim$rC <- coef(best.mod.fit.full.c)[["rC"]]
param_sim$rS <- coef(best.mod.fit.full.c)[["rS"]]
param_sim$rT <- coef(best.mod.fit.full.c)[["rT"]]
param_sim$rST <- coef(best.mod.fit.full.c)[["rST"]]
param_sim$KC.log <- coef(best.mod.fit.full.c)[["KC.log"]]
param_sim$KS.log <- coef(best.mod.fit.full.c)[["KS.log"]]
param_sim$KT.log <- coef(best.mod.fit.full.c)[["KT.log"]]
param_sim$KST.log <- coef(best.mod.fit.full.c)[["KST.log"]]

param_sim <- param_sim %>%
  as.matrix()


N_C <- 1:max(FR_predictions$Resource_start_density[FR_predictions$SEL_environment == "C"])
N_S <- 1:max(FR_predictions$Resource_start_density[FR_predictions$SEL_environment == "S"])
N_T <- 1:max(FR_predictions$Resource_start_density[FR_predictions$SEL_environment == "T"])
N_ST <- 1:max(FR_predictions$Resource_start_density[FR_predictions$SEL_environment == "ST"])


param_sim_c <- param_sim
save(param_sim_c, N_C, N_S, N_T, N_ST,
     file = "FR data/Model predictions with simulated parameters/ANOVA-like fits/1000 simulations/FR_ANOVA_simulation_parameters_c.RData")

load("FR data/Model predictions with simulated parameters/ANOVA-like fits/1000 simulations/FR_ANOVA_simulation_parameters_c.RData")

# Compute the number of prey eaten for all resource densities and varying parameters
### Output will be a matrix where each column is the estimated consumption of
### one of the parameter combinations
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

predictions_dd_C <- data.frame(predictions_C) %>%
  pivot_longer(cols = 1:n_sim, names_to = "Simulation_nr", values_to = "Ndead")
predictions_dd_C$N <- rep(N_C, each = n_sim)
predictions_dd_C$b.log.estimate <- rep(param_sim[,"bC.log"], times = length(N_C))
predictions_dd_C$h.estimate <- rep(param_sim[,"hC"], times = length(N_C))
predictions_dd_C$m.estimate <- rep(param_sim[,"mC"], times = length(N_C))

predictions_dd_S <- data.frame(predictions_S) %>%
  pivot_longer(cols = 1:n_sim, names_to = "Simulation_nr", values_to = "Ndead")
predictions_dd_S$N <- rep(N_S, each = n_sim)
predictions_dd_S$b.log.estimate <- rep(param_sim[,"bC.log"] + param_sim[,"bS.log"], times = length(N_S))
predictions_dd_S$h.estimate <- rep(param_sim[,"hC"] + param_sim[,"hS"], times = length(N_S))
predictions_dd_S$m.estimate <- rep(param_sim[,"mC"] + param_sim[,"mS"], times = length(N_S))

predictions_dd_T <- data.frame(predictions_T) %>%
  pivot_longer(cols = 1:n_sim, names_to = "Simulation_nr", values_to = "Ndead")
predictions_dd_T$N <- rep(N_T, each = n_sim)
predictions_dd_T$b.log.estimate <- rep(param_sim[,"bC.log"] + param_sim[,"bT.log"], times = length(N_T))
predictions_dd_T$h.estimate <- rep(param_sim[,"hC"] + param_sim[,"hT"], times = length(N_T))
predictions_dd_T$m.estimate <- rep(param_sim[,"mC"] + param_sim[,"mT"], times = length(N_T))

predictions_dd_ST <- data.frame(predictions_ST) %>%
  pivot_longer(cols = 1:n_sim, names_to = "Simulation_nr", values_to = "Ndead")
predictions_dd_ST$N <- rep(N_ST, each = n_sim)
predictions_dd_ST$b.log.estimate <- rep(param_sim[,"bC.log"] + param_sim[,"bS.log"] + param_sim[,"bT.log"] + param_sim[,"bST.log"], times = length(N_ST))
predictions_dd_ST$h.estimate <- rep(param_sim[,"hC"] + param_sim[,"hS"] + param_sim[,"hT"] + param_sim[,"hST"], times = length(N_ST))
predictions_dd_ST$m.estimate <- rep(param_sim[,"mC"] + param_sim[,"mS"] + param_sim[,"mT"] + param_sim[,"mST"], times = length(N_ST))





############## Checking the predictions and filtering
# ------------------ C ------------------
ggplot(predictions_dd_C, aes(x = N, y = Ndead, group = Simulation_nr)) +
  geom_line(linewidth = .5, alpha = .3)

# There are no predictions that look weird

# All columns that have biologically impossible predictions will be NAs
predictions_filtered_C <- apply(predictions_C, 2, function(col) {
  if (sum(col < 0) > 0) {
    return(rep(NA, length(col)))
  } else if (sum(col > 7500) > 0) { # Value used to exclude only really weird fits
    return(rep(NA, length(col)))
  } else {
    return(col)
  }
})

predictions_filtered_dd_C <- data.frame(predictions_filtered_C) %>%
  pivot_longer(cols = 1:n_sim, names_to = "Simulation_nr", values_to = "Ndead")
predictions_filtered_dd_C$N <- rep(N_C, each = n_sim)

predictions_filtered_dd_C$b.log.estimate <- rep(param_sim[,"bC.log"], times = length(N_C))
predictions_filtered_dd_C$h.estimate <- rep(param_sim[,"hC"], times = length(N_C))
predictions_filtered_dd_C$m.estimate <- rep(param_sim[,"mC"], times = length(N_C))

ggplot(predictions_filtered_dd_C, aes(x = N, y = Ndead, group = Simulation_nr)) +
  geom_line(linewidth = .5, alpha = .3)


# Compute mean and confidence intervals
response_mean_C <- apply(predictions_filtered_C, 1, mean, na.rm = T)
response_lower_C <- apply(predictions_filtered_C, 1, quantile, probs = 0.025, na.rm = T)
response_upper_C <- apply(predictions_filtered_C, 1, quantile, probs = 0.975, na.rm = T)

# Visualize
plot(N_C, response_mean_C, type = "l", lwd = 2, col = "blue",
     xlab = "Resource Density (N)", ylab = "Functional Response",
     ylim = range(c(response_lower_C, response_upper_C)),
     main = "Functional Response with Confidence Intervals - SEL S")
polygon(c(N_C, rev(N_C)),
        c(response_lower_C, rev(response_upper_C)),
        col = rgb(0, 0, 1, 0.2), border = NA)


# Manually checking the biologically weird predictions
weird_fits_IDs_C <- c()
for(i in 1:dim(predictions_filtered_C)[2]) {
  
  if(sum(is.na(predictions_filtered_C[,i])) == length(predictions_filtered_C[,i])) {
    weird_fits_IDs_C <- c(weird_fits_IDs_C, i)
  }
  
} # NULL because no fits are weird

predictions_C_c <- predictions_C
predictions_filtered_C_c <- predictions_filtered_C
response_mean_C_c <- response_mean_C
response_lower_C_c <- response_lower_C
response_upper_C_c <- response_upper_C

save(predictions_C_c,
     predictions_filtered_C_c, #predictions_filtered_dd_C,
     response_mean_C_c, response_lower_C_c, response_upper_C_c,
     file = "FR data/Model predictions with simulated parameters/ANOVA-like fits/1000 simulations/FR_ANOVA_multiple_model_predictions_C_c.RData")

load("FR data/Model predictions with simulated parameters/ANOVA-like fits/1000 simulations/FR_ANOVA_multiple_model_predictions_C_c.RData")











# Compute the end resource density for all starting resource densities and varying parameters
### Output will be a matrix where each column is the estimated end resource density of
### one of the parameter combinations
# time1 <- Sys.time()
# predictionLV_C <- apply(param_sim, 1, function(p) {
#   Nend <- unlist(lapply(N_C, FUN = function(n) {
#     eq.odeint.general.pred.t2.mP(start.v = c(n, 5),
#                                  b = exp(p[["bC.log"]]),
#                                  h = p[["hC"]],
#                                  r = p[["rC"]],
#                                  K = exp(p[["KC.log"]]),
#                                  c = p[["c"]],
#                                  m = p[["mC"]],
#                                  Tt = 1,
#                                  timesteplength = 1/100)[100+1,2]
#   }))
# })
# Sys.time() - time1

# In parallel
library(doParallel)
bigstatsr::nb_cores()

my.cluster <- parallel::makeCluster(
  bigstatsr::nb_cores(), 
  type = "PSOCK"
)

doParallel::registerDoParallel(cl = my.cluster)
foreach::getDoParRegistered()
foreach::getDoParWorkers()

predictionLV_C <- NA # name of file where parallel runs will be saved

(t0 <- Sys.time()) # timer
mcoptions <- list(preschedule=FALSE)

predictionLV_C <- foreach(i = 1:nrow(param_sim),
                          .combine = 'cbind',
                          .options.multicore = mcoptions,
                          .packages = c('bbmle', 'odeintr', 'tidyverse', 'R.utils'),
                          .errorhandling = "pass") %dopar%
  {
    compile_sys("FRII_gen_pred_mP", FRII.gen.pred.mP, pars = c("b","h","r","K","c","m"), method = "rk54")
    
    p <- param_sim[i,]
    
    predictionLV_column <- unlist(lapply(N_C, FUN = function(n) {
      eq.odeint.general.pred.t2.mP(start.v = c(n, 5),
                                   b = exp(p[["bC.log"]]),
                                   h = p[["hC"]],
                                   r = p[["rC"]],
                                   K = exp(p[["KC.log"]]),
                                   c = p[["c"]],
                                   m = p[["mC"]],
                                   Tt = 1,
                                   timesteplength = 1/100)[100+1,2]
    }))
    
    return(predictionLV_column)
  }

parallel::stopCluster(cl = my.cluster)
stopImplicitCluster()

tend<-Sys.time()
tend-t0 
# beep(sound = 2, expr = NULL)






predictionLV_dd_C <- data.frame(predictionLV_C) %>%
  pivot_longer(cols = 1:n_sim, names_to = "Simulation_nr", values_to = "Nend")
predictionLV_dd_C$N <- rep(N_C, each = n_sim)
predictionLV_dd_C$b.log.estimate <- rep(param_sim[,"bC.log"], times = length(N_C))
predictionLV_dd_C$h.estimate <- rep(param_sim[,"hC"], times = length(N_C))
predictionLV_dd_C$m.estimate <- rep(param_sim[,"mC"], times = length(N_C))


predictionLV_dd_C$Simulation_nr <- gsub(predictionLV_dd_C$Simulation_nr, pattern = "result.", replacement = "X")

ggplot() +
  geom_line(aes(x = N, y = Nend, group = Simulation_nr, color = b.log.estimate), filter(predictionLV_dd_C, !(Simulation_nr %in% paste0("X", weird_fits_IDs_C))), size = .5, alpha = .8) +
  geom_point(aes(x = Resource_start_density, y = Resource_end_density), data = filter(FR_all, SEL_environment == "C", Consumer_presence == 1))


# No weird fits = no filtering
predictionLV_filtered_C <- predictionLV_C


# Compute mean and confidence intervals
LV_mean_C <- apply(predictionLV_filtered_C, 1, mean, na.rm = T)
LV_lower_C <- apply(predictionLV_filtered_C, 1, quantile, probs = 0.025, na.rm = T)
LV_upper_C <- apply(predictionLV_filtered_C, 1, quantile, probs = 0.975, na.rm = T)

predictionLV_C_c <- predictionLV_C
predictionLV_filtered_C_c <- predictionLV_filtered_C
LV_mean_C_c <- LV_mean_C
LV_lower_C_c <- LV_lower_C
LV_upper_C_c <- LV_upper_C

save(predictionLV_C_c,
     predictionLV_filtered_C_c,
     LV_mean_C_c, LV_lower_C_c, LV_upper_C_c,
     file = "FR data/Model predictions with simulated parameters/ANOVA-like fits/1000 simulations/FR_ANOVA_multiple_LV_predictions_C_c.RData")
# load("FR data/Model predictions with simulated parameters/ANOVA-like fits/1000 simulations/FR_ANOVA_multiple_LV_predictions_C_c.RData")

predictionLV_dd_C_c <- predictionLV_dd_C
predictions_dd_C_c <- predictions_dd_C

save(predictionLV_dd_C_c, predictions_dd_C_c,
     file = "FR data/Model predictions with simulated parameters/ANOVA-like fits/1000 simulations/FR_ANOVA_sim_plot_data_C_c.RData")
load("FR data/Model predictions with simulated parameters/ANOVA-like fits/1000 simulations/FR_ANOVA_sim_plot_data_C_c.RData")




# ------------------ S ------------------
ggplot(predictions_dd_S, aes(x = N, y = Ndead, group = Simulation_nr)) +
  geom_line(linewidth = .5, alpha = .3)

# There are no predictions that look weird

# All columns that have biologically impossible predictions will be NAs
predictions_filtered_S <- apply(predictions_S, 2, function(col) {
  if (sum(col < 0) > 0) {
    return(rep(NA, length(col)))
  } else if (sum(col > 12000) > 0) { # Value used to exclude only really weird fits
    return(rep(NA, length(col)))
  } else {
    return(col)
  }
})

predictions_filtered_dd_S <- data.frame(predictions_filtered_S) %>%
  pivot_longer(cols = 1:n_sim, names_to = "Simulation_nr", values_to = "Ndead")
predictions_filtered_dd_S$N <- rep(N_S, each = n_sim)

predictions_filtered_dd_S$b.log.estimate <- rep(param_sim[,"bC.log"] + param_sim[,"bS.log"] + param_sim[,"bS.log"], times = length(N_S))
predictions_filtered_dd_S$h.estimate <- rep(param_sim[,"hC"] + param_sim[,"hS"], times = length(N_S))
predictions_filtered_dd_S$m.estimate <- rep(param_sim[,"mC"] + param_sim[,"mS"], times = length(N_S))

ggplot(predictions_filtered_dd_S, aes(x = N, y = Ndead, group = Simulation_nr)) +
  geom_line(linewidth = .5, alpha = .3)


# Compute mean and confidence intervals
response_mean_S <- apply(predictions_filtered_S, 1, mean, na.rm = T)
response_lower_S <- apply(predictions_filtered_S, 1, quantile, probs = 0.025, na.rm = T)
response_upper_S <- apply(predictions_filtered_S, 1, quantile, probs = 0.975, na.rm = T)

# Visualize
plot(N_S, response_mean_S, type = "l", lwd = 2, col = "blue",
     xlab = "Resource Density (N)", ylab = "Functional Response",
     ylim = range(c(response_lower_S, response_upper_S)),
     main = "Functional Response with Confidence Intervals - SEL S")
polygon(c(N_S, rev(N_S)),
        c(response_lower_S, rev(response_upper_S)),
        col = rgb(0, 0, 1, 0.2), border = NA)


# Manually checking the biologically weird predictions
weird_fits_IDs_S <- c()
for(i in 1:dim(predictions_filtered_S)[2]) {
  
  if(sum(is.na(predictions_filtered_S[,i])) == length(predictions_filtered_S[,i])) {
    weird_fits_IDs_S <- c(weird_fits_IDs_S, i)
  }
  
} # NULL


predictions_S_c <- predictions_S
predictions_filtered_S_c <- predictions_filtered_S
predictions_filtered_dd_S_c <- predictions_filtered_dd_S
response_mean_S_c <- response_mean_S
response_lower_S_c <- response_lower_S
response_upper_S_c <- response_upper_S

save(predictions_S_c, #predictions_dd_S_c,
     predictions_filtered_S_c, predictions_filtered_dd_S_c,
     response_mean_S_c, response_lower_S_c, response_upper_S_c,
     file = "FR data/Model predictions with simulated parameters/ANOVA-like fits/1000 simulations/FR_ANOVA_multiple_model_predictions_S_c.RData")
# load("FR data/Model predictions with simulated parameters/ANOVA-like fits/1000 simulations/FR_ANOVA_multiple_model_predictions_S_c.RData")











# Compute the end resource density for all starting resource densities and varying parameters
### Output will be a matrix where each column is the estimated end resource density of
### one of the parameter combinations
# time1 <- Sys.time()
# predictionLV_S <- apply(param_sim, 1, function(p) {
#   Nend <- unlist(lapply(N_S, FUN = function(n) {
#     eq.odeint.general.pred.t2.mP(start.v = c(n, 5),
#                                  b = exp(p[["bC.log"]]),
#                                  h = p[["hC"]],
#                                  r = p[["rC"]],
#                                  K = exp(p[["KC.log"]]),
#                                  c = p[["c"]],
#                                  m = p[["mC"]],
#                                  Tt = 1,
#                                  timesteplength = 1/100)[100+1,2]
#   }))
# })
# Sys.time() - time1

# In parallel
library(doParallel)
bigstatsr::nb_cores()

my.cluster <- parallel::makeCluster(
  bigstatsr::nb_cores(), 
  type = "PSOCK"
)

doParallel::registerDoParallel(cl = my.cluster)
foreach::getDoParRegistered()
foreach::getDoParWorkers()

predictionLV_S <- NA # name of file where parallel runs will be save

(t0 <- Sys.time()) # timer
mcoptions <- list(preschedule=FALSE)

predictionLV_S <- foreach(i = 1:nrow(param_sim),
                          .combine = 'cbind',
                          .options.multicore = mcoptions,
                          .packages = c('bbmle', 'odeintr', 'tidyverse', 'R.utils'),
                          .errorhandling = "pass") %dopar%
  {
    compile_sys("FRII_gen_pred_mP", FRII.gen.pred.mP, pars = c("b","h","r","K","c","m"), method = "rk54")
    
    p <- param_sim[i,]
    
    predictionLV_column <- unlist(lapply(N_S, FUN = function(n) {
      eq.odeint.general.pred.t2.mP(start.v = c(n, 5),
                                   b = exp(p[["bC.log"]] + p[["bS.log"]]),
                                   h = p[["hC"]] + p[["hS"]],
                                   r = p[["rS"]],
                                   K = exp(p[["KS.log"]]),
                                   c = p[["c"]],
                                   m = p[["mC"]] + p[["mS"]],
                                   Tt = 1,
                                   timesteplength = 1/100)[100+1,2]
    }))
    
    return(predictionLV_column)
  }

parallel::stopCluster(cl = my.cluster)
stopImplicitCluster()

tend<-Sys.time()
tend-t0 
# beep(sound = 2, expr = NULL)






predictionLV_dd_S <- data.frame(predictionLV_S) %>%
  pivot_longer(cols = 1:n_sim, names_to = "Simulation_nr", values_to = "Nend")
predictionLV_dd_S$N <- rep(N_S, each = n_sim)
predictionLV_dd_S$b.log.estimate <- rep(param_sim[,"bC.log"] + param_sim[,"bS.log"], times = length(N_S))
predictionLV_dd_S$h.estimate <- rep(param_sim[,"hC"] + param_sim[,"hS"], times = length(N_S))
predictionLV_dd_S$m.estimate <- rep(param_sim[,"mC"] + param_sim[,"mS"], times = length(N_S))


predictionLV_dd_S$Simulation_nr <- gsub(predictionLV_dd_S$Simulation_nr, pattern = "result.", replacement = "X")

ggplot() +
  geom_line(aes(x = N, y = Nend, group = Simulation_nr, color = b.log.estimate), filter(predictionLV_dd_S, !(Simulation_nr %in% paste0("X", weird_fits_IDs_S))), size = .5, alpha = .8) +
  geom_point(aes(x = Resource_start_density, y = Resource_end_density), data = filter(FR_all, SEL_environment == "C", Consumer_presence == 1))


# There are two lines that don't fit the data well, but they're continuous so no exclusion
ggplot() +
  geom_line(aes(x = N, y = Nend, group = Simulation_nr, color = b.log.estimate), filter(predictionLV_dd_S, Simulation_nr %in% paste0("X", weird_fits_IDs_S)), size = .5, alpha = .8) +
  geom_point(aes(x = Resource_start_density, y = Resource_end_density), data = filter(FR_all, SEL_environment == "C", Consumer_presence == 1))


# Exclude weird fits from this dataset too
predictionLV_filtered_S <- predictionLV_S


# Compute mean and confidence intervals
LV_mean_S <- apply(predictionLV_filtered_S, 1, mean, na.rm = T)
LV_lower_S <- apply(predictionLV_filtered_S, 1, quantile, probs = 0.025, na.rm = T)
LV_upper_S <- apply(predictionLV_filtered_S, 1, quantile, probs = 0.975, na.rm = T)


predictionLV_S_c <- predictionLV_S
predictionLV_filtered_S_c <- predictionLV_filtered_S
predictionLV_dd_S_c <- predictionLV_dd_S
LV_mean_S_c <- LV_mean_S
LV_lower_S_c <- LV_lower_S
LV_upper_S_c <- LV_upper_S

save(predictionLV_S_c, predictionLV_filtered_S_c, predictionLV_dd_S_c,
     LV_mean_S_c, LV_lower_S_c, LV_upper_S_c,
     file = "FR data/Model predictions with simulated parameters/ANOVA-like fits/1000 simulations/FR_ANOVA_multiple_LV_predictions_S_c.RData")
# load("FR data/Model predictions with simulated parameters/ANOVA-like fits/1000 simulations/FR_ANOVA_multiple_LV_predictions_S_c.RData")

predictionLV_dd_S_c <- predictionLV_dd_S
predictions_dd_S_c <- predictions_dd_S
predictions_dd_S_c <- predictions_dd_S

save(predictionLV_dd_S_c, predictions_dd_S_c, predictions_dd_S_c,
     file = "FR data/Model predictions with simulated parameters/ANOVA-like fits/1000 simulations/FR_ANOVA_sim_plot_data_S_c.RData")
# load("FR data/Model predictions with simulated parameters/ANOVA-like fits/1000 simulations/FR_ANOVA_sim_plot_data_S_c.RData")

# ------------------ T ------------------
ggplot(predictions_dd_T, aes(x = N, y = Ndead, group = Simulation_nr)) +
  geom_line(linewidth = .5, alpha = .3)

# There are no weird fits

# All columns that have biologically impossible predictions will be NAs
predictions_filtered_T <- apply(predictions_T, 2, function(col) {
  if (sum(col < 0) > 0) {
    return(rep(NA, length(col)))
  } else if (sum(col > 2000) > 0) { # Value used to exclude only really weird fits
    return(rep(NA, length(col)))
  } else {
    return(col)
  }
})

predictions_filtered_dd_T <- data.frame(predictions_filtered_T) %>%
  pivot_longer(cols = 1:n_sim, names_to = "Simulation_nr", values_to = "Ndead")
predictions_filtered_dd_T$N <- rep(N_T, each = n_sim)

predictions_filtered_dd_T$b.log.estimate <- rep(param_sim[,"bC.log"] + param_sim[,"bT.log"], times = length(N_T))
predictions_filtered_dd_T$h.estimate <- rep(param_sim[,"hC"] + param_sim[,"hT"], times = length(N_T))
predictions_filtered_dd_T$m.estimate <- rep(param_sim[,"mC"]+ param_sim[,"mT"], times = length(N_T))

# Compute mean and confidence intervals
response_mean_T <- apply(predictions_filtered_T, 1, mean, na.rm = T)
response_lower_T <- apply(predictions_filtered_T, 1, quantile, probs = 0.025, na.rm = T)
response_upper_T <- apply(predictions_filtered_T, 1, quantile, probs = 0.975, na.rm = T)

# Visualize
plot(N_T, response_mean_T, type = "l", lwd = 2, col = "blue",
     xlab = "Resource Density (N)", ylab = "Functional Response",
     ylim = range(c(response_lower_T, response_upper_T)),
     main = "Functional Response with Confidence Intervals - SEL S")
polygon(c(N_T, rev(N_T)),
        c(response_lower_T, rev(response_upper_T)),
        col = rgb(0, 0, 1, 0.2), border = NA)


# Manually checking the biologically weird predictions
weird_fits_IDs_T <- c()
for(i in 1:dim(predictions_filtered_T)[2]) {
  
  if(sum(is.na(predictions_filtered_T[,i])) == length(predictions_filtered_T[,i])) {
    weird_fits_IDs_T <- c(weird_fits_IDs_T, i)
  }
  
} # NULL object - no weird fits

predictions_T_c <- predictions_T
predictions_dd_T_c <- predictions_dd_T
predictions_filtered_T_c <- predictions_filtered_T
predictions_filtered_dd_T_c <- predictions_filtered_dd_T
response_mean_T_c <- response_mean_T
response_lower_T_c <- response_lower_T
response_upper_T_c <- response_upper_T

save(predictions_T_c, predictions_dd_T_c,
     predictions_filtered_T_c, predictions_filtered_dd_T_c,
     response_mean_T_c, response_lower_T_c, response_upper_T_c,
     file = "FR data/Model predictions with simulated parameters/ANOVA-like fits/1000 simulations/FR_ANOVA_multiple_model_predictions_T_c.RData")
# load("FR data/Model predictions with simulated parameters/ANOVA-like fits/1000 simulations/FR_ANOVA_multiple_model_predictions_T_c.RData")












# Compute the end resource density for all starting resource densities and varying parameters
### Output will be a matrix where each column is the estimated end resource density of
### one of the parameter combinations
# time1 <- Sys.time()
# predictionLV_T <- apply(param_sim, 1, function(p) {
#   Nend <- unlist(lapply(N_T, FUN = function(n) {
#     eq.odeint.general.pred.t2.mP(start.v = c(n, 5),
#                                  b = exp(p[["bC.log"]]),
#                                  h = p[["hC"]],
#                                  r = p[["rC"]],
#                                  K = exp(p[["KC.log"]]),
#                                  c = p[["c"]],
#                                  m = p[["mC"]],
#                                  Tt = 1,
#                                  timesteplength = 1/100)[100+1,2]
#   }))
# })
# Sys.time() - time1

# In parallel
library(doParallel)
bigstatsr::nb_cores()

my.cluster <- parallel::makeCluster(
  bigstatsr::nb_cores(), 
  type = "PSOCK"
)

doParallel::registerDoParallel(cl = my.cluster)
foreach::getDoParRegistered()
foreach::getDoParWorkers()

predictionLV_T <- NA # name of file where parallel runs will be save

(t0 <- Sys.time()) # timer
mcoptions <- list(preschedule=FALSE)

predictionLV_T <- foreach(i = 1:nrow(param_sim),
                          .combine = 'cbind',
                          .options.multicore = mcoptions,
                          .packages = c('bbmle', 'odeintr', 'tidyverse', 'R.utils'),
                          .errorhandling = "pass") %dopar%
  {
    compile_sys("FRII_gen_pred_mP", FRII.gen.pred.mP, pars = c("b","h","r","K","c","m"), method = "rk54")
    
    p <- param_sim[i,]
    
    predictionLV_column <- unlist(lapply(N_T, FUN = function(n) {
      eq.odeint.general.pred.t2.mP(start.v = c(n, 5),
                                   b = exp(p[["bC.log"]] + p[["bT.log"]]),
                                   h = p[["hC"]] + p[["hT"]],
                                   r = p[["rT"]],
                                   K = exp(p[["KT.log"]]),
                                   c = p[["c"]],
                                   m = p[["mC"]] + p[["mT"]],
                                   Tt = 1,
                                   timesteplength = 1/100)[100+1,2]
    }))
    
    return(predictionLV_column)
  }

parallel::stopCluster(cl = my.cluster)
stopImplicitCluster()

tend<-Sys.time()
tend-t0 
# beep(sound = 2, expr = NULL)






predictionLV_dd_T <- data.frame(predictionLV_T) %>%
  pivot_longer(cols = 1:n_sim, names_to = "Simulation_nr", values_to = "Nend")
predictionLV_dd_T$N <- rep(N_T, each = n_sim)
predictionLV_dd_T$b.log.estimate <- rep(param_sim[,"bC.log"] + param_sim[,"bT.log"], times = length(N_T))
predictionLV_dd_T$h.estimate <- rep(param_sim[,"hC"] + param_sim[,"hT"], times = length(N_T))
predictionLV_dd_T$m.estimate <- rep(param_sim[,"mC"] + param_sim[,"mT"], times = length(N_T))


predictionLV_dd_T$Simulation_nr <- gsub(predictionLV_dd_T$Simulation_nr, pattern = "result.", replacement = "X")

ggplot() +
  geom_line(aes(x = N, y = Nend, group = Simulation_nr, color = b.log.estimate), filter(predictionLV_dd_T, !(Simulation_nr %in% paste0("X", weird_fits_IDs_T))), size = .5, alpha = .8) +
  geom_point(aes(x = Resource_start_density, y = Resource_end_density), data = filter(FR_all, SEL_environment == "T", Consumer_presence == 1))


# Exclude weird fits from this dataset too
predictionLV_filtered_T <- predictionLV_T

# Compute mean and confidence intervals
LV_mean_T <- apply(predictionLV_filtered_T, 1, mean, na.rm = T)
LV_lower_T <- apply(predictionLV_filtered_T, 1, quantile, probs = 0.025, na.rm = T)
LV_upper_T <- apply(predictionLV_filtered_T, 1, quantile, probs = 0.975, na.rm = T)

predictionLV_T_c <- predictionLV_T
predictionLV_dd_T_c <- predictionLV_dd_T
predictionLV_filtered_T_c <- predictionLV_filtered_T
LV_mean_T_c <- LV_mean_T
LV_lower_T_c <- LV_lower_T
LV_upper_T_c <- LV_upper_T

save(predictionLV_T_c, predictionLV_filtered_T_c, predictionLV_dd_T_c,
     LV_mean_T_c, LV_lower_T_c, LV_upper_T_c,
     file = "FR data/Model predictions with simulated parameters/ANOVA-like fits/1000 simulations/FR_ANOVA_multiple_LV_predictions_T_c.RData")
# load("FR data/Model predictions with simulated parameters/ANOVA-like fits/1000 simulations/FR_ANOVA_multiple_LV_predictions_T_c.RData")

predictionLV_dd_T_c <- predictionLV_dd_T
predictions_dd_T_c <- predictions_dd_T

save(predictionLV_dd_T_c, predictions_dd_T_c,
     file = "FR data/Model predictions with simulated parameters/ANOVA-like fits/1000 simulations/FR_ANOVA_sim_plot_data_T_c.RData")
# load("FR data/Model predictions with simulated parameters/ANOVA-like fits/1000 simulations/FR_ANOVA_sim_plot_data_T_c.RData")




# ------------------ ST ------------------
ggplot(predictions_dd_ST, aes(x = N, y = Ndead, group = Simulation_nr)) +
  geom_line(linewidth = .5, alpha = .3)

# There are some predictions that look weird

# All columns that have biologically impossible predictions will be NAs
predictions_filtered_ST <- apply(predictions_ST, 2, function(col) {
  if (sum(col < 0) > 0) {
    return(rep(NA, length(col)))
  } else if (sum(col > 2000) > 0) { # Value used to exclude only really weird fits
    return(rep(NA, length(col)))
  } else {
    return(col)
  }
})

predictions_filtered_dd_ST <- data.frame(predictions_filtered_ST) %>%
  pivot_longer(cols = 1:n_sim, names_to = "Simulation_nr", values_to = "Ndead")
predictions_filtered_dd_ST$N <- rep(N_ST, each = n_sim)

predictions_filtered_dd_ST$b.log.estimate <- rep(param_sim[,"bC.log"] + param_sim[,"bS.log"] + param_sim[,"bT.log"] + param_sim[,"bST.log"], times = length(N_ST))
predictions_filtered_dd_ST$h.estimate <- rep(param_sim[,"hC"] + param_sim[,"hS"] + param_sim[,"hT"] + param_sim[,"hST"], times = length(N_ST))
predictions_filtered_dd_ST$m.estimate <- rep(param_sim[,"mC"] + param_sim[,"mS"] + param_sim[,"mT"] + param_sim[,"mST"], times = length(N_ST))

ggplot(predictions_filtered_dd_ST, aes(x = N, y = Ndead, group = Simulation_nr)) +
  geom_line(linewidth = .5, alpha = .3)


# Compute mean and confidence intervals
response_mean_ST <- apply(predictions_filtered_ST, 1, mean, na.rm = T)
response_lower_ST <- apply(predictions_filtered_ST, 1, quantile, probs = 0.025, na.rm = T)
response_upper_ST <- apply(predictions_filtered_ST, 1, quantile, probs = 0.975, na.rm = T)

# Visualize
plot(N_ST, response_mean_ST, type = "l", lwd = 2, col = "blue",
     xlab = "Resource Density (N)", ylab = "Functional Response",
     ylim = range(c(response_lower_ST, response_upper_ST)),
     main = "Functional Response with Confidence Intervals - SEL S")
polygon(c(N_ST, rev(N_ST)),
        c(response_lower_ST, rev(response_upper_ST)),
        col = rgb(0, 0, 1, 0.2), border = NA)


# Manually checking the biologically weird predictions
weird_fits_IDs_ST <- c()
for(i in 1:dim(predictions_filtered_ST)[2]) {
  
  if(sum(is.na(predictions_filtered_ST[,i])) == length(predictions_filtered_ST[,i])) {
    weird_fits_IDs_ST <- c(weird_fits_IDs_ST, i)
  }
  
} # NULL

predictions_ST_c <- predictions_ST
predictions_dd_ST_c <- predictions_dd_ST
predictions_filtered_ST_c <- predictions_filtered_ST
predictions_filtered_dd_ST_c <- predictions_filtered_dd_ST
response_mean_ST_c <- response_mean_ST
response_lower_ST_c <- response_lower_ST
response_upper_ST_c <- response_upper_ST

save(predictions_ST_c, predictions_dd_ST_c,
     predictions_filtered_ST_c, predictions_filtered_dd_ST_c,
     response_mean_ST_c, response_lower_ST_c, response_upper_ST_c,
     file = "FR data/Model predictions with simulated parameters/ANOVA-like fits/1000 simulations/FR_ANOVA_multiple_model_predictions_ST_c.RData")
# load("FR data/Model predictions with simulated parameters/ANOVA-like fits/1000 simulations/FR_ANOVA_multiple_model_predictions_ST_c.RData")













# Compute the end resource density for all starting resource densities and varying parameters
### Output will be a matrix where each column is the estimated end resource density of
### one of the parameter combinations
# time1 <- Sys.time()
# predictionLV_ST <- apply(param_sim, 1, function(p) {
#   Nend <- unlist(lapply(N_ST, FUN = function(n) {
#     eq.odeint.general.pred.t2.mP(start.v = c(n, 5),
#                                  b = exp(p[["bC.log"]]),
#                                  h = p[["hC"]],
#                                  r = p[["rC"]],
#                                  K = exp(p[["KC.log"]]),
#                                  c = p[["c"]],
#                                  m = p[["mC"]],
#                                  Tt = 1,
#                                  timesteplength = 1/100)[100+1,2]
#   }))
# })
# Sys.time() - time1

# In parallel
library(doParallel)
bigstatsr::nb_cores()

my.cluster <- parallel::makeCluster(
  bigstatsr::nb_cores(), 
  type = "PSOCK"
)

doParallel::registerDoParallel(cl = my.cluster)
foreach::getDoParRegistered()
foreach::getDoParWorkers()

predictionLV_ST <- NA # name of file where parallel runs will be save

(t0 <- Sys.time()) # timer
mcoptions <- list(preschedule=FALSE)

predictionLV_ST <- foreach(i = 1:nrow(param_sim),
                           .combine = 'cbind',
                           .options.multicore = mcoptions,
                           .packages = c('bbmle', 'odeintr', 'tidyverse', 'R.utils'),
                           .errorhandling = "pass") %dopar%
  {
    compile_sys("FRII_gen_pred_mP", FRII.gen.pred.mP, pars = c("b","h","r","K","c","m"), method = "rk54")
    
    p <- param_sim[i,]
    
    predictionLV_column <- unlist(lapply(N_ST, FUN = function(n) {
      eq.odeint.general.pred.t2.mP(start.v = c(n, 5),
                                   b = exp(p[["bC.log"]] + p[["bS.log"]] + p[["bT.log"]] + p[["bST.log"]]),
                                   h = p[["hC"]] + p[["hS"]] + p[["hT"]] + p[["hST"]],
                                   r = p[["rST"]],
                                   K = exp(p[["KST.log"]]),
                                   c = p[["c"]],
                                   m = p[["mC"]] + p[["mS"]] + p[["mT"]] + p[["mST"]],
                                   Tt = 1,
                                   timesteplength = 1/100)[100+1,2]
    }))
    
    return(predictionLV_column)
  }

parallel::stopCluster(cl = my.cluster)
stopImplicitCluster()

tend<-Sys.time()
tend-t0 
# beep(sound = 2, expr = NULL)






predictionLV_dd_ST <- data.frame(predictionLV_ST) %>%
  pivot_longer(cols = 1:n_sim, names_to = "Simulation_nr", values_to = "Nend")
predictionLV_dd_ST$N <- rep(N_ST, each = n_sim)
predictionLV_dd_ST$b.log.estimate <- rep(param_sim[,"bC.log"] + param_sim[,"bS.log"] + param_sim[,"bT.log"] + param_sim[,"bST.log"], times = length(N_ST))
predictionLV_dd_ST$h.estimate <- rep(param_sim[,"hC"] + param_sim[,"hS"] + param_sim[,"hT"] + param_sim[,"hST"], times = length(N_ST))
predictionLV_dd_ST$m.estimate <- rep(param_sim[,"mC"] + param_sim[,"mS"] + param_sim[,"mT"] + param_sim[,"mST"], times = length(N_ST))


predictionLV_dd_ST$Simulation_nr <- gsub(predictionLV_dd_ST$Simulation_nr, pattern = "result.", replacement = "X")

ggplot() +
  geom_line(aes(x = N, y = Nend, group = Simulation_nr, color = b.log.estimate), filter(predictionLV_dd_ST, !(Simulation_nr %in% paste0("X", weird_fits_IDs_ST))), size = .5, alpha = .8) +
  geom_point(aes(x = Resource_start_density, y = Resource_end_density), data = filter(FR_all, SEL_environment == "C", Consumer_presence == 1))


# Exclude weird fits from this dataset too
predictionLV_filtered_ST <- predictionLV_ST


# Compute mean and confidence intervals
LV_mean_ST <- apply(predictionLV_filtered_ST, 1, mean, na.rm = T)
LV_lower_ST <- apply(predictionLV_filtered_ST, 1, quantile, probs = 0.025, na.rm = T)
LV_upper_ST <- apply(predictionLV_filtered_ST, 1, quantile, probs = 0.975, na.rm = T)

predictionLV_ST_c <- predictionLV_ST
predictionLV_filtered_ST_c <- predictionLV_filtered_ST
predictionLV_dd_ST_c <- predictionLV_dd_ST
LV_mean_ST_c <- LV_mean_ST
LV_lower_ST_c <- LV_lower_ST
LV_upper_ST_c <- LV_upper_ST

save(predictionLV_ST_c, predictionLV_filtered_ST_c, predictionLV_dd_ST_c,
     LV_mean_ST_c, LV_lower_ST_c, LV_upper_ST_c,
     file = "FR data/Model predictions with simulated parameters/ANOVA-like fits/1000 simulations/FR_ANOVA_multiple_LV_predictions_ST_c.RData")
# load("FR data/Model predictions with simulated parameters/ANOVA-like fits/1000 simulations/FR_ANOVA_multiple_LV_predictions_ST_c.RData")

predictions_dd_ST_c <- predictions_dd_ST

save(predictionLV_dd_ST_c, predictions_dd_ST_c,
     file = "FR data/Model predictions with simulated parameters/ANOVA-like fits/1000 simulations/FR_ANOVA_sim_plot_data_ST_c.RData")
# load("FR data/Model predictions with simulated parameters/ANOVA-like fits/1000 simulations/FR_ANOVA_sim_plot_data_ST_c.RData")



