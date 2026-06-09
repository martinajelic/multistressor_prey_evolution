# ----------------------------------------------------------------------------------------------------------------------------------
# Data analysis and results
# ----------------------------------------------------------------------------------------------------------------------------------
# 
# The following code yields the two manuscript figures:
# 
# - Figure 3: prey growth rates and traits
# - Figure 4: functional response parameters and fits
# 
# See figure labels in the main text for detailed descriptions.
# 
# NOTE: this script was run on R version 4.5.2 (provided in folder
# "R versions/R version for data analysis"). It will not run correctly on R version 4.3.2 due to 
# package compatibility issues.
# 
# ----------------------------------------------------------------------------------------------------------------------------------

# --------------------------------------------------------------------------------------------------
# Figure 3: prey traits and growth rate
# --------------------------------------------------------------------------------------------------

rm(list = ls())

library(tidyverse)
library(generics)
library(broom)
library(patchwork)
library(cowplot)

load("Data/Traits_assay.RData")
load("Data/Marginal_means_estimates.RData")

# Marginal means and sd for each trait at each selection x assay environment combination
dd_GR <- as.data.frame(summary(emm_GR)) %>%
  dplyr::select(Assay_environment, Selection_environment, emmean, SE)
dd_size <- as.data.frame(summary(emm_size)) %>%
  dplyr::select(Assay_environment, Selection_environment, emmean, SE)
dd_speed <- as.data.frame(summary(emm_speed)) %>%
  dplyr::select(Assay_environment, Selection_environment, emmean, SE)

# Plotting the data for each response variable
pos <- position_jitterdodge(jitter.width = 0.1, dodge.width = 0.75)

P_GR <- ggplot() +
  geom_point(data = traits_assay, aes(x = Assay_environment, y = GR,
                                      color = Selection_environment, fill = Selection_environment, shape = Biorep),
             alpha = 0.2, size = 1, position = pos) +
  scale_shape_manual(values = 22:26) +
  geom_errorbar(data = dd_GR, aes(x = Assay_environment, y = emmean, ymin = emmean - SE, ymax = emmean + SE,
                                  color = Selection_environment),
                linewidth = 0.5, position = pos, width = 0.4) +
  geom_point(data = dd_GR, aes(x = Assay_environment, y = emmean,
                               color = Selection_environment),
             size = 1.5, position = pos) +
  theme_light() +
  labs(x = "Assay environment",
       y = latex2exp::TeX("Prey growth rate ($day^{-1}$)"),
       color = "Selection environment",
       fill = "Selection environment",
       shape = "Biological replicate") +
  scale_color_manual(values = c("#000000", "#E69F00", "#56B4E9", "#009E73")) +
  scale_fill_manual(values = c("#000000", "#E69F00", "#56B4E9", "#009E73")) +
  theme(legend.position = "none") +
  geom_hline(yintercept = 0, linetype = 2, color = "gray") +
  theme(strip.background = element_blank())

P_size <- ggplot() +
  geom_point(data = traits_assay, aes(x = Assay_environment, y = Mean_cell_area,
                                      color = Selection_environment, fill = Selection_environment, shape = Biorep),
             alpha = 0.2, size = 1, position = pos) +
  scale_shape_manual(values = 22:26) +
  geom_errorbar(data = dd_size, aes(x = Assay_environment, y = emmean, ymin = emmean - SE, ymax = emmean + SE,
                                    color = Selection_environment),
                linewidth = 0.5, position = pos, width = 0.4) +
  geom_point(data = dd_size, aes(x = Assay_environment, y = emmean,
                                 color = Selection_environment),
             size = 1.5, position = pos) +
  theme_light() +
  labs(x = "Assay environment",
       y = latex2exp::TeX("Prey cell area $(µm^2)$"),
       color = "Selection environment",
       fill = "Selection environment",
       shape = "Biological replicate") +
  scale_color_manual(values = c("#000000", "#E69F00", "#56B4E9", "#009E73")) +
  scale_fill_manual(values = c("#000000", "#E69F00", "#56B4E9", "#009E73")) +
  theme(legend.position = "none") +
  theme(strip.background = element_blank())

P_speed <- ggplot() +
  geom_point(data = traits_assay, aes(x = Assay_environment, y = Mean_movement_speed,
                                      color = Selection_environment, fill = Selection_environment, shape = Biorep),
             alpha = 0.2, size = 1, position = pos) +
  scale_shape_manual(values = 22:26) +
  geom_errorbar(data = dd_speed, aes(x = Assay_environment, y = emmean, ymin = emmean - SE, ymax = emmean + SE,
                                     color = Selection_environment),
                linewidth = 0.5, position = pos, width = 0.4) +
  geom_point(data = dd_speed, aes(x = Assay_environment, y = emmean,
                                  color = Selection_environment),
             size = 1.5, position = pos) +
  theme_light() +
  labs(x = "Assay environment",
       y = latex2exp::TeX("Prey mov. speed ($µm$ $sec^{-1}$)"),
       color = "Selection environment",
       fill = "Selection environment",
       shape = "Biological replicate") +
  scale_color_manual(values = c("#000000", "#E69F00", "#56B4E9", "#009E73")) +
  scale_fill_manual(values = c("#000000", "#E69F00", "#56B4E9", "#009E73")) +
  theme(legend.position = "none") +
  theme(strip.background = element_blank())

# Extract the common legend
legend_growth <- get_legend(
  P_GR + theme(legend.position = "bottom",
               # legend.text = element_text(size = 12),
               legend.title = element_text(size = 14, hjust = 0.5),
               legend.title.position = "top")
)
legend_row <- plot_grid(legend_growth,
                        ncol = 1) +
  theme(plot.margin = margin(0, 0, 0, 0))

# Arrange the grid
facet_plots <- (P_GR +
                  theme(axis.title.y = element_text(hjust = 0.5, vjust = 0.5)) | P_size | P_speed) +
    plot_layout(ncol = 3, nrow = 1)  +
    plot_annotation(tag_levels = "a") &
    theme(plot.tag = element_text(face = "bold",
                                  margin = margin(t = -20, r = 0, b = 0, l = 0))) &
    theme(plot.margin = margin(10, 20, 0, 0)) &  # top, right, bottom, left
    theme(axis.text = element_text(size = 12),
          axis.title.x = element_text(size = 12),
          axis.title.y = element_text(size = 12))
facet_row <- plot_grid(facet_plots, ncol = 1)

# Combine the two
plot_grid(legend_row, facet_row, ncol = 1, rel_heights = c(0.18, 1)) +
    theme(plot.margin = margin(0, 20, 0, 40))  # top, right, bottom, left

# ggsave(filename = "Figures/Figure_3.pdf",
#        width = 10, height = 3.25)
































# --------------------------------------------------------------------------------------------------
# Figure 4: FR parameters and fits
# --------------------------------------------------------------------------------------------------

rm(list = ls())

library(tidyverse)
library(generics)
library(broom)

# ---------------------- Loading relevant datasets ----------------------
# Model predictions from best model fit (created in 03_FR_fits_simulations.R)
load("Data/FR assay/04 - FR fits and simulations/Best_fit_curves_and_95%CIs.RData")

# Predicted values from best model fit (created in 04_Data_analysis_and_results.R)
load("Data/FR assay/04 - FR fits and simulations/Best_fit_estimates_per_environment.RData")

# ---------------------- Computing estimated parameters and uncertainty ----------------------
# Creating datasets with parameter estimates and 95% CIs for plotting
param_df <- rbind(FRC_coefs %>% dplyr::select(term, estimate, std.error) %>% mutate(SEL_environment = factor("C", levels = c("C", "S", "T", "ST"))),
                  FRS_coefs %>% dplyr::select(term, estimate, std.error) %>% mutate(SEL_environment = factor("S", levels = c("C", "S", "T", "ST"))),
                  FRT_coefs %>% dplyr::select(term, estimate, std.error) %>% mutate(SEL_environment = factor("T", levels = c("C", "S", "T", "ST"))),
                  FRST_coefs %>% dplyr::select(term, estimate, std.error) %>% mutate(SEL_environment = factor("ST", levels = c("C", "S", "T", "ST")))) %>%
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

# Creating a new dataset with the b.log variable changed to exp(b.log) = b
param_df_exp <- param_df

# Exp the b.log values
param_df_exp$Estimate <- ifelse(param_df_exp$Parameter == "b.log", exp(param_df_exp$Estimate), param_df_exp$Estimate)

# Change b.log label to b
param_df_exp$Parameter <- ifelse(param_df_exp$Parameter == "b.log", "b", param_df_exp$Parameter)

# Compute upper and lower values of Std. error bars for exp(b.log) values
param_df_exp$StdUpper <- ifelse(param_df_exp$Parameter == "b",
                                exp(param_df_exp$StdUpper),
                                param_df_exp$StdUpper)
param_df_exp$StdLower <- ifelse(param_df_exp$Parameter == "b",
                                exp(param_df_exp$StdLower),
                                param_df_exp$StdLower)

# Compute upper and lower values of CIs for exp(b.log) values
param_df_exp$CI_lower <- ifelse(param_df_exp$Parameter == "b",
                                exp(param_df_exp$CI_lower),
                                param_df_exp$CI_lower)
param_df_exp$CI_upper <- ifelse(param_df_exp$Parameter == "b",
                                exp(param_df_exp$CI_upper),
                                param_df_exp$CI_upper)

# ---------------------- Figure 1 ----------------------
Pb <- ggplot(filter(param_df_exp, Parameter == "b"), aes(x = SEL_environment, y = Estimate, color = SEL_environment)) +
  geom_point(size = 3, position = position_dodge(width = 0.7)) +
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper), width = 0.2, size = 0.7,
                position = position_dodge(width = 0.7)) +
  theme_light() +
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
  theme_light() +
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
  theme_light() +
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
  theme_light() +
  scale_color_manual(values = c("#716868", "#E69F00")) +
  scale_fill_manual(values = c("#716868", "#E69F00")) +
  labs(x = latex2exp::TeX("Prey start density ($mL^{-1}$)"),
       y = latex2exp::TeX("Prey eaten ($predator^{-1} day^{-1}$)"),
       color = "Selection environment",
       fill = "Selection environment") +
  theme(strip.background = element_blank(),
        legend.position = "bottom") +
  theme(plot.margin = margin(t = 10, r = 5, b = 0, l = 20, "pt")) +
  ylim(0, 85)


PFR_CT <- ggplot(filter(FR_predictions, Prey_start_density < 1000,
                        SEL_environment %in% c("C", "T")), aes(x = Prey_start_density, y = Prey_consumed, color = SEL_environment, group = SEL_environment)) +
  geom_line(size = 1.1) +
  geom_ribbon(aes(x = Prey_start_density, y = Prey_consumed, ymin = CI_lower, ymax = CI_upper, fill = SEL_environment, color = NULL), alpha = 0.15) +
  theme_light() +
  scale_color_manual(values = c("#716868", "#56B4E9")) +
  scale_fill_manual(values = c("#716868", "#56B4E9")) +
  labs(x = latex2exp::TeX("Prey start density ($mL^{-1}$)"),
       y = latex2exp::TeX("Prey eaten ($predator^{-1} day^{-1}$)"),
       color = "Selection environment",
       fill = "Selection environment") +
  theme(strip.background = element_blank(),
        legend.position = "bottom") +
  theme(plot.margin = margin(t = 10, r = 5, b = 0, l = 20, "pt")) +
  ylim(0, 85)

PFR_CST <- ggplot(filter(FR_predictions, Prey_start_density < 1000,
                         SEL_environment %in% c("C", "ST")), aes(x = Prey_start_density, y = Prey_consumed, color = SEL_environment, group = SEL_environment)) +
  geom_line(size = 1.1) +
  geom_ribbon(aes(x = Prey_start_density, y = Prey_consumed, ymin = CI_lower, ymax = CI_upper, fill = SEL_environment, color = NULL), alpha = 0.15) +
  theme_light() +
  scale_color_manual(values = c("#716868", "#009E73")) +
  scale_fill_manual(values = c("#716868", "#009E73")) +
  labs(x = latex2exp::TeX("Prey start density ($mL^{-1}$)"),
       y = latex2exp::TeX("Prey eaten ($predator^{-1} day^{-1}$)"),
       color = "Selection environment",
       fill = "Selection environment") +
  theme(strip.background = element_blank(),
        legend.position = "bottom") +
  theme(plot.margin = margin(t = 10, r = 5, b = 0, l = 20, "pt")) +
  ylim(0, 85)

ggpubr::ggarrange(Pb, Ph, Pm, PFR_CS, PFR_CT, PFR_CST, ncol = 3, nrow = 2, labels = "auto", common.legend = T)

# ggsave(filename = "Figures/Figure_4.pdf",
#        width = 10, height = 6)




