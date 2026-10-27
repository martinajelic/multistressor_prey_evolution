# ----------------------------------------------------------------------------------------------------------------------------------
# Supplement S.6: Correlation of prey traits and functional response parameters
# ----------------------------------------------------------------------------------------------------------------------------------
# 
# The following code answers the additional research question: do evolutionary trait changes in the
# prey species correlate to changes in functional response parameters?
# 
# NOTE: this script was run on R version 4.5.2 (provided in folder
# "R versions/R version for data analysis"). It will not run on R version 4.3.2 due to package
# compatibility issues.
# 
# --------------------------------------------------------------------------------------------------


rm(list = ls())

library(tidyverse)
library(generics)
library(broom)
library(patchwork)
library(cowplot)

# Raw data
load("Data/Traits_assay.RData")

# Model fitting to FR assay data was done in "R scripts/02_FR_search_grid.R"
load("Data/FR assay/03 - FR search grid/FR_ANOVA_like_best_fit.RData")

# Coefficient computation per selection environment was done in
# "R scripts/04_Data_analysis_and_results.R"
load("Data/FR assay/04 - FR fits and simulations/Best_fit_estimates_per_environment.RData")

(FR_coefs <- generics::tidy(best_FR_fit))

# Merging all coefficient estimates for regressions
FR_all <- rbind(FRC_coefs %>% mutate(Selection_environment = "C"),
                FRS_coefs %>% mutate(Selection_environment = "S"),
                FRT_coefs %>% mutate(Selection_environment = "T"),
                FRST_coefs %>% mutate(Selection_environment = "ST"))
FR_all$Selection_environment <- factor(FR_all$Selection_environment,
                                       levels = c("C", "S", "T", "ST"))

FR_all$term <- sub("S", "", FR_all$term)
FR_all$term <- sub("T", "", FR_all$term)
FR_all$term <- sub("C", "", FR_all$term)

FR_all_wide <- pivot_wider(FR_all,
                           id_cols = c(Selection_environment),
                           values_from = estimate,
                           names_from = c(term))

# Keep only trait data from assay environment C
relevant_envs <- dplyr::filter(traits_assay, Assay_environment == "C")

# Mean cell area and mean movement speed per selection environment
FR_traits_comparison <- merge(relevant_envs, FR_all_wide, by = "Selection_environment") %>%
  group_by(Selection_environment, Selection_T, Selection_S,
           Assay_T, Assay_S, Assay_environment, b.log, h) %>%
  summarise(Mean_cell_area = mean(Mean_cell_area),
            Mean_movement_speed = mean(Mean_movement_speed))


# --------------- Movement speed and space clearance rate
blog_speed_mod <- lm(b.log ~ Mean_movement_speed,
                     FR_traits_comparison)
summary(blog_speed_mod)

cor.test(FR_traits_comparison$Mean_movement_speed, FR_traits_comparison$b.log)


# --------------- Cell size and space clearance rate
blog_size_mod <- lm(b.log ~ Mean_cell_area,
                    FR_traits_comparison)
summary(blog_size_mod)

cor.test(FR_traits_comparison$Mean_cell_area, FR_traits_comparison$b.log)


# --------------- Cell size and handling time
h_size_mod <- lm(h ~ Mean_cell_area,
                 FR_traits_comparison)
summary(h_size_mod)

cor.test(FR_traits_comparison$Mean_cell_area, FR_traits_comparison$h)




# --------------------------------------------------------------------------------------------------
# Figure S4: FR parameters and prey traits
# --------------------------------------------------------------------------------------------------
newdat <- data.frame(
  Mean_movement_speed = seq(min(FR_traits_comparison$Mean_movement_speed),
                            max(FR_traits_comparison$Mean_movement_speed),
                            length.out = 200),
  Mean_cell_area = seq(min(FR_traits_comparison$Mean_cell_area),
                       max(FR_traits_comparison$Mean_cell_area),
                       length.out = 200)
)

Pred_blog_speed <- cbind(newdat,
                         predict(blog_speed_mod, newdat, interval = "confidence", level = 0.95))
Pred_blog_size <- cbind(newdat,
                        predict(blog_size_mod, newdat, interval = "confidence", level = 0.95))
Pred_h_size <- cbind(newdat,
                     predict(h_size_mod, newdat, interval = "confidence", level = 0.95))

# One plot for each regression
plot_b_speed <- ggplot() +
    geom_point(data = FR_traits_comparison, aes(x = Mean_movement_speed, y = b.log, color = Selection_environment)) +
    scale_color_manual(values = c("#716868", "#E69F00", "#56B4E9", "#009E73")) +
    geom_line(data = Pred_blog_speed, aes(x = Mean_movement_speed, y = fit)) +
  geom_ribbon(data = Pred_blog_speed,
                aes(x = Mean_movement_speed, y = fit, ymin = lwr, ymax = upr),
                alpha = 0.1) +
    theme_light() +
    labs(x = latex2exp::TeX("Movement speed ($µm$ $sec^{-1}$)"),
         y = latex2exp::TeX("$ln(b)$"),
         color = "Selection environment",
         shape = "Biological replicate") +
    theme(legend.position = "none")

plot_b_size <- ggplot() +
    geom_point(data = FR_traits_comparison, aes(x = Mean_cell_area, y = b.log, color = Selection_environment)) +
    scale_color_manual(values = c("#716868", "#E69F00", "#56B4E9", "#009E73")) +
    geom_line(data = Pred_blog_size, aes(x = Mean_cell_area, y = fit)) +
    geom_ribbon(data = Pred_blog_size,
                aes(x = Mean_cell_area, y = fit, ymin = lwr, ymax = upr),
                alpha = 0.1) +
    theme_light() +
    labs(x = latex2exp::TeX("Cell area $(µm^2)$"),
         y = latex2exp::TeX("$ln(b)$"),
         color = "Selection environment",
         shape = "Biological replicate") +
    theme(legend.position = "none")

plot_h_size <- ggplot() +
    geom_point(data = FR_traits_comparison, aes(x = Mean_cell_area, y = h, color = Selection_environment)) +
    scale_color_manual(values = c("#716868", "#E69F00", "#56B4E9", "#009E73")) +
    geom_line(data = Pred_h_size, aes(x = Mean_cell_area, y = fit)) +
    geom_ribbon(data = Pred_h_size,
                aes(x = Mean_cell_area, y = fit, ymin = lwr, ymax = upr),
                alpha = 0.1) +
    theme_light() +
    labs(x = latex2exp::TeX("Cell area $(µm^2)$"),
         y = latex2exp::TeX("$h$"),
         color = "Selection environment",
         shape = "Biological replicate") +
    theme(legend.position = "none")

# Top row with plots
facet_plots <- (plot_b_speed | plot_b_size | plot_h_size) +
  plot_layout(ncol = 3, nrow = 1)  +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold", margin = margin(t = -20, r = 0, b = 0, l = 0))) &
  theme(plot.margin = margin(10, 0, 10, 10)) &  # top, right, bottom, left
  theme(axis.text = element_text(size = 12))

facet_row <- plot_grid(facet_plots, ncol = 1)

# Bottom row with legends
legend_FR_traits <- cowplot::get_legend(
  plot_b_speed + theme(legend.position = "bottom",
                       # legend.text = element_text(size = 12),
                       legend.title = element_text(hjust = 0.5),
                       legend.title.position = "top")
)

legend_row <- plot_grid(legend_FR_traits,
                    ncol = 1) +
  theme(plot.margin = margin(0, 0, 0, 0))

# Figure 4
combined <- plot_grid(facet_row, legend_row, ncol = 1, rel_heights = c(1, 0.1)) +
  theme(plot.margin = margin(0, 10, 0, -10))  # top, right, bottom, left

cor_b_speed <- cor.test(FR_traits_comparison$Mean_movement_speed, FR_traits_comparison$b.log)
cor_b_size <- cor.test(FR_traits_comparison$Mean_cell_area, FR_traits_comparison$b.log)
cor_h_size <- cor.test(FR_traits_comparison$Mean_cell_area, FR_traits_comparison$h)

ggdraw(combined) +
    draw_label(paste0("r = ", round(cor_b_speed$estimate, digits = 4), ", p = ", round(cor_b_speed$p.value, digits = 4), " *"),
               x = 0.19, y = 0.88, vjust = 0, size = 10) +
    draw_label(paste0("r = ", round(cor_b_size$estimate, digits = 4), ", p = ", round(cor_b_size$p.value, digits = 4), ""),
               x = 0.515, y = 0.88, vjust = 0, size = 10) +
    draw_label(paste0("r = ", round(cor_h_size$estimate, digits = 4), ", p = ", round(cor_h_size$p.value, digits = 4), ""),
               x = 0.837, y = 0.88, vjust = 0, size = 10)

# ggsave(filename = "Figures/Figure_S4.pdf",
#        height = 3.6,
#        width = 9.6)
