# ----------------------------------------------------------------------------------------------------------------------------------
# Data analysis
# ----------------------------------------------------------------------------------------------------------------------------------


# --------------------------------------------------------------------------------------------------
# Traits assay
# --------------------------------------------------------------------------------------------------

rm(list = ls())

library(tidyverse)

load("Data/Traits_assay.RData")

ggplot() +
  geom_boxplot(data = traits_assay, aes(x = Assay_environment, y = GR,
                                                          color = Selection_environment, fill = Selection_environment),
               alpha = 0.3, linewidth = 0.2) +
  geom_point(data = traits_assay, aes(x = Assay_environment, y = GR,
                                                        color = Selection_environment, fill = Selection_environment, shape = Biorep),
             alpha = 1, size = 1, position = position_jitterdodge(jitter.width = 0.09, dodge.width = 0.75)) +
  theme_bw() +
  labs(x = "Assay environment",
       y = latex2exp::TeX("Growth rate $ln(N_{t}/N_{0})/\\Delta t$"),
       color = "Selection environment",
       fill = "Selection environment",
       shape = "Replicate") +
  scale_color_manual(values = c("#000000", "#E69F00", "#56B4E9", "#009E73")) +
  scale_fill_manual(values = c("#000000", "#E69F00", "#56B4E9", "#009E73")) +
  geom_hline(yintercept = 0, linetype = 2, color = "gray") +
  theme(strip.background = element_blank())
