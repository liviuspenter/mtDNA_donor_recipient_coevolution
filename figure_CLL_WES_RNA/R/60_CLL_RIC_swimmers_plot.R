library(dplyr)
library(ggplot2)
library(tidyr)

swimmer.data <- readxl::read_excel("./data/CLL_RIC/20231229_CLL_RIC_mtDNA_meta.xlsx", sheet = 3)

swimmer.data$sample.to.sample <- swimmer.data$Sample_to_HSCT + swimmer.data$HSCT_to_relapse_sample
swimmer.data <- swimmer.data[order(swimmer.data$sample.to.sample, decreasing = F), ]

swimmer.data$Patient <- factor(as.character(swimmer.data$Patient), levels = as.character(swimmer.data$Patient))

ggplot() +
  geom_col(data = swimmer.data, aes(x = Patient, y = sample.to.sample), fill = "grey50") +
  geom_col(data = swimmer.data, aes(x = Patient, y = Sample_to_HSCT), fill = "grey90") +
  geom_col(data = swimmer.data, aes(x = Patient, y = sample.to.sample), fill = NA, color = "black") +
  geom_point(data = swimmer.data, aes(x = Patient, y = 0), color = "blue") +
  geom_point(data = swimmer.data, aes(x = Patient, y = sample.to.sample), color = "firebrick") +
  scale_x_discrete("RIC CLL") +
  scale_y_continuous("Days from baseline sample") +
  coord_flip() +
  theme_classic() +
  theme(
    axis.title = element_text("Arial", size = 10, color = "black"),
    axis.text = element_text("Arial", size = 10, color = "black")
  )
ggsave("./figure_CLL_WES_RNA/figures/CLL_RIC/20240126_CLL_RIC_swimmers_plot.svg", width = 2.5, height = 1.7)
