chimerism.corr <- readxl::read_excel("./data/10026/mtDNA/20240604_correlation_chimerism.xlsx")

cor.test(chimerism.corr$Clinical.chimerism, chimerism.corr$mtDNA.chimerism)

ggsave("./figure_10026/figures/plots/20240604_correlation_chimerism.svg", width = 1.8, height = 1.8)
ggplot(data = chimerism.corr, aes(x = 100 * Clinical.chimerism, y = 100 * mtDNA.chimerism)) +
  geom_abline(slope = 1) +
  geom_point(size = 0.5) +
  scale_x_continuous("% clinical T cell chimerism") +
  scale_y_continuous("% mtDNA single T cell chimerism") +
  theme_classic() +
  theme(
    axis.title = element_text("Arial", size = 10, color = "black"),
    axis.text = element_text("Arial", size = 10, color = "black")
  )
dev.off()
