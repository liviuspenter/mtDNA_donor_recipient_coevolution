library(ggplot2)

mito.mutation.order <- read.table("./data/AML_coevolution/AML1010/20240624_AML1010_combined_markers.csv", sep = "\t")

# perform sanity check by looking at VMR and strand coordination in ASAP-seq data
high.confidence.variants.asap.1010 <- data.table::fread("./data/10026/mtDNA/20210429_AML1010_high_confidence_variants.csv")
rownames(high.confidence.variants.asap.1010) <- high.confidence.variants.asap.1010$variant
high.confidence.variants.asap.1010$highlight <- ifelse(high.confidence.variants.asap.1010$variant %in% mito.mutation.order$x, "yes", "no")
high.confidence.variants.asap.1010$highlight <- factor(high.confidence.variants.asap.1010$highlight, levels = c("no", "yes"))

ggplot() +
  geom_vline(xintercept = 0.65) +
  geom_hline(yintercept = 0.01) +
  geom_point(
    data = high.confidence.variants.asap.1010[which(high.confidence.variants.asap.1010$highlight == "no"), ],
    aes(x = strand_correlation, y = vmr, color = highlight), size = 0.5
  ) +
  geom_point(
    data = high.confidence.variants.asap.1010[which(high.confidence.variants.asap.1010$highlight == "yes"), ],
    aes(x = strand_correlation, y = vmr, color = highlight), size = 0.5
  ) +
  scale_x_continuous("strand correlation", limits = c(-1, 1)) +
  scale_y_log10("VMR", limits = c(0.001, 3)) +
  scale_color_manual(values = c("yes" = "firebrick", "no" = "grey")) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black")
  )
ggsave("./figure_AML_coevolution/figures/AML1010/20240624_sanity_check_mutations.svg", width = 2, height = 2)
