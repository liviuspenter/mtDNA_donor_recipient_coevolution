library(ggplot2)

mito.mutation.order <- read.table("./data/AML_coevolution/AML1026/20240624_AML1026_combined_markers.csv", sep = "\t")

# perform sanity check by looking at VMR and strand coordination in ASAP-seq data
high.confidence.variants.asap.1026 <- data.table::fread("./data/10026/mtDNA/20210429_AML1026_high_confidence_variants.csv")
rownames(high.confidence.variants.asap.1026) <- high.confidence.variants.asap.1026$variant
high.confidence.variants.asap.1026$highlight <- ifelse(high.confidence.variants.asap.1026$variant %in% mito.mutation.order$x, "yes", "no")
high.confidence.variants.asap.1026$highlight <- factor(high.confidence.variants.asap.1026$highlight, levels = c("no", "yes"))

mito.mutation.order$x %in% high.confidence.variants.asap.1026$variant

ggplot() +
  geom_vline(xintercept = 0.65) +
  geom_hline(yintercept = 0.01) +
  geom_point(
    data = high.confidence.variants.asap.1026[which(high.confidence.variants.asap.1026$highlight == "no"), ],
    aes(x = strand_correlation, y = vmr, color = highlight), size = 0.5
  ) +
  geom_point(
    data = high.confidence.variants.asap.1026[which(high.confidence.variants.asap.1026$highlight == "yes"), ],
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
ggsave("./figure_AML_coevolution/figures/AML1026/20240624_sanity_check_mutations.svg", width = 2, height = 2)
