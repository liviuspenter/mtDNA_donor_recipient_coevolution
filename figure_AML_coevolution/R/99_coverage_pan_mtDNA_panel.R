library(dplyr)
library(ggplot2)

GenePos.tib <- tibble(
  Names = c(
    "ATP6", "ATP8", "COX1", "COX2", "COX3", "CYTB", "ND1", "ND2", "ND3",
    "ND4", "ND4L", "ND5", "ND6", "RNR1", "RNR2"
  ),
  start = c(8527, 8366, 5904, 7586, 9207, 14747, 3307, 4470, 10059, 10760, 10470, 12337, 14149, 648, 1671),
  end = c(9207, 8572, 7445, 8269, 9990, 15887, 4262, 5511, 10404, 12137, 10766, 14148, 14673, 1601, 3229)
)
GenePos.tib$Names <- factor(GenePos.tib$Names, levels = GenePos.tib$Names)

coverage.df <- data.table::fread("./data/AML_coevolution/AML1/AML_1.coverage.txt.gz")
colnames(coverage.df) <- c("chrM", "bc", "depth")

boo <- coverage.df %>%
  group_by(chrM) %>%
  summarize(mean.depth = mean(depth))

y.max <- max(boo$mean.depth)
GenePos.tib <- GenePos.tib %>%
  arrange(start) %>%
  mutate(mid = round((end - start) / 2 + start, 0), ycoord = rep(c(
    1.3 * y.max,
    1.2 * y.max,
    1.1 * y.max
  ), length.out = 15))

for (i in seq(1, nrow(GenePos.tib))) {
  boo[which(between(boo$chrM, GenePos.tib$start[i], GenePos.tib$end[i])), "gene"] <- as.character(GenePos.tib$Names[i])
}
boo$gene <- factor(boo$gene, levels = GenePos.tib$Names)

col.list <- as.character(BuenColors::jdb_palette(name = "corona", n = 15))
names(col.list) <- levels(GenePos.tib$Names)

p <- ggplot(boo, aes(x = chrM, y = mean.depth)) +
  ggrastr::rasterize(geom_col(aes(fill = gene)), dpi = 600) +
  geom_segment(data = GenePos.tib, aes(x = start, y = ycoord, xend = end, yend = ycoord, color = Names)) +
  geom_text(data = GenePos.tib, aes(x = mid, y = 1.03 * ycoord, label = Names), size = 2) +
  scale_y_continuous("mean coverage") +
  scale_fill_manual(values = col.list, na.value = "grey70") +
  scale_color_manual(values = col.list, na.value = "grey70") +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black")
  )
ggsave("./figure_AML_coevolution/figures/20240528_mtDNA_coverage_AML1.svg", width = 3, height = 2, plot = p)
