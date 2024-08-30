library(ArchR)
library(gplots)
library(grid)
library(Seurat)

addArchRGenome("hg38")

AML.1007.mito <- loadArchRProject("./data/10026/AML.1007.mito/")
AML.1007.mito$manual.clusters <- c("AML")

AML.1007.mito$manual.clusters[which(AML.1007.mito$Clusters %in% c("C2"))] <- c("HSC")
AML.1007.mito$manual.clusters[which(AML.1007.mito$Clusters %in% c("C1", "C5"))] <- c("Mono")
AML.1007.mito$manual.clusters[which(AML.1007.mito$Clusters %in% c("C3", "C4"))] <- c("GMP")
AML.1007.mito$manual.clusters[which(AML.1007.mito$Clusters %in% c("C6"))] <- c("erythroid")
AML.1007.mito$manual.clusters[which(AML.1007.mito$Clusters %in% "C8")] <- c("CD4")
AML.1007.mito$manual.clusters[which(AML.1007.mito$Clusters %in% "C7")] <- c("CD8")


AML.1007.mito.colors <- c(
  "HSC" = RColorBrewer::brewer.pal(n = 5, name = "Reds")[5],
  "GMP" = RColorBrewer::brewer.pal(n = 5, name = "Reds")[4],
  "Mono" = RColorBrewer::brewer.pal(n = 5, name = "Reds")[3],
  "erythroid" = RColorBrewer::brewer.pal(n = 5, name = "Reds")[2],
  "CD4" = RColorBrewer::brewer.pal(n = 3, name = "Blues")[2],
  "CD8" = RColorBrewer::brewer.pal(n = 3, name = "Blues")[3]
)

p <- plotEmbedding(AML.1007.mito, name = "manual.clusters", pal = c(AML.1007.mito.colors))
plotPDF(ArchRProj = AML.1007.mito, name = "manual.clusters", plotList = list(p), addDOC = F)
saveArchRProject(AML.1007.mito)

# cluster over time
cells.mat <- data.frame(sample = as.character(), cluster = as.character(), cells = as.numeric())
for (s in c("AML1007_1", "AML1007_5", "AML1007_6")) {
  for (cluster in names(AML.1007.mito.colors)) {
    boo <- length(which(AML.1007.mito$Sample == s & AML.1007.mito$manual.clusters == cluster))
    cells.mat <- rbind(cells.mat, data.frame(
      sample = s, cluster = cluster,
      cells = boo,
      freq = boo / length(which(AML.1007.mito$Sample == s))
    ))
  }
}
cells.mat$cluster <- factor(cells.mat$cluster, levels = names(AML.1007.mito.colors))
cells.mat$sample <- factor(cells.mat$sample, levels = c("AML1007_1", "AML1007_5", "AML1007_6"))
ggplot(data = cells.mat, aes(x = sample, y = 100 * freq, group = cluster, fill = cluster)) +
  geom_col() +
  scale_fill_manual(values = AML.1007.mito.colors) +
  scale_x_discrete(labels = c("Screening", "EOT", "Relapse\npost-HSCT2")) +
  scale_y_continuous("%cells") +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.title = element_text("Arial", size = 10, color = "black"),
    axis.text.x = element_text("Arial", size = 10, color = "black", angle = 90, hjust = 1, vjust = 0.5),
    axis.text = element_text("Arial", size = 10, color = "black")
  )
ggsave("./figure_10026/figures/plots/20240523_AML1007_cluster_kinetics.svg", width = 1.2, height = 2)


ggplot(data = cells.mat[which(cells.mat$sample %in% c("AML1007_1", "AML1007_5")), ], aes(x = sample, y = 100 * freq, group = cluster, fill = cluster)) +
  geom_col() +
  scale_fill_manual(values = AML.1007.mito.colors) +
  scale_x_discrete(labels = c("Screening", "EOT")) +
  scale_y_continuous("%cells") +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.title = element_text("Arial", size = 10, color = "black"),
    axis.text.x = element_text("Arial", size = 10, color = "black", angle = 90, hjust = 1, vjust = 0.5),
    axis.text = element_text("Arial", size = 10, color = "black")
  )
ggsave("./figure_10026/figures/plots/20240604_AML1007_cluster_kinetics.svg", width = 1.0, height = 2)
