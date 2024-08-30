library(ArchR)
library(gplots)
library(grid)
library(Seurat)

addArchRGenome("hg38")

AML.1012.mito <- loadArchRProject("./data/10026/AML.1012.all.mito/")
AML.1012.mito$manual.clusters <- c("AML")
AML.1012.mito$manual.clusters[which(AML.1012.mito$Clusters %in% c("C11", "C12"))] <- c("B cell")
AML.1012.mito$manual.clusters[which(AML.1012.mito$Clusters %in% "C10")] <- c("CD4")
AML.1012.mito$manual.clusters[which(AML.1012.mito$Clusters %in% "C9")] <- c("CD8")
AML.1012.mito$manual.clusters[which(AML.1012.mito$Clusters %in% "C8")] <- c("NK")

AML.1012.mito$manual.clusters[which(AML.1012.mito$Clusters %in% c("C1"))] <- c("erythroid")
AML.1012.mito$manual.clusters[which(AML.1012.mito$Clusters %in% c("C3"))] <- c("HSC")
AML.1012.mito$manual.clusters[which(AML.1012.mito$Clusters %in% c("C2", "C4"))] <- c("GMP")
AML.1012.mito$manual.clusters[which(AML.1012.mito$Clusters %in% c("C5", "C6", "C7"))] <- c("Mono")

AML.1012.mito.colors <- c(
  "HSC" = RColorBrewer::brewer.pal(n = 5, name = "Reds")[5],
  "GMP" = RColorBrewer::brewer.pal(n = 5, name = "Reds")[4],
  "Mono" = RColorBrewer::brewer.pal(n = 5, name = "Reds")[3],
  "erythroid" = RColorBrewer::brewer.pal(n = 5, name = "Reds")[2],
  "CD4" = RColorBrewer::brewer.pal(n = 3, name = "Blues")[2],
  "CD8" = RColorBrewer::brewer.pal(n = 3, name = "Blues")[3],
  "NK" = "purple",
  "B cell" = "orange"
)

p <- plotEmbedding(AML.1012.mito, name = "manual.clusters", pal = c(AML.1012.mito.colors))
plotPDF(ArchRProj = AML.1012.mito, name = "manual.clusters", plotList = list(p), addDOC = F)
saveArchRProject(AML.1012.mito)

# cluster over time
cells.mat <- data.frame(sample = as.character(), cluster = as.character(), cells = as.numeric())
for (s in c("AML1012_1", "AML1012_2", "AML1012_3", "AML1012_4", "AML1012_A")) {
  for (cluster in names(AML.1012.mito.colors)) {
    boo <- length(which(AML.1012.mito$Sample.hash == s & AML.1012.mito$manual.clusters == cluster))
    cells.mat <- rbind(cells.mat, data.frame(
      sample = s, cluster = cluster,
      cells = boo,
      freq = boo / length(which(AML.1012.mito$Sample.hash == s))
    ))
  }
}
cells.mat$cluster <- factor(cells.mat$cluster, levels = names(AML.1012.mito.colors))
cells.mat$sample <- factor(cells.mat$sample, levels = c("AML1012_A", "AML1012_1", "AML1012_2", "AML1013_3", "AML1012_4"))
ggplot(data = cells.mat[which(cells.mat$sample != "AML1012_3"), ], aes(x = sample, y = 100 * freq, group = cluster, fill = cluster)) +
  geom_col() +
  scale_fill_manual(values = AML.1012.mito.colors) +
  scale_x_discrete(labels = c("Relapse pre-HSCT", "Screening", "Decitabine", "EOT")) +
  scale_y_continuous("%cells") +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.title = element_text("Arial", size = 10, color = "black"),
    axis.text.x = element_text("Arial", size = 10, color = "black", angle = 90, hjust = 1, vjust = 0.5),
    axis.text = element_text("Arial", size = 10, color = "black")
  )
ggsave("./figure_10026/figures/plots/20240523_AML1012_cluster_kinetics.svg", width = 1.5, height = 2)
