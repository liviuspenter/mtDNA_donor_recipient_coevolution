library(ArchR)
library(dplyr)
library(gplots)
library(grid)
library(Seurat)

addArchRGenome("hg38")

AML.1011.mito <- loadArchRProject("./data/10026//AML.1011.mito/")
AML.1011.mito$manual.clusters <- c("AML")
AML.1011.mito$manual.clusters[which(AML.1011.mito$Clusters %in% "C5")] <- c("B cell")
AML.1011.mito$manual.clusters[which(AML.1011.mito$Clusters %in% "C6")] <- c("CD4")
AML.1011.mito$manual.clusters[which(AML.1011.mito$Clusters %in% "C7")] <- c("CD8")
AML.1011.mito$manual.clusters[which(AML.1011.mito$Clusters %in% "C8")] <- c("NK")

AML.1011.mito$manual.clusters[which(AML.1011.mito$Clusters %in% c("C3", "C4", "C9"))] <- c("Mono")
AML.1011.mito$manual.clusters[which(AML.1011.mito$Clusters %in% c("C1"))] <- c("erythroid")
AML.1011.mito$manual.clusters[which(AML.1011.mito$Clusters %in% c("C12"))] <- c("HSC")
AML.1011.mito$manual.clusters[which(AML.1011.mito$Clusters %in% c("C10", "C11", "C2"))] <- c("GMP")


AML.1011.mito.colors <- c(
  "HSC" = RColorBrewer::brewer.pal(n = 5, name = "Reds")[5],
  "GMP" = RColorBrewer::brewer.pal(n = 5, name = "Reds")[4],
  "Mono" = RColorBrewer::brewer.pal(n = 5, name = "Reds")[3],
  "erythroid" = RColorBrewer::brewer.pal(n = 5, name = "Reds")[2],
  "CD4" = RColorBrewer::brewer.pal(n = 3, name = "Blues")[2],
  "CD8" = RColorBrewer::brewer.pal(n = 3, name = "Blues")[3],
  "NK" = "purple",
  "B cell" = "orange"
)

p <- plotEmbedding(AML.1011.mito, name = "manual.clusters", pal = c(AML.1011.mito.colors))
plotPDF(ArchRProj = AML.1011.mito, name = "manual.clusters", plotList = list(p), addDOC = F)
saveArchRProject(AML.1011.mito)

donor.variants <- gtools::mixedsort(c(
  "73A>G", "295C>T", "462C>T", "489T>C", "3010G>A", "4216T>C", "10398A>G", "11251A>G", "11719G>A", "12612A>G",
  "13708G>A", "14766C>T", "14798T>C", "15452C>A", "16069C>T", "16126T>C", "16163A>G", "16519T>C"
))
recipient.variants <- gtools::mixedsort(c("72T>C", "195T>C", "4580G>A", "4639T>C", "5263C>T", "7299A>G", "8869A>G", "15904C>T", "16298T>C"))

### identify recipient and donor-derived cells
# HSC, GMP, Mono and erythroid cells are 99.9% recipient
combined.frequencies <- readRDS("./data/10026/mtDNA/20240516_AML1011_combined_mutation_frequencies.rds")
combined.frequencies <- combined.frequencies[-which(rownames(combined.frequencies) == "310T>C"), ]
combined.frequencies <- combined.frequencies[c(recipient.variants, donor.variants), ]

df <- data.frame(
  recipient.variants = colMeans(combined.frequencies[recipient.variants, ]),
  donor.variants = colMeans(combined.frequencies[donor.variants, ]),
  bc = colnames(combined.frequencies)
)
df$individual <- "none"
df$individual[which(df$recipient.variants > 0.8 & df$donor.variants < 0.2)] <- "recipient"
df$individual[which(df$recipient.variants < 0.2 & df$donor.variants > 0.8)] <- "donor"
df <- df[intersect(df$bc, AML.1011.mito$cellNames), ]
df$celltype <- AML.1011.mito[df$bc]$manual.clusters

View(df %>% group_by(celltype, individual) %>% tally())

write.csv2(df, file = "./data/10026/mtDNA/20240522_AML1011_donor_recipient.csv")

# cluster over time
cells.mat <- data.frame(sample = as.character(), cluster = as.character(), cells = as.numeric())
for (s in c("AML1011_1", "AML1011_2", "AML1011_A", "AML1011_B")) {
  for (cluster in names(AML.1011.mito.colors)) {
    boo <- length(which(AML.1011.mito$Sample == s & AML.1011.mito$manual.clusters == cluster))
    cells.mat <- rbind(cells.mat, data.frame(
      sample = s, cluster = cluster,
      cells = boo,
      freq = boo / length(which(AML.1011.mito$Sample == s))
    ))
  }
}
cells.mat$cluster <- factor(cells.mat$cluster, levels = names(AML.1011.mito.colors))
cells.mat$sample <- factor(cells.mat$sample, levels = c("AML1011_A", "AML1011_B", "AML1011_1", "AML1011_2"))
ggplot(data = cells.mat, aes(x = sample, y = 100 * freq, group = cluster, fill = cluster)) +
  geom_col() +
  scale_fill_manual(values = AML.1011.mito.colors) +
  scale_x_discrete(labels = c("Relapse\npre-HSCT", "Remission\npre-HSCT", "Screening", "C2")) +
  scale_y_continuous("%cells") +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.title = element_text("Arial", size = 10, color = "black"),
    axis.text.x = element_text("Arial", size = 10, color = "black", angle = 90, hjust = 1, vjust = 0.5),
    axis.text = element_text("Arial", size = 10, color = "black")
  )
ggsave("./figure_10026/figures/plots/20240523_AML1011_cluster_kinetics.svg", width = 1.5, height = 2)


ggplot(data = cells.mat[which(cells.mat$sample %in% c("AML1011_1", "AML1011_2")), ], aes(x = sample, y = 100 * freq, group = cluster, fill = cluster)) +
  geom_col() +
  scale_fill_manual(values = AML.1011.mito.colors) +
  scale_x_discrete(labels = c("Screening", "C2")) +
  scale_y_continuous("%cells") +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.title = element_text("Arial", size = 10, color = "black"),
    axis.text.x = element_text("Arial", size = 10, color = "black", angle = 90, hjust = 1, vjust = 0.5),
    axis.text = element_text("Arial", size = 10, color = "black")
  )
ggsave("./figure_10026/figures/plots/20240604_AML1011_cluster_kinetics.svg", width = 1.0, height = 2)
