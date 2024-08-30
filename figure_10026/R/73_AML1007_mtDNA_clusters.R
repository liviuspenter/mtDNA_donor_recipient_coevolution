library(ArchR)
library(dplyr)
library(gplots)
library(parallel)
library(Seurat)
library(ComplexHeatmap)

source("./figure_10026/R/Tcell.colors.R")
source("./R/20200328_mtscatac_seq.R")

AML.1007.cluster.colors <- as.vector(BuenColors::jdb_palette(name = "corona", n = 16))
names(AML.1007.cluster.colors) <- seq(0, 15)

AML.1007.mito <- loadArchRProject("./data/10026/AML.1007.mito/")

### 1007
s <- "AML1007"

# read mtDNA data
variants <- read.table(paste0("./data/10026/mtDNA/20240516_", s, "_high_confidence_variants.csv"), header = T)
combined.frequencies <- readRDS(paste0("./data/10026/mtDNA/20240516_", s, "_combined_mutation_frequencies.rds"))
combined.frequencies <- combined.frequencies[-which(rownames(combined.frequencies) == "310T>C"), ]

# cells with mtDNA information
cells <- AML.1007.mito$cellNames[which(grepl(s, AML.1007.mito$Sample) & !AML.1007.mito$manual.clusters %in% c("CD4", "CD8"))]
cells <- cells[which(cells %in% colnames(combined.frequencies))]
combined.frequencies <- combined.frequencies[, cells]

mtDNA.so <- CreateSeuratObject(combined.frequencies)
mtDNA.so <- FindVariableFeatures(mtDNA.so)
mtDNA.so <- NormalizeData(mtDNA.so)
mtDNA.so <- ScaleData(mtDNA.so, features = rownames(mtDNA.so))
mtDNA.so <- RunPCA(mtDNA.so)
mtDNA.so <- RunUMAP(mtDNA.so, dims = 1:10)
mtDNA.so <- FindNeighbors(mtDNA.so)
mtDNA.so <- FindClusters(mtDNA.so, resolution = 0.4)
saveRDS(mtDNA.so, file = "./data/objects/20240523_AML1007_mtDNA.so.rds")

mtDNA.so <- readRDS(file = "./data/10026/objects/20240523_AML1007_mtDNA.so.rds")
mtDNA.markers <- FindAllMarkers(mtDNA.so, min.pct = 0.05, logfc.threshold = 0.05)
mtDNA.markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC) -> top10
DoHeatmap(subset(mtDNA.so, downsample = 100),
  features = top10$gene, slot = "counts", disp.max = 0.1, group.colors =
    AML.1007.cluster.colors, label = F, raster = T
) + # NoLegend() +
  scale_fill_gradientn(colours = c("white", BuenColors::jdb_palette(name = "solar_rojos")), na.value = "black") +
  theme(axis.text = element_text("Arial", size = 8, color = "black"))
ggsave("./figure_10026/figures/heatmaps/20240523_AML1007_mtDNA_clusters.svg", width = 3, height = 2.5)
ggsave("./figure_10026/figures/heatmaps/20240523_AML1007_mtDNA_clusters.png", width = 3, height = 2.5, dpi = 600)

# mtDNA clusters versus celltypes
AML.1007.mito$mito.cluster <- mtDNA.so$seurat_clusters[AML.1007.mito$cellNames]
df <- data.frame(
  bc = AML.1007.mito$cellNames,
  sample = AML.1007.mito$Sample,
  manual.cluster = AML.1007.mito$manual.clusters,
  mito.cluster = AML.1007.mito$mito.cluster
)
df <- df[which(df$bc %in% cells), ]
df$UMAP1 <- getEmbedding(AML.1007.mito)[df$bc, 1]
df$UMAP2 <- getEmbedding(AML.1007.mito)[df$bc, 2]

df2 <- df %>%
  group_by(sample, manual.cluster, mito.cluster) %>%
  tally()
df2 <- df2 %>%
  group_by(sample, manual.cluster) %>%
  mutate(freq = n / sum(n))
df2$manual.cluster <- factor(df2$manual.cluster, levels = c("HSC", "GMP", "Mono", "erythroid", "CD4", "CD8", "NK", "B cell"))
df2$mito.cluster <- factor(df2$mito.cluster, levels = names(AML.1007.cluster.colors))

for (s in c("AML1007_1", "AML1007_5", "AML1007_6")) {
  p <- ggplot(
    df2[which(df2$sample == s & df2$manual.cluster %in% c("HSC", "GMP", "Mono", "erythroid")), ],
    aes(x = manual.cluster, y = 100 * freq, fill = mito.cluster)
  ) +
    geom_col() +
    scale_fill_manual(values = AML.1007.cluster.colors) +
    scale_x_discrete(labels = c("HSC", "GMP", "Mono", "erythroid")) +
    scale_y_continuous("% cells") +
    theme_classic() +
    theme(
      legend.position = "none",
      axis.title = element_text("Arial", size = 10, color = "black"),
      axis.text = element_text("Arial", size = 10, color = "black"),
      axis.title.x = element_blank(),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
    )
  ggsave(paste0("./figure_10026/figures/plots/AML1007/20240523_mito_cluster_dynamics_", s, ".svg"), width = 1.1, height = 2, plot = p)
}

# test for changes in individual mtDNA mutations between pre and post-HSCT
AML.1007.1 <- AML.1007.mito$cellNames[which(AML.1007.mito$Sample == "AML1007_1" & AML.1007.mito$manual.clusters %in% c("HSC", "GMP", "Mono"))]
AML.1007.5 <- AML.1007.mito$cellNames[which(AML.1007.mito$Sample == "AML1007_5" & AML.1007.mito$manual.clusters %in% c("HSC", "GMP", "Mono"))]
AML.1007.6 <- AML.1007.mito$cellNames[which(AML.1007.mito$Sample == "AML1007_6" & AML.1007.mito$manual.clusters %in% c("HSC", "GMP", "Mono"))]

boo <- data.frame(
  AML.1007.1 = rowMeans(combined.frequencies[, AML.1007.1]),
  AML.1007.5 = rowMeans(combined.frequencies[, AML.1007.5]),
  AML.1007.6 = rowMeans(combined.frequencies[, AML.1007.6])
)
boo <- boo[which(boo$HSC.pre > 0.001 | boo$HSC.post > 0.001), ]


df <- getEmbedding(AML.1007.mito)
colnames(df) <- c("UMAP1", "UMAP2")
df <- df[colnames(combined.frequencies), ]
df <- cbind(df, t(combined.frequencies[, rownames(df)]))
df$sample <- stringr::str_split_fixed(rownames(df), pattern = "#", n = 2)[, 1]
df$`5668G>C`[which(df$`5668G>C` > 0.1)] <- 0.1

df <- df[order(df$`5668G>C`), ]

background <- getEmbedding(AML.1007.mito)
colnames(background) <- c("UMAP1", "UMAP2")

ggplot() +
  geom_point(data = background, aes(x = UMAP1, y = UMAP2), color = "lightgrey") +
  geom_point(data = df[which(df$sample == "AML1007_1"), ], aes(x = UMAP1, y = UMAP2, color = `5668G>C`)) +
  scale_color_gradientn(colours = BuenColors::jdb_palette(name = "solar_rojos"), limits = c(0, 0.1)) +
  theme_classic() +
  theme(legend.position = "none") +
  NoAxes()
ggsave("./figure_10026/figures/umaps/AML1007/20240531_AML1007_1_5668G>C.png", width = 4, height = 4, dpi = 600)

ggplot() +
  geom_point(data = background, aes(x = UMAP1, y = UMAP2), color = "lightgrey") +
  geom_point(data = df[which(df$sample == "AML1007_5"), ], aes(x = UMAP1, y = UMAP2, color = `5668G>C`)) +
  scale_color_gradientn(colours = BuenColors::jdb_palette(name = "solar_rojos"), limits = c(0, 0.1)) +
  theme_classic() +
  theme(legend.position = "none") +
  NoAxes()
ggsave("./figure_10026/figures/umaps/AML1007/20240531_AML1007_5_5668G>C.png", width = 4, height = 4, dpi = 600)

ggplot() +
  geom_point(data = background, aes(x = UMAP1, y = UMAP2), color = "lightgrey") +
  geom_point(data = df[which(df$sample == "AML1007_6"), ], aes(x = UMAP1, y = UMAP2, color = `5668G>C`)) +
  scale_color_gradientn(colours = BuenColors::jdb_palette(name = "solar_rojos"), limits = c(0, 0.1)) +
  theme_classic() +
  theme(legend.position = "none") +
  NoAxes()
ggsave("./figure_10026/figures/umaps/AML1007/20240531_AML1007_6_5668G>C.png", width = 4, height = 4, dpi = 600)


boo <- data.frame(bc = AML.1007.mito$cellNames, cluster = AML.1007.mito$manual.clusters, Sample = AML.1007.mito$Sample)
rownames(boo) <- boo$bc
mtDNA.so$manual.cluster <- boo[colnames(mtDNA.so), "cluster"]
mtDNA.so$Sample <- boo[colnames(mtDNA.so), "Sample"]
boo <- as.data.frame(prop.table(table(Idents(mtDNA.so), mtDNA.so$Sample), margin = 2))
boo$Freq[which(boo$Freq == 0)] <- min(boo$Freq[which(boo$Freq != 0)])

ggplot(
  boo[which(boo$Var2 != "none"), ],
  aes(x = Var2, y = 100 * Freq, color = Var1)
) +
  geom_line(aes(group = Var1)) +
  scale_x_discrete(
    labels = c("Screening", "EOT", "2nd HSCT"),
    limits = c("AML1007_1", "AML1007_5", "AML1007_6")
  ) +
  scale_y_log10("% cells") +
  scale_color_manual(values = ) +
  scale_color_manual(values = BuenColors::jdb_palette(name = "corona", n = 13)) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.title = element_text("Arial", size = 10, color = "black"),
    axis.text.x = element_text("Arial", size = 10, color = "black", angle = 90, hjust = 1, vjust = 0.5),
    axis.text = element_text("Arial", size = 10, color = "black")
  )
ggsave("./figure_10026/figures/plots/20240523_AML1007_mtDNA_clusters.svg", width = 1.2, height = 2)
