library(ArchR)
library(dplyr)
library(gplots)
library(parallel)
library(Seurat)
library(ComplexHeatmap)

source("./figure_10026/R/Tcell.colors.R")
source("./R/20200328_mtscatac_seq.R")

AML.1011.cluster.colors <- as.vector(BuenColors::jdb_palette(name = "corona", n = 16))
names(AML.1011.cluster.colors) <- seq(0, 15)

AML.1011.mito <- loadArchRProject("./data/10026/AML.1011.mito/")

donor.variants <- gtools::mixedsort(c(
  "73A>G", "295C>T", "462C>T", "489T>C", "3010G>A", "4216T>C", "10398A>G", "11251A>G", "11719G>A", "12612A>G",
  "13708G>A", "14766C>T", "14798T>C", "15452C>A", "16069C>T", "16126T>C", "16163A>G", "16519T>C"
))
recipient.variants <- gtools::mixedsort(c("72T>C", "195T>C", "4580G>A", "4639T>C", "5263C>T", "7299A>G", "8869A>G", "15904C>T", "16298T>C"))


### 1011
s <- "AML1011"

# read mtDNA data
variants <- read.table(paste0("./data/10026/mtDNA/20240516_", s, "_high_confidence_variants.csv"), header = T)
combined.frequencies <- readRDS(paste0("./data/10026/mtDNA/20240516_", s, "_combined_mutation_frequencies.rds"))
combined.frequencies <- combined.frequencies[-which(rownames(combined.frequencies) == "310T>C"), ]
combined.frequencies <- combined.frequencies[-which(rownames(combined.frequencies) %in% c(donor.variants, recipient.variants)), ]

# cells with mtDNA information
chimerism.df <- read.csv2("./data/10026/mtDNA/20240522_AML1011_donor_recipient.csv", row.names = 1)
cells <- AML.1011.mito$cellNames[which(grepl(s, AML.1011.mito$Sample) & !AML.1011.mito$manual.clusters %in% c("CD4", "CD8", "NK", "B cell"))]
cells <- intersect(cells, chimerism.df$bc[which(chimerism.df$individual == "recipient")])
cells <- cells[which(cells %in% colnames(combined.frequencies))]
combined.frequencies <- combined.frequencies[, cells]

recipient.non.T.cells <- AML.1011.mito$cellNames[which(AML.1011.mito$manual.clusters %in% c("HSCT", "GMP", "Mono", "erythroid"))]

mtDNA.so <- CreateSeuratObject(combined.frequencies)
mtDNA.so <- FindVariableFeatures(mtDNA.so)
mtDNA.so <- NormalizeData(mtDNA.so)
mtDNA.so <- ScaleData(mtDNA.so, features = rownames(mtDNA.so))
mtDNA.so <- RunPCA(mtDNA.so)
mtDNA.so <- RunUMAP(mtDNA.so, dims = 1:10)
mtDNA.so <- FindNeighbors(mtDNA.so)
mtDNA.so <- FindClusters(mtDNA.so, resolution = 0.04)
saveRDS(mtDNA.so, file = "./data/objects/20240523_AML1011_mtDNA.so.rds")

mtDNA.so <- readRDS(file = "./data/10026/objects/20240523_AML1011_mtDNA.so.rds")
mtDNA.markers <- FindAllMarkers(mtDNA.so, min.pct = 0.05, logfc.threshold = 0.05)
mtDNA.markers %>%
  group_by(cluster) %>%
  top_n(n = 4, wt = avg_log2FC) -> top10
DoHeatmap(subset(mtDNA.so, downsample = 100),
  features = top10$gene, slot = "counts", disp.max = 0.1, group.colors =
    AML.1011.cluster.colors, label = F, raster = T
) + # NoLegend() +
  scale_fill_gradientn(colours = c("white", BuenColors::jdb_palette(name = "solar_rojos")), na.value = "black") +
  theme(axis.text = element_text("Arial", size = 8, color = "black"))
ggsave("./figure_10026/figures/heatmaps/20240523_AML1011_mtDNA_clusters.svg", width = 3, height = 2.5)
ggsave("./figure_10026/figures/heatmaps/20240523_AML1011_mtDNA_clusters.png", width = 3, height = 2.5, dpi = 600)

# mtDNA clusters versus celltypes
AML.1011.mito$mito.cluster <- mtDNA.so$seurat_clusters[AML.1011.mito$cellNames]
df <- data.frame(
  bc = AML.1011.mito$cellNames,
  sample = AML.1011.mito$Sample,
  manual.cluster = AML.1011.mito$manual.clusters,
  mito.cluster = AML.1011.mito$mito.cluster
)
df <- df[which(df$bc %in% cells), ]
df$UMAP1 <- getEmbedding(AML.1011.mito)[df$bc, 1]
df$UMAP2 <- getEmbedding(AML.1011.mito)[df$bc, 2]

df2 <- df %>%
  group_by(sample, manual.cluster, mito.cluster) %>%
  tally()
df2 <- df2 %>%
  group_by(sample, manual.cluster) %>%
  mutate(freq = n / sum(n))
df2$manual.cluster <- factor(df2$manual.cluster, levels = c("HSC", "GMP", "Mono", "erythroid", "CD4", "CD8", "NK", "B cell"))
df2$mito.cluster <- factor(df2$mito.cluster, levels = names(AML.1011.cluster.colors))

for (s in c("AML1011_A", "AML1011_B", "AML1011_1", "AML1011_2")) {
  p <- ggplot(
    df2[which(df2$sample == s & df2$manual.cluster %in% c("GMP", "Mono", "erythroid")), ],
    aes(x = manual.cluster, y = 100 * freq, fill = mito.cluster)
  ) +
    geom_col() +
    scale_fill_manual(values = AML.1011.cluster.colors) +
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
  ggsave(paste0("./figure_10026/figures/plots/AML1011/20240523_mito_cluster_dynamics_", s, ".svg"), width = 1.1, height = 2, plot = p)
}

for (s in c("AML1011_1", "AML1011_2", "AML1011_A", "AML1011_B")) {
  p <- ggplot(df[which(df$sample == s), ], aes(x = UMAP1, y = UMAP2, color = mito.cluster)) +
    geom_point(size = 0.5) +
    scale_color_manual(values = AML.1011.cluster.colors) +
    theme_classic() +
    NoAxes() +
    NoLegend()
  ggsave(paste0("./figure_10026/figures/umaps/AML1011/20240523_", s, "_UMAP_mitoclusters.png"), width = 3, height = 3, dpi = 600, plot = p)
}

# small plot for clinical annotation
boo <- data.frame(
  label = c("Relapse pre-HSCT", "Remission post-HSCT", "Screening", "C2"),
  blasts = c(21, 1, 28, 12)
)
boo$label <- factor(boo$label, levels = boo$label)
ggplot(boo, aes(x = label, y = blasts)) +
  geom_line(group = 1) +
  geom_point(size = 0.5) +
  scale_y_continuous("% blasts", limits = c(0, 30)) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title = element_text("Arial", size = 10, color = "black"),
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )
ggsave("./figure_10026/figures/plots/AML1011/20240524_clinical_course.svg", width = 1.3, height = 3)

boo <- data.frame(bc = AML.1011.mito$cellNames, cluster = AML.1011.mito$manual.clusters, Sample = AML.1011.mito$Sample)
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
    labels = c("Relapse pre-HSCT", "Remission pre-HSCT", "Screening", "C2"),
    limits = c("AML1011_A", "AML1011_B", "AML1011_1", "AML1011_2")
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
ggsave("./figure_10026/figures/plots/20240523_AML1011_mtDNA_clusters.svg", width = 1.5, height = 2.7)


saveArchRProject(AML.1011.mito)


df <- getEmbedding(AML.1011.mito)
colnames(df) <- c("UMAP1", "UMAP2")
df$origin <- chimerism.df[rownames(df), "individual"]
df[AML.1011.mito$cellNames, "celltype"] <- AML.1011.mito$manual.clusters
df[AML.1011.mito$cellNames, "Sample"] <- AML.1011.mito$Sample
df <- df[which(df$Sample != "AML1011_A"), ]

ggplot() +
  geom_point(data = df[which(df$origin == "recipient"), ], aes(x = UMAP1, y = UMAP2), color = "lightgrey", size = 0.5) +
  geom_point(data = df[which(df$origin == "donor" & df$celltype %in% c("CD4", "CD8", "NK", "B cell")), ], aes(x = UMAP1, y = UMAP2), color = "lightblue", size = 1) +
  geom_point(data = df[which(df$origin == "donor" & (!df$celltype %in% c("CD4", "CD8", "NK", "B cell"))), ], aes(x = UMAP1, y = UMAP2), color = "black", size = 1) +
  theme_classic() +
  theme(
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank()
  )
ggsave("./figure_10026/figures/umaps/20240531_AML1011_donor_derived_cells.png", width = 2.5, height = 2.5)
