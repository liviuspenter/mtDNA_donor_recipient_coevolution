library(ComplexHeatmap)
library(dplyr)
library(Seurat)


so.all <- readRDS("./data/AML_coevolution/objects/20240416_AML_so.rds")

cluster.colors <- as.character(BuenColors::jdb_palette(name = "corona", n = 16))
names(cluster.colors) <- seq(0, 15)

AML1026.1.recipient.df <- data.table::fread(file = "./data/AML_coevolution/AML1026/AML1026.1.recipient.csv", header = T, data.table = F, dec = ",")
rownames(AML1026.1.recipient.df) <- AML1026.1.recipient.df$V1
AML1026.1.recipient.df <- AML1026.1.recipient.df[, -1]

AML1026.3.recipient.df <- data.table::fread(file = "./data/AML_coevolution/AML1026/AML1026.3.recipient.csv", header = T, data.table = F, dec = ",")
rownames(AML1026.3.recipient.df) <- AML1026.3.recipient.df$V1
AML1026.3.recipient.df <- AML1026.3.recipient.df[, -1]

AML1026.1.donor.df <- data.table::fread(file = "./data/AML_coevolution/AML1026/AML1026.1.donor.csv", header = T, data.table = F, dec = ",")
rownames(AML1026.1.donor.df) <- AML1026.1.donor.df$V1
AML1026.1.donor.df <- AML1026.1.donor.df[, -1]

AML1026.3.donor.df <- data.table::fread(file = "./data/AML_coevolution/AML1026/AML1026.3.donor.csv", header = T, data.table = F, dec = ",")
rownames(AML1026.3.donor.df) <- AML1026.3.donor.df$V1
AML1026.3.donor.df <- AML1026.3.donor.df[, -1]

AML1026.leukemia.cells <- intersect(
  colnames(so.all)[which(so.all$manual.cluster %in% c("HSC", "Progenitor", "Mono", "Erythroid"))],
  c(rownames(AML1026.1.recipient.df), rownames(AML1026.3.recipient.df))
)

AML1026.leukemia.cells <- names(sort(so.all$manual.cluster[AML1026.leukemia.cells]))

AML1026.recipient.immune.cells <- intersect(
  colnames(so.all)[which(!so.all$manual.cluster %in% c("HSC", "Progenitor", "Mono", "Erythroid"))],
  c(rownames(AML1026.1.recipient.df), rownames(AML1026.3.recipient.df))
)
AML1026.donor.immune.cells <- intersect(
  colnames(so.all)[which(!so.all$manual.cluster %in% c("HSC", "Progenitor", "Mono", "Erythroid"))],
  c(rownames(AML1026.1.donor.df), rownames(AML1026.3.donor.df))
)


AML1026.heteroplasmy <- data.table::fread("./data/AML_coevolution/AML1026/20240421_AML1026_heteroplasmy_mgatk.csv") %>% as.data.frame()
rownames(AML1026.heteroplasmy) <- AML1026.heteroplasmy$V1
AML1026.heteroplasmy <- AML1026.heteroplasmy[, -1]
AML1026.leukemia.cells <- intersect(AML1026.leukemia.cells, rownames(AML1026.heteroplasmy))

mito.mutations <- intersect(
  colnames(AML1026.1.recipient.df)[which(grepl("chrM", colnames(AML1026.1.recipient.df)))],
  colnames(AML1026.3.recipient.df)[which(grepl("chrM", colnames(AML1026.3.recipient.df)))]
)
mito.mutations <- setdiff(mito.mutations, names(which(colMeans(AML1026.1.recipient.df[, mito.mutations]) > 90)))
mito.mutations <- setdiff(mito.mutations, names(which(colMeans(AML1026.3.recipient.df[, mito.mutations]) > 90)))
mito.mutations <- setdiff(mito.mutations, names(which(colMeans(AML1026.1.recipient.df[, mito.mutations]) < 1)))
mito.mutations <- setdiff(mito.mutations, names(which(colMeans(AML1026.3.recipient.df[, mito.mutations]) < 1)))

auto.mutations <- intersect(
  colnames(AML1026.1.recipient.df)[which(!grepl("chrM", colnames(AML1026.1.recipient.df)))],
  colnames(AML1026.3.recipient.df)[which(!grepl("chrM", colnames(AML1026.3.recipient.df)))]
)
auto.mutations <- setdiff(auto.mutations, names(which(colMeans(AML1026.1.recipient.df[, auto.mutations]) > 90)))
auto.mutations <- setdiff(auto.mutations, names(which(colMeans(AML1026.3.recipient.df[, auto.mutations]) > 90)))
auto.mutations <- setdiff(auto.mutations, names(which(colMeans(AML1026.1.recipient.df[, auto.mutations]) < 1)))
auto.mutations <- setdiff(auto.mutations, names(which(colMeans(AML1026.3.recipient.df[, auto.mutations]) < 1)))

combined.df <- rbind(
  AML1026.1.recipient.df[, auto.mutations], AML1026.3.recipient.df[, auto.mutations],
  AML1026.1.donor.df[, auto.mutations], AML1026.3.donor.df[, auto.mutations]
)
combined.df <- cbind(combined.df, 100 * AML1026.heteroplasmy[rownames(combined.df), ])

citeseq.df <- GetAssayData(so.all, layer = "scale.data")[c("CD34", "CD117", "CD33", "CD14", "CD16", "CD123", "CD11b", "CD11c", "CD71", "CD3", "CD4", "CD8", "CD56", "CD19", "CD22", "CD138"), ]


# find clusters with auto and mito mutations
relevant.auto.mutations <- c(
  "chr3:39119933:C/T", "chr21:44637679:T/A", "chr12:106376431:C/G", "chrX:101041322:G/A", "chr11:5545344:C/G", "chrX:49880642:G/A",
  "chr12:57489913:C/T", "chr6:70294211:C/A",
  "chr15:40294323:G/A", "chr11:64343825:G/C", "chr16:4600149:G/C", "GAK:chr4:896481:G/A", "chr1:237784006:C/G"
)

mito.mutations <- colnames(AML1026.heteroplasmy)

combined.so <- CreateSeuratObject(t(combined.df[AML1026.leukemia.cells, c(relevant.auto.mutations, mito.mutations)]))
combined.so <- FindVariableFeatures(combined.so)
combined.so <- NormalizeData(combined.so)
combined.so <- ScaleData(combined.so, features = rownames(combined.so))
combined.so <- RunPCA(combined.so)
combined.so <- RunUMAP(combined.so, dims = 1:10)
combined.so <- FindNeighbors(combined.so)
combined.so <- FindClusters(combined.so, resolution = 0.1)
saveRDS(file = "./data/AML_coevolution/objects/20240423_AML1026_combined.rds", combined.so)

combined.so <- readRDS(file = "./data/AML_coevolution/objects/20240423_AML1026_combined.rds")
combined.markers <- FindAllMarkers(combined.so, min.pct = 0.05, logfc.threshold = 0.05)
combined.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10

cluster.order <- c(2, 1, 4, 0, 3)
mito.mutation.order <- top10[which(grepl(">", top10$gene)), ] %>% arrange(factor(cluster, levels = cluster.order))
mito.mutation.order <- unique(mito.mutation.order$gene)

cells <- AML1026.leukemia.cells[which(grepl("AML2#", AML1026.leukemia.cells))]
ha <- columnAnnotation(
  cluster = factor(combined.so$seurat_clusters[cells], levels = cluster.order),
  col = list("cluster" = cluster.colors),
  simple_anno_size = unit(5, "pt"),
  border = T, show_legend = F, show_annotation_name = F
)
h0 <- Heatmap(citeseq.df[c("CD34", "CD33", "CD14", "CD123", "CD11b", "CD11c"), cells],
  top_annotation = ha,
  column_split = factor(combined.so$seurat_clusters[cells], levels = cluster.order),
  cluster_rows = F, cluster_columns = F,
  show_column_names = F, show_row_dend = F, row_names_side = "left", row_names_gp = gpar(fontsize = 6),
  col = col_fun_cite,
  border = T, use_raster = T, raster_quality = 10
)

h1 <- Heatmap(t(combined.df[cells, c("chr17:76736877:G/T", "chr1:114716124:C/G", "chr2:197401788:C/T", relevant.auto.mutations)]), # top_annotation = ha,
  column_split = factor(combined.so$seurat_clusters[cells], levels = cluster.order),
  cluster_rows = F, cluster_columns = F,
  show_column_names = F, show_row_dend = F, row_names_side = "left", row_names_gp = gpar(fontsize = 6),
  col = col_fun2, border = T, use_raster = T, raster_quality = 10
)
h2 <- Heatmap(t(combined.df[cells, mito.mutation.order]),
  column_split = factor(combined.so$seurat_clusters[cells], levels = cluster.order),
  cluster_rows = F, cluster_columns = F,
  show_column_names = F, row_names_side = "left", row_names_gp = gpar(fontsize = 6),
  col = col_fun, border = T, use_raster = T, raster_quality = 10
)

svglite::svglite("./figure_AML_coevolution/figures/AML1026/20240423_AML1026.1_mtDNA_clustering.svg", width = 5, height = 5)
h0 %v% h1 %v% h2
dev.off()

write.table(mito.mutation.order, file = "./data/AML_coevolution/AML1026/20240624_AML1026_combined_markers.csv", sep = "\t")


cells <- AML1026.leukemia.cells[which(grepl("AML3#", AML1026.leukemia.cells))]
ha <- columnAnnotation(
  cluster = factor(combined.so$seurat_clusters[cells], levels = cluster.order),
  col = list("cluster" = cluster.colors),
  simple_anno_size = unit(5, "pt"),
  border = T, show_legend = F, show_annotation_name = F
)
h0 <- Heatmap(citeseq.df[c("CD34", "CD117", "CD33", "CD14", "CD123", "CD11b", "CD11c", "CD71"), cells],
  top_annotation = ha,
  column_split = factor(combined.so$seurat_clusters[cells], levels = cluster.order),
  cluster_rows = F, cluster_columns = F,
  show_column_names = F, show_row_dend = F, row_names_side = "left", row_names_gp = gpar(fontsize = 6),
  col = col_fun_cite,
  border = T, use_raster = T, raster_quality = 10
)

h1 <- Heatmap(t(combined.df[cells, c("chr17:76736877:G/T", "chr1:114716124:C/G", "chr2:197401788:C/T", relevant.auto.mutations)]), # top_annotation = ha,
  column_split = factor(combined.so$seurat_clusters[cells], levels = cluster.order),
  cluster_rows = F, cluster_columns = F,
  show_column_names = F, show_row_dend = F, row_names_side = "left", row_names_gp = gpar(fontsize = 6),
  col = col_fun2, border = T, use_raster = T, raster_quality = 10
)
h2 <- Heatmap(t(combined.df[cells, mito.mutation.order]),
  column_split = factor(combined.so$seurat_clusters[cells], levels = cluster.order),
  cluster_rows = F, cluster_columns = F,
  show_column_names = F, row_names_side = "left", row_names_gp = gpar(fontsize = 6),
  col = col_fun, border = T, use_raster = T, raster_quality = 10
)

svglite::svglite("./figure_AML_coevolution/figures/AML1026/20240423_AML1026.3_mtDNA_clustering.svg", width = 5, height = 5)
h0 %v% h1 %v% h2
dev.off()


# changes in subclusters
df <- data.frame(
  cb = AML1026.leukemia.cells,
  cluster = combined.so$seurat_clusters[AML1026.leukemia.cells]
)
df$sample <- ifelse(grepl("AML2#", df$cb), "AML1026.1", "AML1026.3")
df$celltype <- so.all$manual.cluster[df$cb]

stats <- df %>%
  group_by(sample, celltype, cluster) %>%
  tally()
stats$main.cluster <- ifelse(stats$cluster %in% c(0, 1, 3), "m1", "m2")
stats$freq <- apply(stats, MARGIN = 1, FUN = function(x) {
  as.numeric(x["n"]) /
    sum(as.numeric(stats$n[which(stats$sample == x["sample"] & stats$celltype == x["celltype"])]))
})
stats$cluster <- factor(stats$cluster, levels = cluster.order)

p <- ggplot(stats[which(stats$sample == "AML1026.1"), ], aes(x = celltype, y = 100 * freq)) +
  geom_col(aes(group = cluster, fill = cluster)) +
  scale_fill_manual(values = cluster.colors) +
  scale_y_continuous("% cells") +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black"),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    axis.title.x = element_blank()
  )
ggsave("./figure_AML_coevolution/figures/AML1026/20240514_AML1026.1_subclones.svg", width = 1.2, height = 2, plot = p)

p <- ggplot(stats[which(stats$sample == "AML1026.3"), ], aes(x = celltype, y = 100 * freq)) +
  geom_col(aes(group = cluster, fill = cluster)) +
  scale_fill_manual(values = cluster.colors) +
  scale_y_continuous("% cells") +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black"),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    axis.title.x = element_blank()
  )
ggsave("./figure_AML_coevolution/figures/AML1026/20240514_AML1026.3_subclones.svg", width = 1.2, height = 2, plot = p)
