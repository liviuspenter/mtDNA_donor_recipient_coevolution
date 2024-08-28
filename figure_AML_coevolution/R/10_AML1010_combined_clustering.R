library(ComplexHeatmap)
library(dplyr)
library(ggplot2)
library(ggpointdensity)
library(Seurat)
library(tidyr)

source("./figure_AML_coevolution/R/00_AML_colors.R")

col_fun_cite <- circlize::colorRamp2(breaks = seq(-3, 3, 6 / 8), colors = BuenColors::jdb_palette(name = "brewer_yes"))
col_fun <- circlize::colorRamp2(breaks = seq(0, 10, 10 / 8), colors = BuenColors::jdb_palette(name = "solar_rojos"))
col_fun2 <- circlize::colorRamp2(breaks = seq(0, 100, 100 / 8), colors = BuenColors::jdb_palette(name = "solar_rojos"))

cluster.colors <- as.character(BuenColors::jdb_palette(name = "corona", n = 16))
names(cluster.colors) <- seq(0, 15)

so.all <- readRDS("./data/AML_coevolution/objects/20240416_AML_so.rds")

AML1010.1.recipient.df <- data.table::fread(file = "./data/AML_coevolution/AML1010/AML1010.1.recipient.csv", header = T, data.table = F, dec = ",")
rownames(AML1010.1.recipient.df) <- AML1010.1.recipient.df$V1
AML1010.1.recipient.df <- AML1010.1.recipient.df[, -1]

AML1010.5.recipient.df <- data.table::fread(file = "./data/AML_coevolution/AML1010/AML1010.5.recipient.csv", header = T, data.table = F, dec = ",")
rownames(AML1010.5.recipient.df) <- AML1010.5.recipient.df$V1
AML1010.5.recipient.df <- AML1010.5.recipient.df[, -1]

AML1010.1.donor.df <- data.table::fread(file = "./data/AML_coevolution/AML1010/AML1010.1.donor.csv", header = T, data.table = F, dec = ",")
rownames(AML1010.1.donor.df) <- AML1010.1.donor.df$V1
AML1010.1.donor.df <- AML1010.1.donor.df[, -1]

AML1010.5.donor.df <- data.table::fread(file = "./data/AML_coevolution/AML1010/AML1010.5.donor.csv", header = T, data.table = F, dec = ",")
rownames(AML1010.5.donor.df) <- AML1010.5.donor.df$V1
AML1010.5.donor.df <- AML1010.5.donor.df[, -1]

AML1010.leukemia.cells <- intersect(
  colnames(so.all)[which(so.all$manual.cluster %in% c("HSC", "Progenitor", "Mono", "Erythroid"))],
  c(rownames(AML1010.1.recipient.df), rownames(AML1010.5.recipient.df))
)

AML1010.leukemia.cells <- names(sort(so.all$manual.cluster[AML1010.leukemia.cells]))

AML1010.recipient.immune.cells <- intersect(
  colnames(so.all)[which(!so.all$manual.cluster %in% c("HSC", "Progenitor", "Mono", "Erythroid"))],
  c(rownames(AML1010.1.recipient.df), rownames(AML1010.5.recipient.df))
)
AML1010.donor.immune.cells <- intersect(
  colnames(so.all)[which(!so.all$manual.cluster %in% c("HSC", "Progenitor", "Mono", "Erythroid"))],
  c(rownames(AML1010.1.donor.df), rownames(AML1010.5.donor.df))
)

AML1010.heteroplasmy <- data.table::fread("./data/AML_coevolution/AML1010/20240421_AML1010_heteroplasmy_mgatk.csv") %>% as.data.frame()
rownames(AML1010.heteroplasmy) <- AML1010.heteroplasmy$V1
AML1010.heteroplasmy <- AML1010.heteroplasmy[, -1]

auto.mutations <- intersect(
  colnames(AML1010.1.recipient.df)[which(!grepl("chrM", colnames(AML1010.1.recipient.df)))],
  colnames(AML1010.5.recipient.df)[which(!grepl("chrM", colnames(AML1010.5.recipient.df)))]
)
auto.mutations <- setdiff(auto.mutations, names(which(colMeans(AML1010.1.recipient.df[, auto.mutations]) > 90)))
auto.mutations <- setdiff(auto.mutations, names(which(colMeans(AML1010.5.recipient.df[, auto.mutations]) > 90)))
auto.mutations <- setdiff(auto.mutations, names(which(colMeans(AML1010.1.recipient.df[, auto.mutations]) < 1)))
auto.mutations <- setdiff(auto.mutations, names(which(colMeans(AML1010.5.recipient.df[, auto.mutations]) < 1)))

mito.mutations <- colnames(AML1010.heteroplasmy)

combined.df <- rbind(
  AML1010.1.recipient.df[, auto.mutations], AML1010.5.recipient.df[, auto.mutations],
  AML1010.1.donor.df[, auto.mutations], AML1010.5.donor.df[, auto.mutations]
)
combined.df <- cbind(combined.df, 100 * AML1010.heteroplasmy[rownames(combined.df), ])

citeseq.df <- GetAssayData(so.all, layer = "scale.data")[c("CD34", "CD117", "CD33", "CD14", "CD11b", "CD11c", "CD71"), ]

# find clusters with mito mutations
combined.so <- CreateSeuratObject(t(combined.df[AML1010.leukemia.cells, ]))
combined.so <- FindVariableFeatures(combined.so)
combined.so <- NormalizeData(combined.so)
combined.so <- ScaleData(combined.so, features = rownames(combined.so))
combined.so <- RunPCA(combined.so)
combined.so <- RunUMAP(combined.so, dims = 1:10)
combined.so <- FindNeighbors(combined.so)
combined.so <- FindClusters(combined.so, resolution = 0.1)

combined.markers <- FindAllMarkers(combined.so, min.pct = 0.05, logfc.threshold = 0.05)
combined.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10

cluster.order <- seq(0, 8)
mito.mutation.order <- top10[which(grepl(">", top10$gene)), ] %>% arrange(factor(cluster, levels = cluster.order))
mito.mutation.order <- unique(mito.mutation.order$gene)

cells <- AML1010.leukemia.cells[which(grepl("AML2#", AML1010.leukemia.cells))]
ha <- columnAnnotation(
  cluster = factor(combined.so$seurat_clusters[cells], levels = cluster.order),
  col = list("cluster" = cluster.colors),
  simple_anno_size = unit(5, "pt"),
  border = T, show_legend = F, show_annotation_name = F
)
h0 <- Heatmap(citeseq.df[, cells],
  top_annotation = ha,
  column_split = factor(combined.so$seurat_clusters[cells], levels = cluster.order),
  cluster_rows = F, cluster_columns = F,
  show_column_names = F, show_row_dend = F, row_names_side = "left", row_names_gp = gpar(fontsize = 6),
  col = col_fun_cite,
  border = T, use_raster = T, raster_quality = 10
)
h1 <- Heatmap(t(combined.df[cells, auto.mutations]),
  top_annotation = ha,
  column_split = factor(combined.so$seurat_clusters[cells], levels = cluster.order),
  cluster_rows = T, cluster_columns = F,
  show_column_names = F, show_row_dend = F, row_names_side = "left", row_names_gp = gpar(fontsize = 6),
  col = col_fun2, border = T, use_raster = T, raster_quality = 10
)
h2 <- Heatmap(t(combined.df[cells, mito.mutation.order]),
  column_split = factor(combined.so$seurat_clusters[cells], levels = cluster.order),
  cluster_rows = F, cluster_columns = F,
  show_column_names = F, row_names_side = "left", row_names_gp = gpar(fontsize = 6),
  col = col_fun, border = T, use_raster = T, raster_quality = 10
)
svglite::svglite("./figure_AML_coevolution/figures/AML1010/20240508_AML1010.1_combined_clustering_raw.svg", width = 6, height = 10)
h0 %v% h1 %v% h2
dev.off()

cluster.order <- seq(0, 8)
mito.mutation.order <- top10[which(grepl(">", top10$gene)), ] %>% arrange(factor(cluster, levels = cluster.order))
mito.mutation.order <- unique(mito.mutation.order$gene)

cells <- AML1010.leukemia.cells[which(grepl("AML1#", AML1010.leukemia.cells))]
ha <- columnAnnotation(
  cluster = factor(combined.so$seurat_clusters[cells], levels = cluster.order),
  col = list("cluster" = cluster.colors),
  simple_anno_size = unit(5, "pt"),
  border = T, show_legend = F, show_annotation_name = F
)
h0 <- Heatmap(citeseq.df[, cells],
  top_annotation = ha,
  column_split = factor(combined.so$seurat_clusters[cells], levels = cluster.order),
  cluster_rows = F, cluster_columns = F,
  show_column_names = F, show_row_dend = F, row_names_side = "left", row_names_gp = gpar(fontsize = 6),
  col = col_fun_cite,
  border = T, use_raster = T, raster_quality = 10
)
h1 <- Heatmap(t(combined.df[cells, auto.mutations]),
  top_annotation = ha,
  column_split = factor(combined.so$seurat_clusters[cells], levels = cluster.order),
  cluster_rows = T, cluster_columns = F,
  show_column_names = F, show_row_dend = F, row_names_side = "left", row_names_gp = gpar(fontsize = 6),
  col = col_fun2, border = T, use_raster = T, raster_quality = 10
)
h2 <- Heatmap(t(combined.df[cells, mito.mutation.order]),
  column_split = factor(combined.so$seurat_clusters[cells], levels = cluster.order),
  cluster_rows = F, cluster_columns = F,
  show_column_names = F, row_names_side = "left", row_names_gp = gpar(fontsize = 6),
  col = col_fun, border = T, use_raster = T, raster_quality = 10
)
svglite::svglite("./figure_AML_coevolution/figures/AML1010/20240508_AML1010.5_combined_clustering_raw.svg", width = 6, height = 10)
h0 %v% h1 %v% h2
dev.off()


# find clusters with relevant auto+mito mutations

# generate curated heatmaps of clustering
relevant.auto.mutations <- c(
  "chr1:114716126:C/T", "chr15:57243509:C/CA", "chr14:23080323:G/A", "chr3:46021532:C/T",
  "chr5:157285369:C/T", "chr10:86900110:A/C", "chr16:82068237:C/T", "chr19:30443831:G/A",
  "DNAH11:chr7:21739637:G/T", "chr6:168653159:G/A", "chr22:50266979:G/C",
  "chr15:57063758:G/GA", "chr19:9186097:T/A",
  "chr4:129102637:G/A", "chr15:40903001:T/C", "SIPA1L3:chr19:38473571:C/T", "chr17:1659937:C/T"
)

combined.so <- CreateSeuratObject(t(combined.df[AML1010.leukemia.cells, c(mito.mutations, relevant.auto.mutations)]))
combined.so <- FindVariableFeatures(combined.so)
combined.so <- NormalizeData(combined.so)
combined.so <- ScaleData(combined.so, features = rownames(combined.so))
combined.so <- RunPCA(combined.so)
combined.so <- RunUMAP(combined.so, dims = 1:10)
combined.so <- FindNeighbors(combined.so)
combined.so <- FindClusters(combined.so, resolution = 0.1)
saveRDS(file = "./data/AML_coevolution/objects/20240508_AML1010_combined.rds", combined.so)
combined.so <- readRDS(file = "./data/AML_coevolution/objects/20240508_AML1010_combined.rds")

combined.markers <- FindAllMarkers(combined.so, min.pct = 0.05, logfc.threshold = 0.05)
combined.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10

# order of clusters
cluster.order <- c(0, 3, 1, 6, 2, 4, 5)
mito.mutation.order <- top10[which(grepl(">", top10$gene)), ] %>% arrange(factor(cluster, levels = cluster.order))
mito.mutation.order <- unique(mito.mutation.order$gene)

# only consider mtDNA mutations that are not homoplasmic
mito.mutation.order <- names(which(colMeans(AML1010.heteroplasmy[AML1010.leukemia.cells, mito.mutation.order]) < 0.95))

cells <- AML1010.leukemia.cells[which(grepl("AML2#", AML1010.leukemia.cells))]
ha <- columnAnnotation(
  cluster = factor(combined.so$seurat_clusters[cells], levels = cluster.order),
  col = list("cluster" = cluster.colors),
  simple_anno_size = unit(5, "pt"),
  border = T, show_legend = F, show_annotation_name = F
)
h0 <- Heatmap(citeseq.df[, cells],
  top_annotation = ha,
  column_split = factor(combined.so$seurat_clusters[cells], levels = cluster.order),
  cluster_rows = F, cluster_columns = F,
  show_column_names = F, show_row_dend = F, row_names_side = "left", row_names_gp = gpar(fontsize = 6),
  col = col_fun_cite,
  border = T, use_raster = T, raster_quality = 10
)
h1 <- Heatmap(t(combined.df[cells, relevant.auto.mutations]), # top_annotation = ha,
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
svglite::svglite("./figure_AML_coevolution/figures/AML1010/20240508_AML1010.1_combined_clustering_curated.svg", width = 6, height = 5)
h0 %v% h1 %v% h2
dev.off()

cells <- AML1010.leukemia.cells[which(grepl("AML1#", AML1010.leukemia.cells))]
ha <- columnAnnotation(
  cluster = factor(combined.so$seurat_clusters[cells], levels = cluster.order),
  col = list("cluster" = cluster.colors),
  simple_anno_size = unit(5, "pt"),
  border = T, show_legend = F, show_annotation_name = F
)
h0 <- Heatmap(citeseq.df[, cells],
  top_annotation = ha,
  column_split = factor(combined.so$seurat_clusters[cells], levels = cluster.order),
  cluster_rows = F, cluster_columns = F,
  show_column_names = F, show_row_dend = F, row_names_side = "left", row_names_gp = gpar(fontsize = 6),
  col = col_fun_cite,
  border = T, use_raster = T, raster_quality = 10
)
h1 <- Heatmap(t(combined.df[cells, relevant.auto.mutations]), # top_annotation = ha,
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
svglite::svglite("./figure_AML_coevolution/figures/AML1010/20240508_AML1010.5_combined_clustering_curated.svg", width = 6, height = 5)
h0 %v% h1 %v% h2
dev.off()

write.table(mito.mutation.order, file = "./data/AML_coevolution/AML1010/20240624_AML1010_combined_markers.csv", sep = "\t")

df <- data.frame(
  cb = AML1010.leukemia.cells,
  cluster = combined.so$seurat_clusters[AML1010.leukemia.cells]
)
df$sample <- ifelse(grepl("AML2#", df$cb), "AML1010.1", "AML1010.5")
df$celltype <- so.all$manual.cluster[df$cb]

stats <- df %>%
  group_by(sample, celltype, cluster) %>%
  tally()
stats$main.cluster <- ifelse(stats$cluster %in% c(0, 1, 3), "m1", "m2")
stats$freq <- apply(stats, MARGIN = 1, FUN = function(x) {
  as.numeric(x["n"]) /
    sum(as.numeric(stats$n[which(stats$sample == x["sample"] & stats$celltype == x["celltype"])]))
})
# stats$cluster = factor(stats$cluster, levels = names(sort(table(df[which(df$sample == 'AML1010.1'), 'cluster']), decreasing = T)))
stats$cluster <- factor(stats$cluster, levels = cluster.order)

p <- ggplot(stats[which(stats$sample == "AML1010.1"), ], aes(x = celltype, y = 100 * freq)) +
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
ggsave("./figure_AML_coevolution/figures/AML1010/20240509_AML1010.1_subclones.svg", width = 1.2, height = 2, plot = p)


p <- ggplot(stats[which(stats$sample == "AML1010.5"), ], aes(x = celltype, y = 100 * freq)) +
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
ggsave("./figure_AML_coevolution/figures/AML1010/20240509_AML1010.5_subclones.svg", width = 1.2, height = 2, plot = p)

# demonstrate NRARwt clone further
NRASwt.subclone <- colnames(combined.so)[which(combined.so$seurat_clusters == "5")]
ggplot() +
  ggrastr::rasterize(geom_point(data = combined.df[AML1010.leukemia.cells, ], aes(x = `chr15:57063758:G/GA`, y = `9820G>C`), color = "grey", size = 0.5), dpi = 600) +
  geom_point(data = combined.df[NRASwt.subclone, ], aes(x = `chr1:114716126:C/T`, y = `9820G>C`), color = cluster.colors[6], size = 0.5) +
  scale_x_continuous("% NRASG12D") +
  scale_y_continuous("% 9820G>C") +
  theme_classic() +
  theme(
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black")
  )
ggsave("./figure_AML_coevolution/figures/AML1010/20240510_NRASwt_clone_NRAS_9820.svg", width = 1.5, height = 1.5)

ggplot() +
  ggrastr::rasterize(geom_point(data = combined.df[AML1010.leukemia.cells, ], aes(x = `chr15:57063758:G/GA`, y = `9820G>C`), color = "grey", size = 0.5), dpi = 600) +
  geom_point(data = combined.df[NRASwt.subclone, ], aes(x = `chr15:57243509:C/CA`, y = `9820G>C`), color = cluster.colors[6], size = 0.5) +
  theme_classic() +
  scale_x_continuous("% TCF12S188fs") +
  scale_y_continuous("% 9820G>C") +
  theme_classic() +
  theme(
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black")
  )
ggsave("./figure_AML_coevolution/figures/AML1010/20240510_NRASwt_clone_TCF12S188fs_9820.svg", width = 1.5, height = 1.5)


ggplot() +
  ggrastr::rasterize(geom_point(data = combined.df[AML1010.leukemia.cells, ], aes(x = `chr15:57063758:G/GA`, y = `9820G>C`), color = "grey", size = 0.5), dpi = 600) +
  geom_point(data = combined.df[NRASwt.subclone, ], aes(x = `chr15:57063758:G/GA`, y = `9820G>C`), color = cluster.colors[6], size = 0.5) +
  theme_classic() +
  scale_x_continuous("% TCF12E65fs") +
  scale_y_continuous("% 9820G>C") +
  theme_classic() +
  theme(
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black")
  )
ggsave("./figure_AML_coevolution/figures/AML1010/20240510_NRASwt_clone_TCF12E65fs_9820.svg", width = 1.5, height = 1.5)







df.embedding <- Embeddings(so.all, reduction = "umap") %>% as.data.frame()
df.embedding$bc <- rownames(df.embedding)
df.embedding$cluster <- combined.so$seurat_clusters[df.embedding$bc]
df.embedding$sample <- ifelse(grepl("AML2#", df.embedding$bc), "AML1010.1", "AML1010.5")

cells.main.cluster1 <- df$cluster %in% c(0, 1, 3)

p <- ggplot(t(citeseq.df[, df$cb[df$cluster %in% c(0, 1, 3)]]) %>% as.data.frame(), aes(x = CD34, y = CD14)) +
  ggrastr::rasterize(geom_pointdensity(size = 0.5), dpi = 600) +
  geom_density_2d(color = "black") +
  scale_color_gradientn(colours = BuenColors::jdb_palette(name = "brewer_celsius")) +
  scale_x_continuous(limits = c(-1.5, 2.5)) +
  scale_y_continuous(limits = c(-1, 5)) +
  theme_classic() +
  theme(
    legend.position = "bottom",
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black")
  )
ggsave("./figure_AML_coevolution/figures/AML1010/20240509_CD34_CD14_main_cluster1.svg", width = 1.5, height = 2.5, plot = p)
p <- ggplot(t(citeseq.df[, df$cb[!df$cluster %in% c(0, 1, 3)]]) %>% as.data.frame(), aes(x = CD34, y = CD14)) +
  ggrastr::rasterize(geom_pointdensity(size = 0.5), dpi = 600) +
  geom_density_2d(color = "black") +
  scale_color_gradientn(colours = BuenColors::jdb_palette(name = "brewer_celsius")) +
  scale_x_continuous(limits = c(-1.5, 2.5)) +
  scale_y_continuous(limits = c(-1, 5)) +
  theme_classic() +
  theme(
    legend.position = "bottom",
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black")
  )
ggsave("./figure_AML_coevolution/figures/AML1010/20240509_CD34_CD14_main_cluster2.svg", width = 1.5, height = 2.5, plot = p)
