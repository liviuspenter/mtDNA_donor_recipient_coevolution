library(ComplexHeatmap)
library(dplyr)
library(Seurat)

so.all <- readRDS("./data/AML_coevolution/objects/20240416_AML_so.rds")

cluster.colors <- as.character(BuenColors::jdb_palette(name = "corona", n = 16))
names(cluster.colors) <- seq(0, 15)

col_fun <- circlize::colorRamp2(breaks = seq(0, 10, 10 / 8), colors = BuenColors::jdb_palette(name = "solar_rojos"))
col_fun2 <- circlize::colorRamp2(breaks = seq(0, 100, 100 / 8), colors = BuenColors::jdb_palette(name = "solar_rojos"))

AML1012.1.recipient.df <- data.table::fread(file = "./data/AML_coevolution/AML1012/AML1012.1.recipient.csv", header = T, data.table = F, dec = ",")
rownames(AML1012.1.recipient.df) <- AML1012.1.recipient.df$V1
AML1012.1.recipient.df <- AML1012.1.recipient.df[, -1]

AML1012.4.recipient.df <- data.table::fread(file = "./data/AML_coevolution/AML1012/AML1012.4.recipient.csv", header = T, data.table = F, dec = ",")
rownames(AML1012.4.recipient.df) <- AML1012.4.recipient.df$V1
AML1012.4.recipient.df <- AML1012.4.recipient.df[, -1]

AML1012.1.donor.df <- data.table::fread(file = "./data/AML_coevolution/AML1012/AML1012.1.donor.csv", header = T, data.table = F, dec = ",")
rownames(AML1012.1.donor.df) <- AML1012.1.donor.df$V1
AML1012.1.donor.df <- AML1012.1.donor.df[, -1]

AML1012.4.donor.df <- data.table::fread(file = "./data/AML_coevolution/AML1012/AML1012.4.donor.csv", header = T, data.table = F, dec = ",")
rownames(AML1012.4.donor.df) <- AML1012.4.donor.df$V1
AML1012.4.donor.df <- AML1012.4.donor.df[, -1]

AML1012.leukemia.cells <- intersect(
  colnames(so.all)[which(so.all$manual.cluster %in% c("HSC", "Progenitor", "Mono", "Erythroid"))],
  c(rownames(AML1012.1.recipient.df), rownames(AML1012.4.recipient.df))
)

AML1012.leukemia.cells <- names(sort(so.all$manual.cluster[AML1012.leukemia.cells]))

AML1012.recipient.immune.cells <- intersect(
  colnames(so.all)[which(!so.all$manual.cluster %in% c("HSC", "Progenitor", "Mono", "Erythroid"))],
  c(rownames(AML1012.1.recipient.df), rownames(AML1012.4.recipient.df))
)
AML1012.donor.immune.cells <- intersect(
  colnames(so.all)[which(!so.all$manual.cluster %in% c("HSC", "Progenitor", "Mono", "Erythroid"))],
  c(rownames(AML1012.1.donor.df), rownames(AML1012.4.donor.df))
)


AML1012.heteroplasmy <- data.table::fread("./data/AML_coevolution/AML1012/20240421_AML1012_heteroplasmy_mgatk.csv") %>% as.data.frame()
rownames(AML1012.heteroplasmy) <- AML1012.heteroplasmy$V1
AML1012.heteroplasmy <- AML1012.heteroplasmy[, -1]
AML1012.leukemia.cells <- intersect(AML1012.leukemia.cells, rownames(AML1012.heteroplasmy))

auto.mutations <- intersect(
  colnames(AML1012.1.recipient.df)[which(!grepl("chrM", colnames(AML1012.1.recipient.df)))],
  colnames(AML1012.4.recipient.df)[which(!grepl("chrM", colnames(AML1012.4.recipient.df)))]
)
auto.mutations <- setdiff(auto.mutations, names(which(colMeans(AML1012.1.recipient.df[, auto.mutations]) > 90)))
auto.mutations <- setdiff(auto.mutations, names(which(colMeans(AML1012.4.recipient.df[, auto.mutations]) > 90)))
auto.mutations <- setdiff(auto.mutations, names(which(colMeans(AML1012.1.recipient.df[, auto.mutations]) < 1)))
auto.mutations <- setdiff(auto.mutations, names(which(colMeans(AML1012.4.recipient.df[, auto.mutations]) < 1)))

mito.mutations <- intersect(
  colnames(AML1012.1.recipient.df)[which(grepl("chrM", colnames(AML1012.1.recipient.df)))],
  colnames(AML1012.4.recipient.df)[which(grepl("chrM", colnames(AML1012.4.recipient.df)))]
)
mito.mutations <- setdiff(mito.mutations, names(which(colMeans(AML1012.1.recipient.df[, mito.mutations]) > 90)))
mito.mutations <- setdiff(mito.mutations, names(which(colMeans(AML1012.4.recipient.df[, mito.mutations]) > 90)))
mito.mutations <- setdiff(mito.mutations, names(which(colMeans(AML1012.1.recipient.df[, mito.mutations]) < 1)))
mito.mutations <- setdiff(mito.mutations, names(which(colMeans(AML1012.4.recipient.df[, mito.mutations]) < 1)))

combined.df <- rbind(
  AML1012.1.recipient.df[, auto.mutations], AML1012.4.recipient.df[, auto.mutations],
  AML1012.1.donor.df[, auto.mutations], AML1012.4.donor.df[, auto.mutations]
)
combined.df <- cbind(combined.df, 100 * AML1012.heteroplasmy[rownames(combined.df), ])

citeseq.df <- GetAssayData(so.all, layer = "scale.data")[c("CD34", "CD117", "CD33", "CD14", "CD16", "CD123", "CD11b", "CD11c", "CD71", "CD3", "CD4", "CD8", "CD56", "CD19", "CD22", "CD138"), ]

# find clusters with auto and mito mutations
relevant.auto.mutations <- auto.mutations
relevant.auto.mutations <- c(
  "MAP4:chr3:47916994:A/G", "chr4:152328233:G/A", "chr4:173292151:G/A", "chr7:122702400:G/A", "chr8:109486749:G/A",
  "THEM6:chr8:143818227:A/C", "chr10:43557531:G/T", "chr11:120420764:C/T", "chr12:2606651:A/G", "chr16:22915253:G/A",
  "chr19:49290212:G/T", "chr22:17727240:AACAAC/A"
)
mito.mutations <- colnames(AML1012.heteroplasmy)

combined.so <- CreateSeuratObject(t(combined.df[AML1012.leukemia.cells, c(relevant.auto.mutations, mito.mutations)]))
combined.so <- FindVariableFeatures(combined.so)
combined.so <- NormalizeData(combined.so)
combined.so <- ScaleData(combined.so, features = rownames(combined.so))
combined.so <- RunPCA(combined.so)
combined.so <- RunUMAP(combined.so, dims = 1:10)
combined.so <- FindNeighbors(combined.so)
combined.so <- FindClusters(combined.so, resolution = 0.1)
saveRDS(file = "./data/AML_coevolution/objects/20240423_AML1012_combined.rds", combined.so)

combined.so <- readRDS(file = "./data/AML_coevolution/objects/20240423_AML1012_combined.rds")
combined.markers <- FindAllMarkers(combined.so, min.pct = 0.05, logfc.threshold = 0.05)
combined.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10

cluster.order <- seq(0, max(as.numeric(top10$cluster)))
mito.mutation.order <- top10[which(grepl(">", top10$gene)), ] %>% arrange(factor(cluster, levels = cluster.order))
mito.mutation.order <- unique(mito.mutation.order$gene)
mito.mutation.order <- setdiff(mito.mutation.order, "2622G>C") # exclude mutation that looks like germline

cells <- AML1012.leukemia.cells[which(grepl("AML1#", AML1012.leukemia.cells))]
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

svglite::svglite("./figure_AML_coevolution/figures/AML1012/20240515_AML1012.1_combined_clustering.svg", width = 5, height = 3.2)
h0 %v% h1 %v% h2
dev.off()


cells <- AML1012.leukemia.cells[which(grepl("AML3#", AML1012.leukemia.cells))]
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

svglite::svglite("./figure_AML_coevolution/figures/AML1012/20240515_AML1012.4_combined_clustering.svg", width = 5, height = 3.2)
h0 %v% h1 %v% h2
dev.off()

write.table(mito.mutation.order, file = "./data/AML_coevolution/AML1012/20240624_AML1012_combined_markers.csv", sep = "\t")

# changes in subclusters
df <- data.frame(
  cb = AML1012.leukemia.cells,
  cluster = combined.so$seurat_clusters[AML1012.leukemia.cells]
)
df$sample <- ifelse(grepl("AML1#", df$cb), "AML1012.1", "AML1012.4")
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

p <- ggplot(stats[which(stats$sample == "AML1012.1"), ], aes(x = celltype, y = 100 * freq)) +
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
ggsave("./figure_AML_coevolution/figures/AML1012/20240514_AML1012.1_subclones.svg", width = 1.2, height = 2, plot = p)

p <- ggplot(stats[which(stats$sample == "AML1012.4"), ], aes(x = celltype, y = 100 * freq)) +
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
ggsave("./figure_AML_coevolution/figures/AML1012/20240514_AML1012.4_subclones.svg", width = 1.2, height = 2, plot = p)

# demonstrate dynamics of individual mtDNA and somatic nuclear DNA mutations
cells.1 <- AML1012.leukemia.cells[which(grepl("AML1#", AML1012.leukemia.cells))]
cells.4 <- AML1012.leukemia.cells[which(grepl("AML3#", AML1012.leukemia.cells))]
dynamics.df <- data.frame(
  mutation = c(relevant.auto.mutations, mito.mutation.order),
  cells.1 = sapply(c(relevant.auto.mutations, mito.mutation.order),
    FUN = function(x) {
      length(which(combined.df[cells.1, x] != 0))
    }
  ),
  cells.4 = sapply(c(relevant.auto.mutations, mito.mutation.order),
    FUN = function(x) {
      length(which(combined.df[cells.4, x] != 0))
    }
  )
)
dynamics.df$cells.1.freq <- dynamics.df$cells.1 / length(cells.1)
dynamics.df$cells.4.freq <- dynamics.df$cells.4 / length(cells.4)
dynamics.df$group <- c(
  rep("auto", length(relevant.auto.mutations)),
  rep("mito", length(mito.mutation.order))
)
dynamics.df <- dynamics.df %>% pivot_longer(names_to = "timepoint", cols = c("cells.1.freq", "cells.4.freq"))
dynamics.df$p.value <- apply(dynamics.df, MARGIN = 1, FUN = function(x) {
  fisher.test(matrix(c(
    length(cells.1), length(cells.4),
    as.numeric(x["cells.1"]), as.numeric(x["cells.4"])
  ), nrow = 2))$p.value
})
dynamics.df$p.adj <- p.adjust(dynamics.df$p.value)
dynamics.df$significant <- ifelse(dynamics.df$p.adj < 0.05, "significant", "non.significant")

ggplot() +
  geom_line(data = dynamics.df[which(dynamics.df$group == "auto"), ], aes(x = timepoint, y = 100 * value, color = significant, group = mutation)) +
  geom_point(data = dynamics.df[which(dynamics.df$group == "auto"), ], aes(x = timepoint, y = 100 * value, color = significant), size = 0.5) +
  scale_color_manual(values = c("non.significant" = "grey", "significant" = "firebrick")) +
  scale_x_discrete(labels = c("Screening", "C4")) +
  scale_y_continuous("% somatic mutation", limits = c(0, 100)) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title = element_text("Arial", size = 10, color = "black"),
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    axis.title.x = element_blank()
  )
ggsave("./figure_AML_coevolution/figures/AML1012/20240515_somatic_mutations_dynamics.svg", width = 1.2, height = 2)

ggplot() +
  geom_line(data = dynamics.df[which(dynamics.df$group == "mito"), ], aes(x = timepoint, y = 100 * value, color = significant, group = mutation)) +
  geom_point(data = dynamics.df[which(dynamics.df$group == "mito"), ], aes(x = timepoint, y = 100 * value, color = significant), size = 0.5) +
  scale_color_manual(values = c("non.significant" = "grey", "significant" = "firebrick")) +
  scale_x_discrete(labels = c("Screening", "C4")) +
  scale_y_continuous("% mtDNA mutation", limits = c(0, 100)) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title = element_text("Arial", size = 10, color = "black"),
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    axis.title.x = element_blank()
  )
ggsave("./figure_AML_coevolution/figures/AML1012/20240515_mtDNA_mutations_dynamics.svg", width = 1.2, height = 2)


df <- cbind(df, combined.df[df$cb, ])

stats2 <- df %>%
  group_by(sample, cluster) %>%
  summarize(
    mutated.FBXW7 = length(which(`chr4:152328233:G/A` != 0)),
    mutated.13708 = length(which(`13708G>C` != 0)),
    mutated.14739 = length(which(`14739G>C` != 0)),
    n = length(cluster)
  )
stats2$freq.FBXW7 <- stats2$mutated.FBXW7 / stats2$n
stats2$freq.13708 <- stats2$mutated.13708 / stats2$n
stats2$freq.14739 <- stats2$mutated.14739 / stats2$n
stats2$sample <- factor(stats2$sample, levels = c("AML1012.1", "AML1012.4"))

ggplot(stats2, aes(x = cluster, y = 100 * freq.FBXW7, group = sample)) +
  geom_col(aes(fill = sample), position = "dodge") +
  scale_fill_manual(values = c("AML1012.1" = "grey", "AML1012.4" = "firebrick")) +
  scale_y_continuous("% FBXW7R465C") +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title = element_text("Arial", size = 10, color = "black"),
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.text.y = element_text(angle = 90, vjust = 1, hjust = 0.5)
  )
ggsave("./figure_AML_coevolution/figures/AML1012/20240515_dynamics_FBXW7.svg", width = 1.5, height = 1.5)

ggplot(stats2, aes(x = cluster, y = 100 * freq.13708, group = sample)) +
  geom_col(aes(fill = sample), position = "dodge") +
  scale_fill_manual(values = c("AML1012.1" = "grey", "AML1012.4" = "firebrick")) +
  scale_y_continuous("% 13708G>C") +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title = element_text("Arial", size = 10, color = "black"),
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.text.y = element_text(angle = 90, vjust = 1, hjust = 0.5)
  )
ggsave("./figure_AML_coevolution/figures/AML1012/20240515_dynamics_13708.svg", width = 1.5, height = 1.5)

ggplot(stats2, aes(x = cluster, y = 100 * freq.14739, group = sample)) +
  geom_col(aes(fill = sample), position = "dodge") +
  scale_fill_manual(values = c("AML1012.1" = "grey", "AML1012.4" = "firebrick")) +
  scale_y_continuous("% 14739G>C") +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title = element_text("Arial", size = 10, color = "black"),
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.text.y = element_text(angle = 90, vjust = 1, hjust = 0.5)
  )
ggsave("./figure_AML_coevolution/figures/AML1012/20240515_dynamics_14739.svg", width = 1.5, height = 1.5)
