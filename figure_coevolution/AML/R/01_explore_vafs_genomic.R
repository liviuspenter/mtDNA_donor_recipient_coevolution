# explore somatic mutations and germline single nucleotide polymorphisms

library(ComplexHeatmap)
library(dplyr)
library(ggplot2)
library(ggExtra)

# read protein data Seurat object
so.13 <- readRDS("./data/coevolution/objects/20220112_AML1026.rds")

# read VAFs of somatic nuclear mutations
vafs.1 <- as.data.frame(data.table::fread("./data/coevolution/1026_1/AF.csv"))
vafs.1$Barcode <- paste0("1026.1#", vafs.1$Barcode)
rownames(vafs.1) <- vafs.1$Barcode
vafs.1 <- vafs.1[, -c(1, 2)]
vafs.3 <- as.data.frame(data.table::fread("./data/coevolution/1026_3/AF.csv"))
vafs.3$Barcode <- paste0("1026.3#", vafs.3$Barcode)
rownames(vafs.3) <- vafs.3$Barcode
vafs.3 <- vafs.3[, -c(1, 2)]

# read read depths of somatic nuclear mutations
depths.1 <- as.data.frame(data.table::fread("./data/coevolution/1026_1/DP.csv"))
depths.1$Barcode <- paste0("1026.1#", depths.1$Barcode)
rownames(depths.1) <- depths.1$Barcode
depths.1 <- depths.1[, -c(1, 2)]
depths.3 <- as.data.frame(data.table::fread("./data/coevolution/1026_3/DP.csv"))
depths.3$Barcode <- paste0("1026.3#", depths.3$Barcode)
rownames(depths.3) <- depths.3$Barcode
depths.3 <- depths.3[, -c(1, 2)]

# consider mutations that are identified for both samples
vafs <- rbind(vafs.1[, intersect(colnames(vafs.1), colnames(vafs.3))], vafs.3[, intersect(colnames(vafs.1), colnames(vafs.3))])
depths <- rbind(depths.1[, intersect(colnames(depths.1), colnames(depths.3))], depths.3[, intersect(colnames(depths.1), colnames(depths.3))])
# so.13 = AddMetaData(so.13, metadata = vafs)

###
# perform donor-recipient deconvolution using germline single nucleotide polymorphisms
###

# variants used for donor-recipient deconvolution with somatic nuclear mutations
variants.donor.recipient <- c(
  "NF1:chr17:29559932:C/A", "TP53:chr17:7579801:G/C", "SF3B1:chr2:198267770:G/GAA",
  "NPM1:chr5:170837457:A/G", "KDM6A:chrX:44938563:G/A", "DNMT3A:chr2:25463483:G/A",
  "BRAF:chr7:140449071:C/G", "FLT3:chr13:28592546:T/C", "NF1:chr17:29483195:G/C"
)

# variants recipient
variants.1 <- c("NF1:chr17:29559932:C/A", "TP53:chr17:7579801:G/C", "TP53:chr17:7578115:T/C")
# variants donor
variants.2 <- c(
  "SF3B1:chr2:198267770:G/GAA",
  "NPM1:chr5:170837457:A/G", "KDM6A:chrX:44938563:G/A", "DNMT3A:chr2:25463483:G/A",
  "BRAF:chr7:140449071:C/G", "FLT3:chr13:28592546:T/C", "NF1:chr17:29483195:G/C"
)

VARIANT_CUTOFF <- 10
variant.df <- data.frame(
  variants.1 = rowMeans(vafs[, variants.1]),
  variants.2 = rowMeans(vafs[, variants.2]),
  barcode = rownames(vafs)
)
variant.df$annotation <- "none"
variant.df$annotation[which(variant.df$variants.1 > VARIANT_CUTOFF & variant.df$variants.2 > VARIANT_CUTOFF)] <- "doublet"
variant.df$annotation[which(variant.df$variants.1 > VARIANT_CUTOFF & variant.df$variants.2 < VARIANT_CUTOFF)] <- "variant1"
variant.df$annotation[which(variant.df$variants.1 < VARIANT_CUTOFF & variant.df$variants.2 > VARIANT_CUTOFF)] <- "variant2"

# scatter plot illustrating donor-recipient deconvolution
p <- ggplot(variant.df, aes(x = variants.1, y = variants.2, color = annotation)) +
  geom_point(size = 0.5) +
  geom_vline(xintercept = 10) +
  geom_hline(yintercept = 10) +
  scale_x_continuous("%VAF variants 1 (n=3)", limits = c(0, 100)) +
  scale_y_continuous("%VAF variants 2 (n=7)", limits = c(0, 100)) +
  scale_color_manual(values = c("none" = "grey", "doublet" = "black", "variant1" = "orange", "variant2" = "purple")) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title = element_text("Arial", size = 10, color = "black"),
    axis.text = element_text("Arial", size = 10, color = "black")
  )
p <- ggMarginal(p, type = "histogram", xparams = list(bins = 100), yparams = list(bins = 100), groupFill = F, groupColour = F)
ggsave("./figure_coevolution/AML/figures/20220112_AML1026_variants_scatter_plot.svg", width = 2, height = 2, plot = p)

so.13 <- AddMetaData(so.13, metadata = variant.df)
saveRDS(object = so.13, file = "./data/objects/20220112_AML1026.rds")

# heatmap demonstrating donor-recipient deconvolution
ha <- HeatmapAnnotation(
  variant = c(
    variant.df$annotation[which(variant.df$annotation == "variant1")],
    variant.df$annotation[which(variant.df$annotation == "variant2")]
  ),
  col = list(variant = c("variant1" = "orange", "variant2" = "purple")),
  annotation_legend_param = list(variant = list(title = "Individual")), border = T
)
col_fun <- circlize::colorRamp2(breaks = seq(0, 1, 1 / 8), colors = BuenColors::jdb_palette(name = "solar_rojos", n = 9))
svglite::svglite("./figure_coevolution/AML/figures/20220112_AML1026_variants_heatmap.svg", width = 3, height = 2)
Heatmap(
  t(vafs[
    c(
      variant.df$barcode[which(variant.df$annotation == "variant1")],
      variant.df$barcode[which(variant.df$annotation == "variant2")]
    ),
    c(variants.1, variants.2)
  ]),
  top_annotation = ha, raster_quality = 10, use_raster = T, border = T,
  cluster_rows = T, cluster_columns = F, show_column_names = F, show_row_names = F, show_column_dend = F, show_row_dend = F
)
dev.off()

# track celltype-specific chimerism across both samples
variant.df$manual.cluster <- so.13$manual.cluster[rownames(variant.df)]
variant.df$sample <- stringr::str_split_fixed(variant.df$barcode, pattern = "#", n = 2)[, 1]

statistics.df <- data.frame()
for (celltype in unique(variant.df$manual.cluster)) {
  for (sample in unique(variant.df$sample)) {
    variant1 <- length(which(variant.df$manual.cluster == celltype & variant.df$sample == sample & variant.df$annotation == "variant1"))
    variant2 <- length(which(variant.df$manual.cluster == celltype & variant.df$sample == sample & variant.df$annotation == "variant2"))
    statistics.df <- rbind(statistics.df, data.frame(
      celltype = celltype, sample = sample, variant1 = variant1, variant2 = variant2,
      chimerism = variant2 / (variant1 + variant2)
    ))
  }
}
ggplot(data = reshape2::dcast(statistics.df, celltype ~ sample, value.var = "chimerism"), aes(x = `1026.1`, y = `1026.3`)) +
  geom_point()

cluster.colors <- c(
  "HSC" = RColorBrewer::brewer.pal(name = "YlOrRd", n = 5)[5],
  "Progenitor" = RColorBrewer::brewer.pal(name = "YlOrRd", n = 5)[4],
  "Mono" = RColorBrewer::brewer.pal(name = "YlOrRd", n = 5)[3],
  "Erythroid" = RColorBrewer::brewer.pal(name = "YlOrRd", n = 5)[2],
  "CD4" = "lightblue",
  "CD8" = "darkblue",
  "Plasma" = "darkgreen"
)

ggplot(statistics.df, aes(x = sample, y = 100 * chimerism, group = celltype, color = celltype)) +
  geom_point(size = 0.5) +
  geom_line() +
  scale_color_manual(values = cluster.colors) +
  scale_x_discrete(labels = c("Screening", "On-treatment")) +
  scale_y_continuous("% chimerism") +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    axis.title = element_text("Arial", size = 10, color = "black"),
    axis.text = element_text("Arial", size = 10, color = "black")
  )
ggsave("./figure_coevolution/AML/figures/20220112_AML1026_variants_kinetics.svg", width = 1, height = 2.5)

###
# track somatic mutations
###
RHP.mutations <- c("ASXL1" = "ASXL1:chr20:31022441:A/G", "NRAS" = "NRAS:chr1:115258745:C/G", "SF3B1" = "SF3B1:chr2:198266512:C/T")
vafs.RHP <- vafs[, RHP.mutations]
colnames(vafs.RHP) <- names(RHP.mutations)
depths.RHP <- depths[, RHP.mutations]

vafs.RHP <- vafs.RHP[complete.cases(vafs.RHP), ]
vafs.RHP$manual.cluster <- so.13$manual.cluster[rownames(vafs.RHP)]
vafs.RHP$annotation <- so.13$annotation[rownames(vafs.RHP)]

vafs.RHP$UMAP1 <- so.13@reductions[["umap"]]@cell.embeddings[rownames(vafs.RHP), "UMAP_1"]
vafs.RHP$UMAP2 <- so.13@reductions[["umap"]]@cell.embeddings[rownames(vafs.RHP), "UMAP_2"]

# visualize NRAS mutation on UMAP
p <- ggplot() +
  geom_point(data = vafs.RHP, aes(x = UMAP1, y = UMAP2), color = "grey90", size = 0.5) +
  geom_point(data = vafs.RHP[which(vafs.RHP$NRAS > 10), ], aes(x = UMAP1, y = UMAP2, color = NRAS), size = 0.5) +
  scale_color_gradientn(colours = BuenColors::jdb_palette(name = "solar_rojos")) +
  theme_classic() +
  NoAxes() +
  NoLegend()
ggsave("./figure_coevolution/AML/figures/20220113_AML1026_NRAS_UMAP.png", width = 2, height = 2, dpi = 600)

# visualize SF3B1 mutation on UMAP
p <- ggplot() +
  geom_point(data = vafs.RHP, aes(x = UMAP1, y = UMAP2), color = "grey90", size = 0.5) +
  geom_point(data = vafs.RHP[which(vafs.RHP$SF3B1 > 10), ], aes(x = UMAP1, y = UMAP2, color = SF3B1), size = 0.5) +
  scale_color_gradientn(colours = BuenColors::jdb_palette(name = "solar_rojos")) +
  theme_classic() +
  NoAxes() +
  NoLegend()
ggsave("./figure_coevolution/AML/figures/20220113_AML1026_SF3B1_UMAP.png", width = 2, height = 2, dpi = 600)

# visualize NRAS and SF3B1 mutation across donor- and recipient-derived cell populations
vafs.RHP$NRAS[which(vafs.RHP$annotation == "variant2")] <- -vafs.RHP$NRAS[which(vafs.RHP$annotation == "variant2")]
vafs.RHP$SF3B1[which(vafs.RHP$annotation == "variant2")] <- -vafs.RHP$SF3B1[which(vafs.RHP$annotation == "variant2")]

labels <- vafs.RHP %>%
  filter(annotation %in% c("variant1", "variant2")) %>%
  group_by(manual.cluster, annotation) %>%
  filter(NRAS != 0) %>%
  count()
labels$y <- ifelse(labels$annotation == "variant1", 120, -100)

p <- ggplot(data = vafs.RHP[which(vafs.RHP$annotation %in% c("variant1", "variant2")), ], aes(x = manual.cluster, y = NRAS)) +
  geom_hline(yintercept = 0) +
  geom_violin(aes(fill = manual.cluster)) +
  geom_jitter(size = 0.5, aes(color = manual.cluster)) +
  scale_y_continuous(limits = c(-100, 120), breaks = c(-100, -50, 0, 50, 100), labels = c(100, 50, 0, 50, 100)) +
  scale_fill_manual(values = cluster.colors) +
  scale_color_manual(values = cluster.colors) +
  geom_label(data = labels, aes(x = manual.cluster, y = y, label = paste0("n=", n)), label.size = 0, size = 2) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    axis.title = element_text("Arial", size = 10, color = "black"),
    axis.text = element_text("Arial", size = 10, color = "black")
  )
ggsave("./figure_coevolution/AML/figures/20220113_AML1026_NRAS.svg", width = 2.6, height = 2.1, plot = p)

labels <- vafs.RHP %>%
  filter(annotation %in% c("variant1", "variant2")) %>%
  group_by(manual.cluster, annotation) %>%
  filter(SF3B1 != 0) %>%
  count()
labels$y <- ifelse(labels$annotation == "variant1", 120, -100)

p <- ggplot(data = vafs.RHP[which(vafs.RHP$annotation %in% c("variant1", "variant2")), ], aes(x = manual.cluster, y = SF3B1)) +
  geom_hline(yintercept = 0) +
  geom_violin(aes(fill = manual.cluster)) +
  geom_jitter(size = 0.5, aes(color = manual.cluster)) +
  scale_y_continuous(limits = c(-100, 120), breaks = c(-100, -50, 0, 50, 100), labels = c(100, 50, 0, 50, 100)) +
  scale_fill_manual(values = cluster.colors) +
  scale_color_manual(values = cluster.colors) +
  geom_label(data = labels, aes(x = manual.cluster, y = y, label = paste0("n=", n)), label.size = 0, size = 2) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    axis.title = element_text("Arial", size = 10, color = "black"),
    axis.text = element_text("Arial", size = 10, color = "black")
  )
ggsave("./figure_coevolution/AML/figures/20220113_AML1026_SF3B1.svg", width = 2.6, height = 2.1, plot = p)

# kinetics of NRAS and SF3B1 mutation across both samples
vafs.RHP$sample <- stringr::str_split_fixed(rownames(vafs.RHP), pattern = "#", n = 2)[, 1]

statistics.df <- data.frame()
for (celltype in unique(vafs.RHP$manual.cluster)) {
  NRAS.1 <- length(which(vafs.RHP$sample == "1026.1" & vafs.RHP$manual.cluster == celltype & vafs.RHP$annotation == "variant1" & vafs.RHP$NRAS != 0)) /
    length(which(vafs.RHP$sample == "1026.1" & vafs.RHP$manual.cluster == celltype & vafs.RHP$annotation == "variant1"))
  NRAS.3 <- length(which(vafs.RHP$sample == "1026.3" & vafs.RHP$manual.cluster == celltype & vafs.RHP$annotation == "variant1" & vafs.RHP$NRAS != 0)) /
    length(which(vafs.RHP$sample == "1026.3" & vafs.RHP$manual.cluster == celltype & vafs.RHP$annotation == "variant1"))
  SF3B1.1 <- length(which(vafs.RHP$sample == "1026.1" & vafs.RHP$manual.cluster == celltype & vafs.RHP$annotation == "variant1" & vafs.RHP$SF3B1 != 0)) /
    length(which(vafs.RHP$sample == "1026.1" & vafs.RHP$manual.cluster == celltype & vafs.RHP$annotation == "variant1"))
  SF3B1.3 <- length(which(vafs.RHP$sample == "1026.3" & vafs.RHP$manual.cluster == celltype & vafs.RHP$annotation == "variant1" & vafs.RHP$SF3B1 != 0)) /
    length(which(vafs.RHP$sample == "1026.3" & vafs.RHP$manual.cluster == celltype & vafs.RHP$annotation == "variant1"))
  statistics.df <- rbind(statistics.df, data.frame(celltype = celltype, NRAS.1 = NRAS.1, NRAS.3 = NRAS.3, SF3B1.1 = SF3B1.1, SF3B1.3 = SF3B1.3))
}

boo <- reshape2::melt(statistics.df, id.vars = "celltype")
boo$gene <- stringr::str_split_fixed(boo$variable, pattern = "\\.", n = 2)[, 1]
boo$gene.manual.cluster <- paste0(boo$gene, ".", boo$celltype)
ggplot(
  boo[which(boo$celltype %in% c("HSC", "Progenitor", "Mono", "Erythroid", "Plasma")), ],
  aes(x = variable, y = 100 * value)
) +
  geom_point(aes(color = celltype)) +
  geom_line(aes(group = gene.manual.cluster, color = celltype)) +
  scale_x_discrete(labels = c("Screening", "On treatment", "Screening", "On treatment")) +
  scale_y_continuous("%recipient positive") +
  scale_color_manual(values = cluster.colors) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    axis.title = element_text("Arial", size = 10, color = "black"),
    axis.text = element_text("Arial", size = 10, color = "black")
  )
ggsave("./figure_coevolution/AML/figures/20220113_AML1026_NRAS_SF3B1_kinetics.svg", width = 1.3, height = 2)

###
# T cell chimerism
###
so.13.Tcell <- readRDS("./data/coevolution/objects/20220112_AML1026_Tcell.rds")

cluster.colors.T <- c(
  "CD4" = RColorBrewer::brewer.pal(name = "Blues", n = 6)[2],
  "CD4 naive" = RColorBrewer::brewer.pal(name = "Blues", n = 6)[3],
  "CD4 CM" = RColorBrewer::brewer.pal(name = "Blues", n = 6)[4],
  "CD4 CM activated" = RColorBrewer::brewer.pal(name = "Blues", n = 6)[5],
  "CD4 EM" = RColorBrewer::brewer.pal(name = "Blues", n = 6)[6],
  "CD8" = RColorBrewer::brewer.pal(name = "Reds", n = 3)[2],
  "CD8 memory" = RColorBrewer::brewer.pal(name = "Reds", n = 3)[3],
  "NK" = "black"
)

T.cell.annotation <- data.frame(
  sample = so.13.Tcell$orig.ident,
  manual.cluster = so.13.Tcell$manual.cluster,
  annotation = so.13.Tcell$annotation
)
T.cell.annotation <- T.cell.annotation %>% filter(annotation %in% c("variant1", "variant2"))

T.cell.statistics <- as.data.frame(T.cell.annotation %>% group_by(manual.cluster) %>% count(annotation))
T.cell.statistics$freq <- 0
for (i in seq(1, nrow(T.cell.statistics))) {
  T.cell.statistics$freq[i] <- T.cell.statistics$n[i] / sum(T.cell.statistics$n[which(T.cell.statistics$manual.cluster == T.cell.statistics$manual.cluster[i])])
}
T.cell.statistics <- T.cell.statistics %>% filter(annotation == "variant2")

ggplot(T.cell.statistics, aes(x = manual.cluster, y = 100 * freq)) +
  geom_col(aes(fill = manual.cluster)) +
  scale_y_continuous("%donor chimerism", limits = c(0, 100)) +
  scale_fill_manual(values = cluster.colors.T) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    axis.title = element_text("Arial", size = 10, color = "black"),
    axis.text = element_text("Arial", size = 10, color = "black")
  )
ggsave("./figure_coevolution/AML/figures/20220113_AML1026_Tcell_chimerism.svg", width = 1.7, height = 2.5)
