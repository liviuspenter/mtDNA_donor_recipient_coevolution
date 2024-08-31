# create heatmap of CLL6 (called CLL2)

library(ComplexHeatmap)
library(dplyr)
library(ggplot2)
library(Seurat)
library(Signac)

`%ni%` <- Negate(`%in%`)

col_fun <- circlize::colorRamp2(breaks = c(-1, seq(0, 100, 100 / 8)), colors = c("grey", BuenColors::jdb_palette(name = "solar_rojos")))

# call mtDNA mutations using mgatk
CLL1.mito <- ReadMGATK("./data/coevolution/CLL1.hg38.mgatk/")
colnames(CLL1.mito$counts) <- gsub(paste0("CLL1#", colnames(CLL1.mito$counts)), pattern = "-1", replacement = "")
rownames(CLL1.mito$depth) <- gsub(paste0("CLL1#", rownames(CLL1.mito$depth)), pattern = "-1", replacement = "")
CLL4.mito <- ReadMGATK("./data/coevolution/CLL4.mgatk/")
colnames(CLL4.mito$counts) <- gsub(paste0("CLL4#", colnames(CLL4.mito$counts)), pattern = "-1", replacement = "")
rownames(CLL4.mito$depth) <- gsub(paste0("CLL4#", rownames(CLL4.mito$depth)), pattern = "-1", replacement = "")
CLL.mito <- CLL1.mito
CLL.mito$depth <- rbind(CLL.mito$depth, CLL4.mito$depth)
CLL.mito$counts <- cbind(CLL.mito$counts, CLL4.mito$counts)

CLL.so <- CreateAssayObject(counts = CLL.mito$counts)
CLL.mgatk <- CreateSeuratObject(counts = CLL.so, project = "CLL.mgatk", assay = "mito")
rm(CLL.so)
CLL.mgatk <- AddMetaData(CLL.mgatk, metadata = CLL.mito$depth, col.name = "mtDNA_depth")

CLL.mgatk <- AlleleFreq(
  object = CLL.mgatk, variants = c(
    "3526G>A", "3532A>C", "3830T>C", "930G>A", "7205C>T", "2623A>G",
    "9554G>A", "13030T>C", "13471G>A", "9540C>T"
  ),
  assay = "mito"
)
CLL.mito.vaf <- as.data.frame(t(as.data.frame(GetAssayData(CLL.mgatk[["alleles"]]))))

so.12 <- readRDS("./data/coevolution/objects/20220504_CLL34_filtered.rds")

vafs.1 <- as.data.frame(data.table::fread("./data/coevolution/CLL1.whitelist/AF.csv"))
vafs.1$Barcode <- gsub(paste0("CLL1#", vafs.1$Barcode), pattern = "-1", replacement = "")
rownames(vafs.1) <- vafs.1$Barcode
vafs.1 <- vafs.1[, -c(1, 2)]
vafs.2 <- as.data.frame(data.table::fread("./data/coevolution/CLL4.whitelist/AF.csv"))
vafs.2$Barcode <- gsub(paste0("CLL4#", vafs.2$Barcode), pattern = "-1", replacement = "")
rownames(vafs.2) <- vafs.2$Barcode
vafs.2 <- vafs.2[, -c(1, 2)]

depth.1 <- as.data.frame(data.table::fread("./data/coevolution/CLL1.whitelist/DP.csv"))
depth.1$Barcode <- gsub(paste0("CLL1#", depth.1$Barcode), pattern = "-1", replacement = "")
rownames(depth.1) <- depth.1$Barcode
depth.1 <- depth.1[, -c(1, 2)]
depth.2 <- as.data.frame(data.table::fread("./data/coevolution/CLL4.whitelist/DP.csv"))
depth.2$Barcode <- gsub(paste0("CLL4#", depth.2$Barcode), pattern = "-1", replacement = "")
rownames(depth.2) <- depth.2$Barcode
depth.2 <- depth.2[, -c(1, 2)]

vafs <- dplyr::bind_rows(vafs.1, vafs.2)
depth <- dplyr::bind_rows(depth.1, depth.2)

annotation1 <- read.csv2("./data/coevolution/CLL1/20220506_deconvolution.csv", row.names = 1)
annotation4 <- read.csv2("./data/coevolution/CLL4/20220506_deconvolution.csv", row.names = 1)
annotation1$bc <- gsub(annotation1$bc, pattern = "-1", replacement = "")
annotation4$bc <- gsub(annotation4$bc, pattern = "-1", replacement = "")

cells.donor.2 <- annotation1$bc[which(annotation1$annotation == "CLL2.donor" & annotation1$bc %in% rownames(vafs))]
cells.recipient.2 <- annotation1$bc[which(annotation1$annotation == "CLL2.recipient" & annotation1$bc %in% rownames(vafs))]

cells.1.CLL <- colnames(so.12)[which(colnames(so.12) %in% annotation4$bc[which(annotation4$annotation == "CLL2.recipient")] &
  so.12$manual.cluster %in% c("CLL", "Cluster 6"))]
cells.1.immune <- colnames(so.12)[which(colnames(so.12) %in% annotation4$bc[which(annotation4$annotation == "CLL2.recipient")] &
  so.12$manual.cluster %ni% c("CLL", "Cluster 6"))]
cells.1 <- c(cells.1.CLL, cells.1.immune)

depth <- depth[c(cells.1, cells.donor.2, cells.recipient.2), ]
vafs <- vafs[c(cells.1, cells.donor.2, cells.recipient.2), ]

vafs.auto <- vafs[, which(!grepl("chrM", colnames(vafs)))]
vafs.auto[vafs.auto < 10] <- 0
depth <- depth[rownames(vafs.auto), ]
depth <- depth[, colnames(vafs.auto)]
vafs.auto[is.na(depth) | depth == 0] <- -1



variant.df <- data.frame(
  vmr = -log10(apply(vafs.auto, 2, FUN = function(x) {
    mean(x) / var(x)
  })),
  vaf = apply(vafs.auto, 2, FUN = function(x) {
    median(x[which(x != 0)])
  }),
  freq = apply(vafs.auto, 2, FUN = function(x) {
    length(which(x != 0))
  })
)

vafs.auto <- vafs.auto[intersect(rownames(vafs.auto), rownames(CLL.mito.vaf)), ]

vafs.mito <- 100 * CLL.mito.vaf[rownames(vafs.auto), ]

vafs <- cbind(vafs.auto, vafs.mito)

relevant.mutations <- unique(c(
  rownames(variant.df)[which(variant.df$vmr > 0 & variant.df$vaf > 25)],
  colnames(vafs.mito)
))

ha <- columnAnnotation(
  fraction = c(
    rep("cells.1.CLL", length(cells.1.CLL)),
    rep("cells.2.CLL", length(cells.recipient.2)),
    rep("cells.1.immune", length(cells.1.immune)),
    rep("cells.2.immune", length(cells.donor.2))
  ),
  cluster = so.12$manual.cluster[c(cells.1, cells.recipient.2, cells.donor.2)],
  col = list(
    "fraction" = c(
      "cells.1.CLL" = "orange",
      "cells.2.CLL" = "orange",
      "cells.1.immune" = "black",
      "cells.2.immune" = "black"
    ),
    "cluster" = c(
      "CLL" = "orange", "CD4 T cell" = "lightblue",
      "CD8 T cell" = "darkblue", "NK" = "purple",
      "Myeloid" = "darkgreen", "Cluster 6" = "orange"
    )
  ),
  simple_anno_size = unit(5, "pt"), border = T
)

svglite::svglite("./figure_coevolution/CLL/figures/CLL6/heatmaps/20230427_CLL6_raw.svg", width = 15, height = 15)
Heatmap(t(vafs[c(cells.1.CLL, cells.recipient.2, cells.1.immune, cells.donor.2), relevant.mutations]),
  show_row_names = T, show_column_names = F,
  show_row_dend = F, show_column_dend = F, row_names_side = "left", cluster_columns = F, cluster_rows = T,
  row_names_gp = gpar(fontsize = 8),
  use_raster = T, raster_quality = 10, column_split = c(
    rep("cells.1.CLL", length(cells.1.CLL)),
    rep("cells.2.CLL", length(cells.recipient.2)),
    rep("cells.1.immune", length(cells.1.immune)),
    rep("cells.2.immune", length(cells.donor.2))
  ),
  col = col_fun, border = T, row_split = factor(ifelse(grepl("chrM", relevant.mutations), "mito", "auto"), levels = c("auto", "mito")),
  top_annotation = ha
)
dev.off()


### make curated heatmap with custom ordering
source("./R/20200328_mtscatac_seq.R")
donor.auto <- c()
donor.mito <- c()
recipient.auto <- c("ESRRG:chr1:216880662:G/A")
CLL.common.auto <- c(
  "ERICH3:chr1:75102115:G/C", "PLEC:chr8:144994678:CCT/C", "WASF3:chr13:27255288:G/T", "DNAJB7:chr22:41257549:A/C",
  "BPIFA2:chr20:31761967:G/C", "ZNF292:chr6:87967797:CA/C", "DCBLD2:chr3:98600507:G/A", "CHD7:chr8:61765460:A/AGCT"
)
CLL.common.mito <- c()
CLL.auto <- c("PKDREJ:chr22:46658486:A/C", "TP53:chr17:7577124:C/T")
CLL.mito <- c("3830T>C", "3526G>A", "930G>A", "2623A>G", "9554G>A")

vafs <- dplyr::bind_rows(vafs.1, vafs.2)
depth <- dplyr::bind_rows(depth.1, depth.2)

depth <- depth[c(cells.1, cells.donor.2, cells.recipient.2), ]
vafs <- vafs[c(cells.1, cells.donor.2, cells.recipient.2), ]

vafs.auto <- vafs[, c(CLL.auto, CLL.common.auto)]
vafs.auto[vafs.auto < 10] <- 0
depth <- depth[rownames(vafs.auto), ]
depth <- depth[, colnames(vafs.auto)]

vafs.auto[is.na(depth) | depth == 0] <- NA
vafs.auto <- vafs.auto[complete.cases(vafs.auto), ]

vafs.auto <- vafs.auto[intersect(rownames(vafs.auto), rownames(CLL.mito.vaf)), ]

# plot from 0 to 10%
vafs.mito <- 1000 * CLL.mito.vaf[rownames(vafs.auto), ]

vafs <- cbind(vafs.auto, vafs.mito)


cells.2.donor <- annotation1$bc[which(annotation1$annotation == "CLL2.donor" & annotation1$bc %in% rownames(vafs))]
cells.2.recipient <- annotation1$bc[which(annotation1$annotation == "CLL2.recipient" & annotation1$bc %in% rownames(vafs))]
cells.1.CLL <- annotation4$bc[which(annotation4$annotation == "CLL2.recipient" & annotation4$bc %in% rownames(vafs) &
  annotation4$bc %in% colnames(so.12)[which(so.12$manual.cluster == "CLL")])]
cells.1.immune <- annotation4$bc[which(annotation4$annotation == "CLL2.recipient" & annotation4$bc %in% rownames(vafs) &
  annotation4$bc %in% colnames(so.12)[which(so.12$manual.cluster %ni% c("CLL", "Cluster 6"))])]

cluster.information.1 <- cluster_relevant_mutations(t(vafs[cells.1.CLL, c(CLL.auto, CLL.mito)]), dims = 2, k_param = 5)
cells.1.CLL <- colnames(cluster.information.1[[1]])
cells.1.CLL <- cells.1.CLL[order(vafs[cells.1.CLL, "TP53:chr17:7577124:C/T"], decreasing = T)]
cluster.information.2 <- cluster_relevant_mutations(t(vafs[cells.2.recipient, c(CLL.auto, CLL.mito)]), dims = 2, k_param = 5)
cells.2.recipient <- colnames(cluster.information.2[[1]])
cells.2.recipient <- cells.2.recipient[order(vafs[cells.2.recipient, "TP53:chr17:7577124:C/T"], decreasing = T)]
cells.1.immune <- names(sort(so.12$manual.cluster[cells.1.immune]))

ha <- columnAnnotation(
  cluster = so.12$manual.cluster[c(cells.1.CLL, cells.1.immune, cells.2.recipient, cells.2.donor)],
  col = list("cluster" = c(
    "CLL" = "orange", "CD4 T cell" = "lightblue", "CD8 T cell" = "darkblue", "NK" = "purple", "Myeloid" = "darkgreen",
    "Cluster 6" = "black"
  )),
  simple_anno_size = unit(5, "pt"), border = T
)

svglite::svglite("./figure_coevolution/CLL/figures/CLL6/heatmaps/20230427_CLL6_curated.svg", width = 7, height = 2)
Heatmap(
  t(vafs[
    c(cells.1.CLL, cells.1.immune, cells.2.recipient, cells.2.donor),
    c(CLL.common.auto, CLL.common.mito, CLL.auto, CLL.mito, donor.mito)
  ]),
  show_row_names = T, show_column_names = F,
  show_row_dend = F, show_column_dend = F, row_names_side = "left", cluster_columns = F,
  row_names_gp = gpar(fontsize = 8), cluster_rows = F,
  use_raster = T, raster_quality = 10,
  row_split = factor(
    c(
      rep("CLL common auto", length(CLL.common.auto)),
      rep("CLL common mito", length(CLL.common.mito)),
      rep("CLL auto", length(CLL.auto)),
      rep("CLL mito", length(CLL.mito)),
      rep("Donor mito", length(donor.mito))
    ),
    levels = c("CLL common auto", "CLL common mito", "CLL auto", "CLL mito", "Donor mito")
  ),
  column_split = factor(
    c(
      rep("CLL pre-FCR", length(cells.1.CLL)),
      rep("Immune cells pre-FCR", length(cells.1.immune)),
      rep("CLL post-RIC", length(cells.2.recipient)),
      rep("Donor", length(cells.2.donor))
    ),
    levels = c("CLL pre-FCR", "CLL post-RIC", "Immune cells pre-FCR", "Donor")
  ),
  column_title_gp = gpar(fontsize = 8),
  row_title_gp = gpar(fontsize = 8),
  col = col_fun, border = T,
  top_annotation = ha
)
dev.off()

### plot individual mutations

ggplot() +
  geom_point(data = vafs[cells.1.CLL, ], aes(x = `TP53:chr17:7577124:C/T`, y = `2623A>G` / 10), color = "blue", size = 0.5) +
  scale_x_continuous("% TP53V272L", limits = c(0, 100)) +
  scale_y_continuous("% 2623A>G", limits = c(0, 100)) +
  theme_classic() +
  theme(
    axis.text = element_text("Arial", size = 8, color = "black"),
    axis.title = element_text("Arial", size = 8, color = "black")
  )
ggsave("./figure_coevolution/CLL/figures/CLL6/plots/20230427_CLL6_pre_FCR_TP53_2623A>G.svg", width = 1.1, height = 1.1)

ggplot() +
  geom_point(data = vafs[cells.2.recipient, ], aes(x = `TP53:chr17:7577124:C/T`, y = `2623A>G` / 10), color = "firebrick", size = 0.5) +
  scale_x_continuous("% TP53V272L", limits = c(0, 100)) +
  scale_y_continuous("% 2623A>G", limits = c(0, 100)) +
  theme_classic() +
  theme(
    axis.text = element_text("Arial", size = 8, color = "black"),
    axis.title = element_text("Arial", size = 8, color = "black")
  )
ggsave("./figure_coevolution/CLL/figures/CLL6/plots/20230427_CLL6_post_HSCT_TP53_2623A>G.svg", width = 1.1, height = 1.1)

ggplot() +
  geom_point(data = vafs[cells.1.CLL, ], aes(x = `PKDREJ:chr22:46658486:A/C`, y = `2623A>G` / 10), color = "blue", size = 0.5) +
  scale_x_continuous("% PKDREJV245P", limits = c(0, 100)) +
  scale_y_continuous("% 2623A>G", limits = c(0, 100)) +
  theme_classic() +
  theme(
    axis.text = element_text("Arial", size = 8, color = "black"),
    axis.title = element_text("Arial", size = 8, color = "black")
  )
ggsave("./figure_coevolution/CLL/figures/CLL6/plots/20230427_CLL6_pre_FCR_PKDREJ_2623A>G.svg", width = 1.1, height = 1.1)

ggplot() +
  geom_point(data = vafs[cells.2.recipient, ], aes(x = `PKDREJ:chr22:46658486:A/C`, y = `2623A>G` / 10), color = "firebrick", size = 0.5) +
  scale_x_continuous("% PKDREJV245P", limits = c(0, 100)) +
  scale_y_continuous("% 2623A>G", limits = c(0, 100)) +
  theme_classic() +
  theme(
    axis.text = element_text("Arial", size = 8, color = "black"),
    axis.title = element_text("Arial", size = 8, color = "black")
  )
ggsave("./figure_coevolution/CLL/figures/CLL6/plots/20230427_CLL6_post_HSCT_PKDREJ_2623A>G.svg", width = 1.1, height = 1.1)




ggplot() +
  geom_point(data = vafs[cells.1.CLL, ], aes(x = `TP53:chr17:7577124:C/T`, y = `3526G>A` / 10), color = "blue", size = 0.5) +
  scale_x_continuous("% TP53V272L", limits = c(0, 100)) +
  scale_y_continuous("% 3526G>A", limits = c(0, 100)) +
  theme_classic() +
  theme(
    axis.text = element_text("Arial", size = 8, color = "black"),
    axis.title = element_text("Arial", size = 8, color = "black")
  )
ggsave("./figure_coevolution/CLL/figures/CLL6/plots/20230427_CLL6_pre_FCR_TP53_3526G>A.svg", width = 1.1, height = 1.1)

ggplot() +
  geom_point(data = vafs[cells.2.recipient, ], aes(x = `TP53:chr17:7577124:C/T`, y = `3526G>A` / 10), color = "firebrick", size = 0.5) +
  scale_x_continuous("% TP53V272L", limits = c(0, 100)) +
  scale_y_continuous("% 3526G>A", limits = c(0, 100)) +
  theme_classic() +
  theme(
    axis.text = element_text("Arial", size = 8, color = "black"),
    axis.title = element_text("Arial", size = 8, color = "black")
  )
ggsave("./figure_coevolution/CLL/figures/CLL6/plots/20230427_CLL6_post_HSCT_TP53_3526G>A.svg", width = 1.1, height = 1.1)

ggplot() +
  geom_point(data = vafs[cells.1.CLL, ], aes(x = `PKDREJ:chr22:46658486:A/C`, y = `3526G>A` / 10), color = "blue", size = 0.5) +
  scale_x_continuous("% PKDREJV245P", limits = c(0, 100)) +
  scale_y_continuous("% 3526G>A", limits = c(0, 100)) +
  theme_classic() +
  theme(
    axis.text = element_text("Arial", size = 8, color = "black"),
    axis.title = element_text("Arial", size = 8, color = "black")
  )
ggsave("./figure_coevolution/CLL/figures/CLL6/plots/20230427_CLL6_pre_FCR_PKDREJ_3526G>A.svg", width = 1.1, height = 1.1)

ggplot() +
  geom_point(data = vafs[cells.2.recipient, ], aes(x = `PKDREJ:chr22:46658486:A/C`, y = `3526G>A` / 10), color = "firebrick", size = 0.5) +
  scale_x_continuous("% PKDREJV245P", limits = c(0, 100)) +
  scale_y_continuous("% 3526G>A", limits = c(0, 100)) +
  theme_classic() +
  theme(
    axis.text = element_text("Arial", size = 8, color = "black"),
    axis.title = element_text("Arial", size = 8, color = "black")
  )
ggsave("./figure_coevolution/CLL/figures/CLL6/plots/20230427_CLL6_post_HSCT_PKDREJ_3526G>A.svg", width = 1.1, height = 1.1)



ggplot() +
  geom_point(data = vafs[cells.1.CLL, ], aes(x = `TP53:chr17:7577124:C/T`, y = `3830T>C` / 10), color = "blue", size = 0.5) +
  scale_x_continuous("% TP53V272L", limits = c(0, 100)) +
  scale_y_continuous("% 3830T>C", limits = c(0, 100)) +
  theme_classic() +
  theme(
    axis.text = element_text("Arial", size = 8, color = "black"),
    axis.title = element_text("Arial", size = 8, color = "black")
  )
ggsave("./figure_coevolution/CLL/figures/CLL6/plots/20230427_CLL6_pre_FCR_TP53_3830T>C.svg", width = 1.1, height = 1.1)

ggplot() +
  geom_point(data = vafs[cells.2.recipient, ], aes(x = `TP53:chr17:7577124:C/T`, y = `3830T>C` / 10), color = "firebrick", size = 0.5) +
  scale_x_continuous("% TP53V272L", limits = c(0, 100)) +
  scale_y_continuous("% 3830T>C", limits = c(0, 100)) +
  theme_classic() +
  theme(
    axis.text = element_text("Arial", size = 8, color = "black"),
    axis.title = element_text("Arial", size = 8, color = "black")
  )
ggsave("./figure_coevolution/CLL/figures/CLL6/plots/20230427_CLL6_post_HSCT_TP53_3830T>C.svg", width = 1.1, height = 1.1)

ggplot() +
  geom_point(data = vafs[cells.1.CLL, ], aes(x = `PKDREJ:chr22:46658486:A/C`, y = `3830T>C` / 10), color = "blue", size = 0.5) +
  scale_x_continuous("% PKDREJV245P", limits = c(0, 100)) +
  scale_y_continuous("% 3830T>C", limits = c(0, 100)) +
  theme_classic() +
  theme(
    axis.text = element_text("Arial", size = 8, color = "black"),
    axis.title = element_text("Arial", size = 8, color = "black")
  )
ggsave("./figure_coevolution/CLL/figures/CLL6/plots/20230427_CLL6_pre_FCR_PKDREJ_3830T>C.svg", width = 1.1, height = 1.1)

ggplot() +
  geom_point(data = vafs[cells.2.recipient, ], aes(x = `PKDREJ:chr22:46658486:A/C`, y = `3830T>C` / 10), color = "firebrick", size = 0.5) +
  scale_x_continuous("% PKDREJV245P", limits = c(0, 100)) +
  scale_y_continuous("% 3830T>C", limits = c(0, 100)) +
  theme_classic() +
  theme(
    axis.text = element_text("Arial", size = 8, color = "black"),
    axis.title = element_text("Arial", size = 8, color = "black")
  )
ggsave("./figure_coevolution/CLL/figures/CLL6/plots/20230427_CLL6_post_HSCT_PKDREJ_3830T>C.svg", width = 1.1, height = 1.1)



ggplot(vafs[cells.1.CLL, ], aes(x = `TP53:chr17:7577124:C/T`, y = `PKDREJ:chr22:46658486:A/C`)) +
  geom_point(size = 0.5, color = "blue") +
  scale_x_continuous("% TP53V272L", limits = c(0, 100)) +
  scale_y_continuous("% PKDREJV245P", limits = c(0, 100)) +
  theme_classic() +
  theme(
    axis.text = element_text("Arial", size = 8, color = "black"),
    axis.title = element_text("Arial", size = 8, color = "black")
  )
ggsave("./figure_coevolution/CLL/figures/CLL6/plots/20230427_CLL6_pre_FCR_TP53_PKDREJ.svg", width = 1.1, height = 1.1)

ggplot(vafs[cells.2.recipient, ], aes(x = `TP53:chr17:7577124:C/T`, y = `PKDREJ:chr22:46658486:A/C`)) +
  geom_point(size = 0.5, color = "firebrick") +
  scale_x_continuous("% TP53V272L", limits = c(0, 100)) +
  scale_y_continuous("% PKDREJV245P", limits = c(0, 100)) +
  theme_classic() +
  theme(
    axis.text = element_text("Arial", size = 8, color = "black"),
    axis.title = element_text("Arial", size = 8, color = "black")
  )
ggsave("./figure_coevolution/CLL/figures/CLL6/plots/20230427_CLL6_post_HSCT_TP53_PKDREJ.svg", width = 1.1, height = 1.1)


ggplot(vafs[cells.1.CLL, ], aes(x = `3830T>C` / 10, y = `3526G>A` / 10)) +
  geom_point(size = 0.5, color = "blue") +
  scale_x_continuous("% 3830T>C", limits = c(0, 100)) +
  scale_y_continuous("% 3526G>A", limits = c(0, 100)) +
  theme_classic() +
  theme(
    axis.text = element_text("Arial", size = 8, color = "black"),
    axis.title = element_text("Arial", size = 8, color = "black")
  )
ggsave("./figure_coevolution/CLL/figures/CLL6/plots/20230427_CLL6_pre_FCR_3830T>C_3526G>A.svg", width = 1.1, height = 1.1)

ggplot(vafs[cells.2.recipient, ], aes(x = `3830T>C` / 10, y = `3526G>A` / 10)) +
  geom_point(size = 0.5, color = "firebrick") +
  scale_x_continuous("% 3830T>C", limits = c(0, 100)) +
  scale_y_continuous("% 3526G>A", limits = c(0, 100)) +
  theme_classic() +
  theme(
    axis.text = element_text("Arial", size = 8, color = "black"),
    axis.title = element_text("Arial", size = 8, color = "black")
  )
ggsave("./figure_coevolution/CLL/figures/CLL6/plots/20230427_CLL6_post_HSCT_3830T>C_3526G>A.svg", width = 1.1, height = 1.1)



TP53.pos <- which(vafs[, "TP53:chr17:7577124:C/T"] != 0)
PKDREJ.pos <- which(vafs[, "PKDREJ:chr22:46658486:A/C"] != 0)

ggplot() +
  geom_histogram(data = vafs[TP53.pos, ], aes(x = `3526G>A` / 10), fill = "firebrick") +
  geom_histogram(data = vafs[PKDREJ.pos, ], aes(x = `3526G>A` / 10), fill = "grey", binwidth = 5) +
  geom_vline(xintercept = 5) +
  scale_x_continuous("% 3526G>A") +
  scale_y_sqrt(breaks = c(1, 20, 100, 200, 500)) +
  theme_classic() +
  theme_classic() +
  theme(
    axis.text = element_text("Arial", size = 8, color = "black"),
    axis.title = element_text("Arial", size = 8, color = "black")
  )
ggsave("./figure_coevolution/CLL/figures/CLL6/plots/20240530_distribution_3526G>A.svg", width = 2, height = 1.1)
length(which(vafs[PKDREJ.pos, "3526G>A"] > 50)) / length(PKDREJ.pos)
length(which(vafs[TP53.pos, "3526G>A"] > 50)) / length(TP53.pos)

fisher.test(matrix(c(
  length(which(vafs[PKDREJ.pos, "3526G>A"] > 50)), length(PKDREJ.pos),
  length(which(vafs[TP53.pos, "3526G>A"] > 50)), length(TP53.pos)
), nrow = 2))

ggplot() +
  geom_histogram(data = vafs[TP53.pos, ], aes(x = `3830T>C` / 10), fill = "firebrick") +
  geom_histogram(data = vafs[PKDREJ.pos, ], aes(x = `3830T>C` / 10), fill = "grey", binwidth = 5) +
  geom_vline(xintercept = 5) +
  scale_x_continuous("% 3830T>C") +
  scale_y_sqrt(breaks = c(1, 20, 100, 200, 500)) +
  theme_classic() +
  theme_classic() +
  theme(
    axis.text = element_text("Arial", size = 8, color = "black"),
    axis.title = element_text("Arial", size = 8, color = "black")
  )
ggsave("./figure_coevolution/CLL/figures/CLL6/plots/20240530_distribution_3830T>C.svg", width = 2, height = 1.1)
length(which(vafs[PKDREJ.pos, "3830T>C"] > 50)) / length(PKDREJ.pos)
length(which(vafs[TP53.pos, "3830T>C"] > 50)) / length(TP53.pos)

fisher.test(matrix(c(
  length(which(vafs[PKDREJ.pos, "3830T>C"] > 50)), length(PKDREJ.pos),
  length(which(vafs[TP53.pos, "3830T>C"] > 50)), length(TP53.pos)
), nrow = 2))
