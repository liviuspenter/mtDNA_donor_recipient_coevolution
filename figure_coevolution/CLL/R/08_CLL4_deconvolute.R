# sample deconvolution of sample CLL3
# CLL3 contains CLL4 pre-FCR (called CLL1) and CLL6 pre-FCR (called CLL2)

library(ComplexHeatmap)
library(dplyr)
library(ggplot2)
library(Seurat)

col_fun <- circlize::colorRamp2(breaks = seq(0, 100, 100 / 8), colors = BuenColors::jdb_palette(name = "solar_rojos"))

so.12 <- readRDS("./data/coevolution/objects/20220504_CLL34.rds")

vafs.1 <- as.data.frame(data.table::fread("./data/coevolution/CLL4/AF.csv"))
vafs.1$Barcode <- paste0("CLL4#", vafs.1$Barcode)
rownames(vafs.1) <- vafs.1$Barcode
vafs.1 <- vafs.1[, -c(1, 2)]

boo <- vafs.1[sample(nrow(vafs.1), size = 1000), ]

ha <- columnAnnotation(
  cluster = so.12$manual.cluster[rownames(boo)],
  col = list("cluster" = c(
    "CLL" = "orange", "CD4 T cell" = "lightblue", "CD8 T cell" = "darkblue", "Myeloid" = "darkgreen",
    "Cluster 6" = "black"
  ))
)

Heatmap(t(boo),
  show_row_names = T, show_column_names = F, show_row_dend = F, show_column_dend = F, row_names_side = "left",
  row_names_gp = gpar(fontsize = 8),
  column_km = 4,
  col = col_fun, border = T, row_split = ifelse(grepl("chrM", colnames(boo)), "mito", "auto"), top_annotation = ha
)

variants.CLL1 <- c("chrM:6414:T/C", "chrM:3541:T/C", "chrM:9633:A/G", "chrM:302:AC/A")
variants.CLL2 <- c(
  "chrM:12706:T/C", "chrM:3614:C/T", "chrM:9541:C/T", "chrM:2708:G/A", "chrM:6000:T/C", "chrM:150:T/C", "chrM:11720:A/G",
  "chrM:16250:C/T"
)

df <- data.frame(
  variants.CLL1 = rowMeans(vafs.1[, variants.CLL1]),
  variants.CLL2 = rowMeans(vafs.1[, variants.CLL2])
)

ggplot(df, aes(x = variants.CLL1, y = variants.CLL2)) +
  geom_point() +
  geom_vline(xintercept = c(20, 80)) +
  geom_hline(yintercept = c(20, 80))

cells.CLL1 <- rownames(df)[which(df$variants.CLL1 > 80 & df$variants.CLL2 < 20)]
cells.CLL2 <- rownames(df)[which(df$variants.CLL2 > 80 & df$variants.CLL1 < 20)]

annotation <- data.frame(
  bc = rownames(vafs.1),
  annotation = NA
)
rownames(annotation) <- annotation$bc
annotation[cells.CLL1, "annotation"] <- "CLL1.recipient"
annotation[cells.CLL2, "annotation"] <- "CLL2.recipient"
annotation$sample <- stringr::str_split_fixed(annotation$annotation, pattern = "\\.", n = 2)[, 1]
annotation$sample[which(annotation$sample == "CLL1")] <- "CLL1_1"
annotation$sample[which(annotation$sample == "CLL2")] <- "CLL2_1"
annotation$individual <- stringr::str_split_fixed(annotation$annotation, pattern = "\\.", n = 2)[, 2]

write.csv2(annotation, file = "./data/coevolution/CLL4/20220506_deconvolution.csv", quote = F)
