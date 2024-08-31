# sample deconvolution of sample CLL3
# CLL3 contains CLL4 post-HSCT (called CLL1) and CLL5 pre-FCR (called CLL3)

library(ComplexHeatmap)
library(dplyr)
library(ggplot2)
library(Seurat)

col_fun <- circlize::colorRamp2(breaks = seq(0, 100, 100 / 8), colors = BuenColors::jdb_palette(name = "solar_rojos"))

so.12 <- readRDS("./data/coevolution/objects/20220504_CLL34.rds")

vafs.1 <- as.data.frame(data.table::fread("./data/coevolution/CLL3/AF.csv"))
vafs.1$Barcode <- paste0("CLL3#", vafs.1$Barcode)
rownames(vafs.1) <- vafs.1$Barcode
vafs.1 <- vafs.1[, -c(1, 2)]

# explore SNPs
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
  column_km = 3,
  col = col_fun, border = T, row_split = ifelse(grepl("chrM", colnames(boo)), "mito", "auto"), top_annotation = ha
)

# variants for CLL4 (called CLL1)
variants.CLL1 <- c("chrM:6414:T/C", "chrM:3541:T/C", "chrM:9633:A/G", "chrM:302:AC/A")
# variants for CLL5 (called CLL3)
variants.CLL3 <- c("chrM:11720:A/G", "chrM:12706:T/C", "chrM:16300:T/C", "chrM:150:T/C", "chrM:9541:C/T", "chrM:14017:G/A")

df <- data.frame(
  variants.CLL1 = rowMeans(vafs.1[, variants.CLL1]),
  variants.CLL3 = rowMeans(vafs.1[, variants.CLL3])
)


cells.CLL1 <- rownames(df)[which(df$variants.CLL1 > 80 & df$variants.CLL3 < 20)]
cells.CLL3 <- rownames(df)[which(df$variants.CLL3 > 80 & df$variants.CLL1 < 20)]

###
# de-convolute donor and recipient CLL4 post-HSCT (CLL1)
###
boo <- vafs.1[cells.CLL1, ]

ha <- columnAnnotation(
  cluster = so.12$manual.cluster[rownames(boo)],
  col = list("cluster" = c(
    "CLL" = "orange", "CD4 T cell" = "lightblue", "CD8 T cell" = "darkblue", "Myeloid" = "darkgreen",
    "Cluster 6" = "black"
  ))
)

Heatmap(t(boo),
  show_row_names = T, show_column_names = F, show_row_dend = F, show_column_dend = F, row_names_side = "left",
  row_names_gp = gpar(fontsize = 8), use_raster = T, raster_quality = 10,
  column_km = 2,
  col = col_fun, border = T, row_split = ifelse(grepl("chrM", colnames(boo)), "mito", "auto"), top_annotation = ha
)

# donor variants
variants.CLL1.donor <- c("ESRRG:chr1:216880662:G/A", "TINAG:chr6:54191736:C/T")
# recipient variants
variants.CLL1.recipient <- c(
  "KCNA6:chr12:4920442:C/T", "SH3RF1:chr4:170043260:C/T", "ADGRL3:chr4:62862035:T/G", "WNK1:chr12:966369:G/T",
  "ZNF215:chr11:6953628:A/C", "DPCD:chr10:103354473:C/A", "SF3B1:chr2:198267371:G/T", "MGA:chr15:42019401:C/T",
  "SH3GL1:chr19:4363776:T/A"
)

df <- data.frame(
  variants.CLL1.donor = rowMeans(vafs.1[cells.CLL1, variants.CLL1.donor]),
  variants.CLL1.recipient = rowMeans(vafs.1[cells.CLL1, variants.CLL1.recipient])
)

ggplot(df, aes(x = variants.CLL1.donor, y = variants.CLL1.recipient)) +
  geom_point() +
  geom_vline(xintercept = c(5)) +
  geom_hline(yintercept = c(5))

cells.CLL1.donor <- rownames(df)[which(df$variants.CLL1.donor > 5 & df$variants.CLL1.recipient < 5)]
cells.CLL1.recipient <- rownames(df)[which(df$variants.CLL1.donor < 5 & df$variants.CLL1.recipient > 5)]

# deconvolute
annotation <- data.frame(
  bc = rownames(vafs.1),
  annotation = NA
)
rownames(annotation) <- annotation$bc
annotation[cells.CLL1.donor, "annotation"] <- "CLL1.donor"
annotation[cells.CLL1.recipient, "annotation"] <- "CLL1.recipient"
annotation[cells.CLL3, "annotation"] <- "CLL3.recipient" # CLL3 sample is pre-transplant
annotation$sample <- stringr::str_split_fixed(annotation$annotation, pattern = "\\.", n = 2)[, 1]
annotation$sample[which(annotation$sample == "CLL1")] <- "CLL1_2"
annotation$sample[which(annotation$sample == "CLL3")] <- "CLL3_1"
annotation$individual <- stringr::str_split_fixed(annotation$annotation, pattern = "\\.", n = 2)[, 2]

write.csv2(annotation, file = "./data/coevolution/CLL3/20220506_deconvolution.csv", quote = F)
