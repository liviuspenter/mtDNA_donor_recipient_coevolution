# donor-recipient deconvolution of sample CLL1

library(ComplexHeatmap)
library(dplyr)
library(ggplot2)
library(Seurat)

col_fun <- circlize::colorRamp2(breaks = seq(0, 100, 100 / 8), colors = BuenColors::jdb_palette(name = "solar_rojos"))

vafs.1 <- as.data.frame(data.table::fread("./data/coevolution/CLL1/AF.csv"))
vafs.1$Barcode <- paste0("CLL1#", vafs.1$Barcode)
rownames(vafs.1) <- vafs.1$Barcode
vafs.1 <- vafs.1[, -c(1, 2)]

# explore germline SNPs
boo <- vafs.1[sample(nrow(vafs.1), size = 1000), ]
Heatmap(t(boo),
  show_row_names = T, show_column_names = F, show_row_dend = F, show_column_dend = F, row_names_side = "left",
  row_names_gp = gpar(fontsize = 8),
  column_km = 2, use_raster = T, raster_quality = 10,
  col = col_fun, border = T, row_split = ifelse(grepl("chrM", colnames(boo)), "mito", "auto")
)

# recipient SNPs
variants.CLL2.recipient <- c(
  "ERICH3:chr1:75102115:G/C", "WASF3:chr13:27255288:G/T", "ESRRG:chr1:216880662:G/A", "DNAJB7:chr22:41257549:A/C",
  "BPIFA2:chr20:31761967:G/C", "ZNF292:chr6:87967797:CA/C", "PLEC:chr8:144994678:CCT/C", "DCBLD2:chr3:98600507:G/A"
)
# donor SNPs
variants.CLL2.donor <- c("LAMA5:chr20:60912683:T/C")

# deconvolute
df <- data.frame(
  variants.CLL2.recipient = rowMeans(vafs.1[, variants.CLL2.recipient]),
  variants.CLL2.donor = vafs.1[, variants.CLL2.donor]
)

cells.CLL2.recipient <- rownames(df)[which(df$variants.CLL2.recipient > 10 & df$variants.CLL2.donor < 10)]
cells.CLL2.donor <- rownames(df)[which(df$variants.CLL2.recipient < 10 & df$variants.CLL2.donor > 10)]

annotation <- data.frame(
  bc = rownames(vafs.1),
  annotation = NA
)
rownames(annotation) <- annotation$bc
annotation[cells.CLL2.recipient, "annotation"] <- "CLL2.recipient"
annotation[cells.CLL2.donor, "annotation"] <- "CLL2.donor"
annotation$sample <- stringr::str_split_fixed(annotation$annotation, pattern = "\\.", n = 2)[, 1]
annotation$sample[which(annotation$sample == "CLL2")] <- "CLL2_2"
annotation$individual <- stringr::str_split_fixed(annotation$annotation, pattern = "\\.", n = 2)[, 2]

write.csv2(annotation, file = "./data/coevolution/CLL1/20220506_deconvolution.csv", quote = F)
