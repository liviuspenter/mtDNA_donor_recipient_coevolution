library(ComplexHeatmap)
library(dplyr)
library(ggplot2)

AML1.cluster <- data.table::fread("./data/AML_coevolution/AML1/clustering.csv") %>%
  as.data.frame()
rownames(AML1.cluster) <- AML1.cluster$Barcode

AML1.df <- data.table::fread("./data/AML_coevolution/AML1/AF.csv") %>% as.data.frame()
rownames(AML1.df) <- AML1.df$Barcode

mito.cols <- grepl("chrM", colnames(AML1.df))

# manually identify maternal mtDNA variants
Heatmap(t(AML1.df[sample(seq(1, nrow(AML1.df)), 1000, replace = F), mito.cols]),
  show_column_names = F,
  cluster_columns = T, row_names_gp = gpar(fontsize = 8)
)

# identify AML1010 and AML1012
AML1010.recipient <- gtools::mixedsort(c(
  "10084T>C", "12612A>G", "15452C>A", "4216T>C", "55T>C", "56A>G",
  "185G>A", "228G>A", "16069C>T", "14798T>C", "16126T>C", "3010G>A",
  "11251A>G"
))

AML1010.donor <- gtools::mixedsort(c("9070T>G", "15530T>C", "15693T>C", "4811A>G", "16179C>T", "1811A>G", "150C>T", "195T>C"))

AML1012 <- gtools::mixedsort(c(
  "16291C>T", "16256C>T", "16270C>T", "15218A>G", "13617T>C", "3197T>C", "15511T>C",
  "9477G>A", "14793A>G", "16399A>G", "12582A>G"
))

AML1.df.mito <- AML1.df[, which(grepl("chrM", colnames(AML1.df)))]
colnames(AML1.df.mito) <- gsub(colnames(AML1.df.mito), pattern = "chrM:", replacement = "")
colnames(AML1.df.mito) <- gsub(colnames(AML1.df.mito), pattern = ":", replacement = "")
colnames(AML1.df.mito) <- gsub(colnames(AML1.df.mito), pattern = "/", replacement = ">")

AML1.df.mito$AML1010.donor <- rowMeans(AML1.df.mito[, intersect(AML1010.donor, colnames(AML1.df.mito))])
AML1.df.mito$AML1010.recipient <- rowMeans(AML1.df.mito[, intersect(AML1010.recipient, colnames(AML1.df.mito))])
AML1.df.mito$AML1012 <- rowMeans(AML1.df.mito[, intersect(AML1012, colnames(AML1.df.mito))])

AML1.df.mito$cluster <- AML1.cluster[rownames(AML1.df.mito), 3]

# compare mito deconvolution to clustering from Mission Bio
ggplot(AML1.df.mito, aes(x = AML1012, y = AML1010.recipient)) +
  geom_point(aes(color = cluster))
ggplot(AML1.df.mito, aes(x = AML1012, y = AML1010.donor)) +
  geom_point(aes(color = cluster))

AML1.df.mito$individual <- "none"
AML1.df.mito$individual[(which(AML1.df.mito$AML1012 < 20 &
  AML1.df.mito$AML1010.donor > 80 &
  AML1.df.mito$AML1010.recipient < 20))] <- "AML1010.donor"
AML1.df.mito$individual[(which(AML1.df.mito$AML1012 < 20 &
  AML1.df.mito$AML1010.donor < 20 &
  AML1.df.mito$AML1010.recipient > 80))] <- "AML1010.recipient"
AML1.df.mito$individual[(which(AML1.df.mito$AML1012 > 80 &
  AML1.df.mito$AML1010.donor < 20 &
  AML1.df.mito$AML1010.recipient < 20))] <- "AML1012"

# sanity check
AML1.df.mito$cluster.individual <- paste0(AML1.df.mito$cluster, ".", AML1.df.mito$individual)
table(AML1.df.mito$cluster.individual)

AML1012.cells <- rownames(AML1.df.mito)[which(AML1.df.mito$individual == "AML1012")]

AML1012.df <- AML1.df[AML1012.cells, seq(3, ncol(AML1.df))]
AML1012.df <- AML1012.df[, which(colMeans(AML1012.df) > 5)]

Heatmap(t(AML1012.df[seq(1, 1000), ]),
  show_column_names = F, use_raster = F,
  cluster_columns = T, row_names_gp = gpar(fontsize = 8)
)

AML1012.df$cluster <- AML1.cluster[rownames(AML1012.df), 3]

AML1012.donor <- c("chr3:39119965:G/A", "WDR48:chr3:39119973:G/A")
AML1012.recipient <- c("chr12:57490018:G/A", "chr11:3228322:C/T")

AML1012.df$AML1012.donor <- rowMeans(AML1012.df[, AML1012.donor])
AML1012.df$AML1012.recipient <- rowMeans(AML1012.df[, AML1012.recipient])

ggplot(AML1012.df, aes(x = AML1012.donor, y = AML1012.recipient)) +
  geom_point(aes(color = cluster))

AML1012.df$individual <- "none"
AML1012.df$individual[which(AML1012.df$AML1012.donor > 80 & AML1012.df$AML1012.recipient < 80)] <- "AML1012.donor"
AML1012.df$individual[which(AML1012.df$AML1012.donor < 80 & AML1012.df$AML1012.recipient > 80)] <- "AML1012.recipient"

write.table(
  file = "./data/AML_coevolution/AML1/AML1010.recipient.barcodes.csv",
  rownames(AML1.df.mito)[which(AML1.df.mito$individual == "AML1010.recipient")],
  quote = F, row.names = F, col.names = F
)
write.table(
  file = "./data/AML_coevolution/AML1/AML1010.donor.barcodes.csv",
  rownames(AML1.df.mito)[which(AML1.df.mito$individual == "AML1010.donor")],
  quote = F, row.names = F, col.names = F
)
write.table(
  file = "./data/AML_coevolution/AML1/AML1012.recipient.barcodes.csv",
  rownames(AML1012.df)[which(AML1012.df$individual == "AML1012.recipient")],
  quote = F, row.names = F, col.names = F
)
write.table(
  file = "./data/AML_coevolution/AML1/AML1012.donor.barcodes.csv",
  rownames(AML1012.df)[which(AML1012.df$individual == "AML1012.donor")],
  quote = F, row.names = F, col.names = F
)

write.table(
  file = "./data/AML_coevolution/AML1/AML1.barcodes.csv",
  c(
    rownames(AML1.df.mito)[which(AML1.df.mito$individual == "AML1010.recipient")],
    rownames(AML1.df.mito)[which(AML1.df.mito$individual == "AML1010.donor")],
    rownames(AML1012.df)[which(AML1012.df$individual == "AML1012.recipient")],
    rownames(AML1012.df)[which(AML1012.df$individual == "AML1012.donor")]
  ),
  quote = F, row.names = F, col.names = F
)
