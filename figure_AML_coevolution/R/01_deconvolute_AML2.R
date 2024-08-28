library(ComplexHeatmap)
library(dplyr)
library(ggplot2)

AML2.cluster <- data.table::fread("./data/AML_coevolution/AML2/clustering.csv") %>%
  as.data.frame()
rownames(AML2.cluster) <- AML2.cluster$Barcode

AML2.df <- data.table::fread("./data/AML_coevolution/AML2/AF.csv") %>% as.data.frame()
rownames(AML2.df) <- AML2.df$Barcode

mito.cols <- grepl("chrM", colnames(AML2.df))

Heatmap(t(AML2.df[sample(seq(1, nrow(AML2.df)), 1000, replace = F), mito.cols]),
  show_column_names = F,
  cluster_columns = T, row_names_gp = gpar(fontsize = 8)
)

# identify AML1010 and AML1026
AML1010.recipient <- gtools::mixedsort(c(
  "3010G>A", "16069C>T", "14798T>C", "185G>A", "228G>A", "55T>C", "56A>G",
  "12612A>G", "10084T>C"
))

AML1010.donor <- gtools::mixedsort(c(
  "1811A>G", "15693T>C", "195T>C", "9070T>G", "11467A>G", "4811A>G", "16179C>T",
  "15530T>C", "12308A>G"
))

AML1026.donor <- gtools::mixedsort(c(
  "16294C>T", "16296C>T", "16153G>A", "13368G>A", "8269G>A", "14905G>A", "11812A>G", "709G>A",
  "8697G>A", "15928G>A", "14233A>G", "1888G>A"
))

AML1026.recipient <- gtools::mixedsort(c("16304T>C", "9722T>C", "93A>G", "8433T>C", "4336T>C"))

AML2.df.mito <- AML2.df[, which(grepl("chrM", colnames(AML2.df)))]
colnames(AML2.df.mito) <- gsub(colnames(AML2.df.mito), pattern = "chrM:", replacement = "")
colnames(AML2.df.mito) <- gsub(colnames(AML2.df.mito), pattern = ":", replacement = "")
colnames(AML2.df.mito) <- gsub(colnames(AML2.df.mito), pattern = "/", replacement = ">")

AML2.df.mito$AML1010.donor <- rowMeans(AML2.df.mito[, intersect(AML1010.donor, colnames(AML2.df.mito))])
AML2.df.mito$AML1010.recipient <- rowMeans(AML2.df.mito[, intersect(AML1010.recipient, colnames(AML2.df.mito))])
AML2.df.mito$AML1026.donor <- rowMeans(AML2.df.mito[, intersect(AML1026.donor, colnames(AML2.df.mito))])
AML2.df.mito$AML1026.recipient <- rowMeans(AML2.df.mito[, intersect(AML1026.recipient, colnames(AML2.df.mito))])

AML2.df.mito$cluster <- AML2.cluster[rownames(AML2.df.mito), 3]

# compare mito deconvolution to clustering from Mission Bio
ggplot(AML2.df.mito, aes(x = AML1026.donor, y = AML1026.recipient)) +
  geom_point(aes(color = cluster))
ggplot(AML2.df.mito, aes(x = AML1010.donor, y = AML1010.recipient)) +
  geom_point(aes(color = cluster))

AML2.df.mito$individual <- "none"
AML2.df.mito$individual[(which(AML2.df.mito$AML1010.recipient > 80 &
  AML2.df.mito$AML1010.donor < 20 &
  AML2.df.mito$AML1026.donor < 20 &
  AML2.df.mito$AML1026.recipient < 20))] <- "AML1010.recipient"
AML2.df.mito$individual[(which(AML2.df.mito$AML1010.recipient < 20 &
  AML2.df.mito$AML1010.donor > 80 &
  AML2.df.mito$AML1026.donor < 20 &
  AML2.df.mito$AML1026.recipient < 20))] <- "AML1010.donor"
AML2.df.mito$individual[(which(AML2.df.mito$AML1010.recipient < 20 &
  AML2.df.mito$AML1010.donor < 20 &
  AML2.df.mito$AML1026.donor > 80 &
  AML2.df.mito$AML1026.recipient < 20))] <- "AML1026.donor"
AML2.df.mito$individual[(which(AML2.df.mito$AML1010.recipient < 20 &
  AML2.df.mito$AML1010.donor < 20 &
  AML2.df.mito$AML1026.donor < 20 &
  AML2.df.mito$AML1026.recipient > 80))] <- "AML1026.recipient"
# sanity check
AML2.df.mito$cluster.individual <- paste0(AML2.df.mito$cluster, ".", AML2.df.mito$individual)
table(AML2.df.mito$cluster.individual)

write.table(
  file = "./data/AML_coevolution/AML2/AML1010.recipient.barcodes.csv",
  rownames(AML2.df.mito)[which(AML2.df.mito$individual == "AML1010.recipient")],
  quote = F, row.names = F, col.names = F
)
write.table(
  file = "./data/AML_coevolution/AML2/AML1010.donor.barcodes.csv",
  rownames(AML2.df.mito)[which(AML2.df.mito$individual == "AML1010.donor")],
  quote = F, row.names = F, col.names = F
)
write.table(
  file = "./data/AML_coevolution/AML2/AML1026.recipient.barcodes.csv",
  rownames(AML2.df.mito)[which(AML2.df.mito$individual == "AML1026.recipient")],
  quote = F, row.names = F, col.names = F
)
write.table(
  file = "./data/AML_coevolution/AML2/AML1026.donor.barcodes.csv",
  rownames(AML2.df.mito)[which(AML2.df.mito$individual == "AML1026.donor")],
  quote = F, row.names = F, col.names = F
)

write.table(
  file = "./data/AML_coevolution/AML2/AML2.barcodes.csv",
  c(
    rownames(AML2.df.mito)[which(AML2.df.mito$individual == "AML1010.recipient")],
    rownames(AML2.df.mito)[which(AML2.df.mito$individual == "AML1010.donor")],
    rownames(AML2.df.mito)[which(AML2.df.mito$individual == "AML1026.recipient")],
    rownames(AML2.df.mito)[which(AML2.df.mito$individual == "AML1026.donor")]
  ),
  quote = F, row.names = F, col.names = F
)
