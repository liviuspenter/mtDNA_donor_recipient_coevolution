library(ComplexHeatmap)
library(dplyr)
library(ggplot2)

AML3.cluster <- data.table::fread("./data/AML_coevolution/AML3/clustering.csv") %>%
  as.data.frame()
rownames(AML3.cluster) <- AML3.cluster$Barcode

AML3.df <- data.table::fread("./data/AML_coevolution/AML3/AF.csv") %>% as.data.frame()
rownames(AML3.df) <- AML3.df$Barcode

mito.cols <- grepl("chrM", colnames(AML3.df))

Heatmap(t(AML3.df[sample(seq(1, nrow(AML3.df)), 1000, replace = F), mito.cols]),
  show_column_names = F,
  cluster_columns = T, row_names_gp = gpar(fontsize = 8)
)

# identify AML1012 and AML1026
AML1012 <- gtools::mixedsort(c(
  "14793A>G", "15218A>G", "13617T>C", "3197T>C", "11467A>G", "16399A>G", "12308A>G",
  "16256C>T", "16270C>T", "16291C>T", "15511T>C", "9477G>A", "12582A>G"
))

AML1026.donor <- gtools::mixedsort(c(
  "15549T>C", "13368G>A", "8269G>A", "709G>A", "14905G>A", "16294C>T", "16296C>T",
  "11812A>G", "15928G>A", "15452C>A", "16519T>C", "4216T>C", "8697G>A", "16126T>C",
  "16153G>A", "150C>T", "14233A>G", "1888G>A", "11251A>G"
))

AML1026.recipient <- gtools::mixedsort(c("16304T>C", "9722T>C", "93A>G", "8433T>C", "4336T>C"))

AML3.df.mito <- AML3.df[, which(grepl("chrM", colnames(AML3.df)))]
colnames(AML3.df.mito) <- gsub(colnames(AML3.df.mito), pattern = "chrM:", replacement = "")
colnames(AML3.df.mito) <- gsub(colnames(AML3.df.mito), pattern = ":", replacement = "")
colnames(AML3.df.mito) <- gsub(colnames(AML3.df.mito), pattern = "/", replacement = ">")

AML3.df.mito$AML1012 <- rowMeans(AML3.df.mito[, intersect(AML1012, colnames(AML3.df.mito))])
AML3.df.mito$AML1026.donor <- rowMeans(AML3.df.mito[, intersect(AML1026.donor, colnames(AML3.df.mito))])
AML3.df.mito$AML1026.recipient <- rowMeans(AML3.df.mito[, intersect(AML1026.recipient, colnames(AML3.df.mito))])

AML3.df.mito$cluster <- AML3.cluster[rownames(AML3.df.mito), 3]

# compare mito deconvolution to clustering from Mission Bio
ggplot(AML3.df.mito, aes(x = AML1012, y = AML1026.recipient)) +
  geom_point(aes(color = cluster))
ggplot(AML3.df.mito, aes(x = AML1012, y = AML1026.donor)) +
  geom_point(aes(color = cluster))

AML3.df.mito$individual <- "none"
AML3.df.mito$individual[(which(AML3.df.mito$AML1012 > 80 &
  AML3.df.mito$AML1026.donor < 20 &
  AML3.df.mito$AML1026.recipient < 20))] <- "AML1012"
AML3.df.mito$individual[(which(AML3.df.mito$AML1012 < 20 &
  AML3.df.mito$AML1026.donor > 80 &
  AML3.df.mito$AML1026.recipient < 20))] <- "AML1026.donor"
AML3.df.mito$individual[(which(AML3.df.mito$AML1012 < 20 &
  AML3.df.mito$AML1026.donor < 20 &
  AML3.df.mito$AML1026.recipient > 80))] <- "AML1026.recipient"
# sanity check
AML3.df.mito$cluster.individual <- paste0(AML3.df.mito$cluster, ".", AML3.df.mito$individual)
table(AML3.df.mito$cluster.individual)

AML1012.cells <- rownames(AML3.df.mito)[which(AML3.df.mito$individual == "AML1012")]

AML1012.df <- AML3.df[AML1012.cells, seq(3, ncol(AML3.df))]
AML1012.df <- AML1012.df[, which(colMeans(AML1012.df) > 5)]


AML1012.donor <- c("chr3:39119965:G/A", "WDR48:chr3:39119973:G/A")
AML1012.recipient <- c("chr12:57490018:G/A", "chr11:3228322:C/T")

AML1012.df$AML1012.donor <- rowMeans(AML1012.df[, AML1012.donor])
AML1012.df$AML1012.recipient <- rowMeans(AML1012.df[, AML1012.recipient])

ggplot(AML1012.df, aes(x = AML1012.donor, y = AML1012.recipient)) +
  geom_point()

AML1012.df$individual <- "none"
AML1012.df$individual[which(AML1012.df$AML1012.donor > 80 & AML1012.df$AML1012.recipient < 80)] <- "AML1012.donor"
AML1012.df$individual[which(AML1012.df$AML1012.donor < 80 & AML1012.df$AML1012.recipient > 80)] <- "AML1012.recipient"


write.table(
  file = "./data/AML_coevolution/AML3/AML1012.recipient.barcodes.csv",
  rownames(AML1012.df)[which(AML1012.df$individual == "AML1012.recipient")],
  quote = F, row.names = F, col.names = F
)
write.table(
  file = "./data/AML_coevolution/AML3/AML1012.donor.barcodes.csv",
  rownames(AML1012.df)[which(AML1012.df$individual == "AML1012.donor")],
  quote = F, row.names = F, col.names = F
)
write.table(
  file = "./data/AML_coevolution/AML3/AML1026.recipient.barcodes.csv",
  rownames(AML3.df.mito)[which(AML3.df.mito$individual == "AML1026.recipient")],
  quote = F, row.names = F, col.names = F
)
write.table(
  file = "./data/AML_coevolution/AML3/AML1026.donor.barcodes.csv",
  rownames(AML3.df.mito)[which(AML3.df.mito$individual == "AML1026.donor")],
  quote = F, row.names = F, col.names = F
)



write.table(
  file = "./data/AML_coevolution/AML3/AML3.barcodes.csv",
  c(
    rownames(AML1012.df)[which(AML1012.df$individual == "AML1012.recipient")],
    rownames(AML1012.df)[which(AML1012.df$individual == "AML1012.donor")],
    rownames(AML3.df.mito)[which(AML3.df.mito$individual == "AML1026.recipient")],
    rownames(AML3.df.mito)[which(AML3.df.mito$individual == "AML1026.donor")]
  ),
  quote = F, row.names = F, col.names = F
)
