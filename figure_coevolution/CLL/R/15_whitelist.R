# create whitelist of mutations to be called by Tapestri software for each sample

library(ComplexHeatmap)
library(dplyr)
library(ggplot2)
library(Seurat)

vafs.1 <- as.data.frame(data.table::fread("./data/coevolution/CLL1/AF.csv"))
vafs.2 <- as.data.frame(data.table::fread("./data/coevolution/CLL2/AF.csv"))
vafs.3 <- as.data.frame(data.table::fread("./data/coevolution/CLL3/AF.csv"))
vafs.4 <- as.data.frame(data.table::fread("./data/coevolution/CLL4/AF.csv"))

# Patient CLL1 in samples 3 and 4
vafs <- dplyr::bind_rows(vafs.3, vafs.4)
mutations.1.bed <- data.frame(
  chr = c(
    stringr::str_split_fixed(colnames(vafs)[which(startsWith(colnames(vafs), "chr"))], pattern = "\\:", n = 3)[, 1],
    stringr::str_split_fixed(colnames(vafs)[which(!startsWith(colnames(vafs), "chr"))], pattern = "\\:", n = 4)[, 2]
  ),
  start = c(
    stringr::str_split_fixed(colnames(vafs)[which(startsWith(colnames(vafs), "chr"))], pattern = "\\:", n = 3)[, 2],
    stringr::str_split_fixed(colnames(vafs)[which(!startsWith(colnames(vafs), "chr"))], pattern = "\\:", n = 4)[, 3]
  ),
  stop = c(
    stringr::str_split_fixed(colnames(vafs)[which(startsWith(colnames(vafs), "chr"))], pattern = "\\:", n = 3)[, 2],
    stringr::str_split_fixed(colnames(vafs)[which(!startsWith(colnames(vafs), "chr"))], pattern = "\\:", n = 4)[, 3]
  )
)
mutations.1.bed <- mutations.1.bed[-which(mutations.1.bed$chr == ""), ]

# Patient CLL2 in samples 1 and 4
vafs <- dplyr::bind_rows(vafs.1, vafs.4)
mutations.2.bed <- data.frame(
  chr = c(
    stringr::str_split_fixed(colnames(vafs)[which(startsWith(colnames(vafs), "chr"))], pattern = "\\:", n = 3)[, 1],
    stringr::str_split_fixed(colnames(vafs)[which(!startsWith(colnames(vafs), "chr"))], pattern = "\\:", n = 4)[, 2]
  ),
  start = c(
    stringr::str_split_fixed(colnames(vafs)[which(startsWith(colnames(vafs), "chr"))], pattern = "\\:", n = 3)[, 2],
    stringr::str_split_fixed(colnames(vafs)[which(!startsWith(colnames(vafs), "chr"))], pattern = "\\:", n = 4)[, 3]
  ),
  stop = c(
    stringr::str_split_fixed(colnames(vafs)[which(startsWith(colnames(vafs), "chr"))], pattern = "\\:", n = 3)[, 2],
    stringr::str_split_fixed(colnames(vafs)[which(!startsWith(colnames(vafs), "chr"))], pattern = "\\:", n = 4)[, 3]
  )
)
mutations.2.bed <- mutations.2.bed[-which(mutations.2.bed$chr == ""), ]

# Patient CLL3 in samples 2 and 3
vafs <- dplyr::bind_rows(vafs.2, vafs.3)
mutations.3.bed <- data.frame(
  chr = c(
    stringr::str_split_fixed(colnames(vafs)[which(startsWith(colnames(vafs), "chr"))], pattern = "\\:", n = 3)[, 1],
    stringr::str_split_fixed(colnames(vafs)[which(!startsWith(colnames(vafs), "chr"))], pattern = "\\:", n = 4)[, 2]
  ),
  start = c(
    stringr::str_split_fixed(colnames(vafs)[which(startsWith(colnames(vafs), "chr"))], pattern = "\\:", n = 3)[, 2],
    stringr::str_split_fixed(colnames(vafs)[which(!startsWith(colnames(vafs), "chr"))], pattern = "\\:", n = 4)[, 3]
  ),
  stop = c(
    stringr::str_split_fixed(colnames(vafs)[which(startsWith(colnames(vafs), "chr"))], pattern = "\\:", n = 3)[, 2],
    stringr::str_split_fixed(colnames(vafs)[which(!startsWith(colnames(vafs), "chr"))], pattern = "\\:", n = 4)[, 3]
  )
)
mutations.3.bed <- mutations.3.bed[-which(mutations.3.bed$chr == ""), ]

# whitelist for sample 1 - CLL2
write.table(mutations.2.bed, file = "./data/coevolution/CLL1/20220508_whitelist.bed", sep = "\t", quote = F, row.names = F, col.names = F)

# whitelist for sample 2 - CLL3
write.table(mutations.3.bed, file = "./data/coevolution/CLL2/20220508_whitelist.bed", sep = "\t", quote = F, row.names = F, col.names = F)

# whitelist for sample 3 - CLL1 + CLL3
write.table(rbind(mutations.1.bed, mutations.3.bed), file = "./data/coevolution/CLL3/20220508_whitelist.bed", sep = "\t", quote = F, row.names = F, col.names = F)

# whitelist for sample 4 - CLL1 + CLL2
write.table(rbind(mutations.1.bed, mutations.2.bed), file = "./data/coevolution/CLL4/20220508_whitelist.bed", sep = "\t", quote = F, row.names = F, col.names = F)
