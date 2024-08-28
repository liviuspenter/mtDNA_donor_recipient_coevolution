library(dplyr)
library(Seurat)

AML1026.1.recipient <- data.table::fread("./data/AML_coevolution/AML2/AML1026.recipient.barcodes.csv", header = F) %>% as.data.frame()
AML1026.1.donor <- data.table::fread("./data/AML_coevolution/AML2/AML1026.donor.barcodes.csv", header = F) %>% as.data.frame()
AML1026.3.recipient <- data.table::fread("./data/AML_coevolution/AML3/AML1026.recipient.barcodes.csv", header = F) %>% as.data.frame()
AML1026.3.donor <- data.table::fread("./data/AML_coevolution/AML3/AML1026.donor.barcodes.csv", header = F) %>% as.data.frame()

AML2.df <- data.table::fread("./data/AML_coevolution/AML2/AF.csv") %>% as.data.frame()
rownames(AML2.df) <- AML2.df$Barcode
AML2.df <- AML2.df[, seq(3, ncol(AML2.df))]
AML1026.1.recipient.df <- AML2.df[AML1026.1.recipient$V1, ]
AML1026.1.donor.df <- AML2.df[AML1026.1.donor$V1, ]
rownames(AML1026.1.recipient.df) <- paste0("AML2#", rownames(AML1026.1.recipient.df))
rownames(AML1026.1.donor.df) <- paste0("AML2#", rownames(AML1026.1.donor.df))
write.csv2(AML1026.1.recipient.df, file = "./data/AML_coevolution/AML1026/AML1026.1.recipient.csv", quote = F)
write.csv2(AML1026.1.donor.df, file = "./data/AML_coevolution/AML1026/AML1026.1.donor.csv", quote = F)

AML3.df <- data.table::fread("./data/AML_coevolution/AML3/AF.csv") %>% as.data.frame()
rownames(AML3.df) <- AML3.df$Barcode
AML3.df <- AML3.df[, seq(3, ncol(AML3.df))]
AML1026.3.recipient.df <- AML3.df[AML1026.3.recipient$V1, ]
AML1026.3.donor.df <- AML3.df[AML1026.3.donor$V1, ]
rownames(AML1026.3.recipient.df) <- paste0("AML3#", rownames(AML1026.3.recipient.df))
rownames(AML1026.3.donor.df) <- paste0("AML3#", rownames(AML1026.3.donor.df))
write.csv2(AML1026.3.recipient.df, file = "./data/AML_coevolution/AML1026/AML1026.3.recipient.csv", quote = F)
write.csv2(AML1026.3.donor.df, file = "./data/AML_coevolution/AML1026/AML1026.3.donor.csv", quote = F)
