library(dplyr)
library(Seurat)

AML1010.1.recipient <- data.table::fread("./data/AML_coevolution/AML2/AML1010.recipient.barcodes.csv", header = F) %>% as.data.frame()
AML1010.1.donor <- data.table::fread("./data/AML_coevolution/AML2/AML1010.donor.barcodes.csv", header = F) %>% as.data.frame()
AML1010.5.recipient <- data.table::fread("./data/AML_coevolution/AML1/AML1010.recipient.barcodes.csv", header = F) %>% as.data.frame()
AML1010.5.donor <- data.table::fread("./data/AML_coevolution/AML1/AML1010.donor.barcodes.csv", header = F) %>% as.data.frame()

AML1.df <- data.table::fread("./data/AML_coevolution/AML1/AF.csv") %>% as.data.frame()
rownames(AML1.df) <- AML1.df$Barcode
AML1.df <- AML1.df[, seq(3, ncol(AML1.df))]
AML1010.5.recipient.df <- AML1.df[AML1010.5.recipient$V1, ]
AML1010.5.donor.df <- AML1.df[AML1010.5.donor$V1, ]
rownames(AML1010.5.recipient.df) <- paste0("AML1#", rownames(AML1010.5.recipient.df))
rownames(AML1010.5.donor.df) <- paste0("AML1#", rownames(AML1010.5.donor.df))
write.csv2(AML1010.5.recipient.df, file = "./data/AML_coevolution/AML1010/AML1010.5.recipient.csv", quote = F)
write.csv2(AML1010.5.donor.df, file = "./data/AML_coevolution/AML1010/AML1010.5.donor.csv", quote = F)

AML2.df <- data.table::fread("./data/AML_coevolution/AML2/AF.csv") %>% as.data.frame()
rownames(AML2.df) <- AML2.df$Barcode
AML2.df <- AML2.df[, seq(3, ncol(AML2.df))]
AML1010.1.recipient.df <- AML2.df[AML1010.1.recipient$V1, ]
AML1010.1.donor.df <- AML2.df[AML1010.1.donor$V1, ]
rownames(AML1010.1.recipient.df) <- paste0("AML2#", rownames(AML1010.1.recipient.df))
rownames(AML1010.1.donor.df) <- paste0("AML2#", rownames(AML1010.1.donor.df))
write.csv2(AML1010.1.recipient.df, file = "./data/AML_coevolution/AML1010/AML1010.1.recipient.csv", quote = F)
write.csv2(AML1010.1.donor.df, file = "./data/AML_coevolution/AML1010/AML1010.1.donor.csv", quote = F)
