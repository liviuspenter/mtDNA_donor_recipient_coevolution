library(dplyr)
library(Seurat)

AML1012.1.recipient <- data.table::fread("./data/AML_coevolution/AML1/AML1012.recipient.barcodes.csv", header = F) %>% as.data.frame()
AML1012.1.donor <- data.table::fread("./data/AML_coevolution/AML1/AML1012.donor.barcodes.csv", header = F) %>% as.data.frame()
AML1012.4.recipient <- data.table::fread("./data/AML_coevolution/AML3/AML1012.recipient.barcodes.csv", header = F) %>% as.data.frame()
AML1012.4.donor <- data.table::fread("./data/AML_coevolution/AML3/AML1012.donor.barcodes.csv", header = F) %>% as.data.frame()

AML1.df <- data.table::fread("./data/AML_coevolution/AML1/AF.csv") %>% as.data.frame()
rownames(AML1.df) <- AML1.df$Barcode
AML1.df <- AML1.df[, seq(3, ncol(AML1.df))]
AML1012.1.recipient.df <- AML1.df[AML1012.1.recipient$V1, ]
AML1012.1.donor.df <- AML1.df[AML1012.1.donor$V1, ]
rownames(AML1012.1.recipient.df) <- paste0("AML1#", rownames(AML1012.1.recipient.df))
rownames(AML1012.1.donor.df) <- paste0("AML1#", rownames(AML1012.1.donor.df))
write.csv2(AML1012.1.recipient.df, file = "./data/AML_coevolution/AML1012/AML1012.1.recipient.csv", quote = F)
write.csv2(AML1012.1.donor.df, file = "./data/AML_coevolution/AML1012/AML1012.1.donor.csv", quote = F)

AML3.df <- data.table::fread("./data/AML_coevolution/AML3/AF.csv") %>% as.data.frame()
rownames(AML3.df) <- AML3.df$Barcode
AML3.df <- AML3.df[, seq(3, ncol(AML3.df))]
AML1012.4.recipient.df <- AML3.df[AML1012.4.recipient$V1, ]
AML1012.4.donor.df <- AML3.df[AML1012.4.donor$V1, ]
rownames(AML1012.4.recipient.df) <- paste0("AML3#", rownames(AML1012.4.recipient.df))
rownames(AML1012.4.donor.df) <- paste0("AML3#", rownames(AML1012.4.donor.df))
write.csv2(AML1012.4.recipient.df, file = "./data/AML_coevolution/AML1012/AML1012.4.recipient.csv", quote = F)
write.csv2(AML1012.4.donor.df, file = "./data/AML_coevolution/AML1012/AML1012.4.donor.csv", quote = F)
