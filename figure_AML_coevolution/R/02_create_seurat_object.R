library(dplyr)
library(ggplot2)
library(Seurat)

barcodes.1 <- data.table::fread("./data/AML_coevolution/AML1/AML1.barcodes.csv", header = F) %>% as.data.frame()
barcodes.1$V1 <- paste0("AML1#", barcodes.1$V1)
barcodes.2 <- data.table::fread("./data/AML_coevolution/AML2/AML2.barcodes.csv", header = F) %>% as.data.frame()
barcodes.2$V1 <- paste0("AML2#", barcodes.2$V1)
barcodes.3 <- data.table::fread("./data/AML_coevolution/AML3/AML3.barcodes.csv", header = F) %>% as.data.frame()
barcodes.3$V1 <- paste0("AML3#", barcodes.3$V1)

protein.1 <- as.data.frame(data.table::fread("./data/AML_coevolution/AML1/protein_reads.csv"))
protein.1$Barcode <- paste0("AML1#", protein.1$Barcode)
rownames(protein.1) <- protein.1$Barcode
protein.1 <- protein.1[, -c(1, 2)]
protein.1 <- t(protein.1)
protein.2 <- as.data.frame(data.table::fread("./data/AML_coevolution/AML2/protein_reads.csv"))
protein.2$Barcode <- paste0("AML2#", protein.2$Barcode)
rownames(protein.2) <- protein.2$Barcode
protein.2 <- protein.2[, -c(1, 2)]
protein.2 <- t(protein.2)
protein.3 <- as.data.frame(data.table::fread("./data/AML_coevolution/AML3/protein_reads.csv"))
protein.3$Barcode <- paste0("AML3#", protein.3$Barcode)
rownames(protein.3) <- protein.3$Barcode
protein.3 <- protein.3[, -c(1, 2)]
protein.3 <- t(protein.3)

so.1 <- CreateSeuratObject(counts = protein.1, assay = "ADT", project = "AML1")
so.2 <- CreateSeuratObject(counts = protein.2, assay = "ADT", project = "AML2")
so.3 <- CreateSeuratObject(counts = protein.3, assay = "ADT", project = "AML3")

so.all <- merge(so.1, list(so.2, so.3))

so.all <- subset(so.all, cells = c(barcodes.1$V1, barcodes.2$V1, barcodes.3$V1))

so.all <- NormalizeData(so.all, normalization.method = "CLR", margin = 2)
so.all <- FindVariableFeatures(so.all)
so.all <- ScaleData(so.all, features = rownames(so.all))
so.all <- RunPCA(so.all)
so.all <- RunUMAP(so.all, dims = 1:15)
so.all <- FindNeighbors(so.all)
so.all <- FindClusters(so.all, resolution = 0.2)

so.all <- JoinLayers(so.all)

pbmc.markers <- FindAllMarkers(so.all, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25) %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC)

so.all$manual.cluster <- "none"
so.all$manual.cluster[which(so.all$seurat_clusters %in% c("0", "3"))] <- "HSC"
so.all$manual.cluster[which(so.all$seurat_clusters %in% c("2"))] <- "Progenitor"
so.all$manual.cluster[which(so.all$seurat_clusters %in% c("1", "4"))] <- "Mono"
so.all$manual.cluster[which(so.all$seurat_clusters %in% c("7"))] <- "CD8+ T cell"
so.all$manual.cluster[which(so.all$seurat_clusters %in% c("6"))] <- "CD4+ T cell"
so.all$manual.cluster[which(so.all$seurat_clusters %in% c("8"))] <- "Plasma cell"
so.all$manual.cluster[which(so.all$seurat_clusters %in% c("9"))] <- "NK"
so.all$manual.cluster[which(so.all$seurat_clusters %in% c("5"))] <- "Erythroid"
so.all$manual.cluster[which(so.all$seurat_clusters %in% c("10"))] <- "B cell"


AML.colors <- c(
  "HSC" = RColorBrewer::brewer.pal(name = "YlOrRd", n = 5)[5],
  "Progenitor" = RColorBrewer::brewer.pal(name = "YlOrRd", n = 5)[4],
  "Mono" = RColorBrewer::brewer.pal(name = "YlOrRd", n = 5)[3],
  "DC" = RColorBrewer::brewer.pal(name = "YlOrRd", n = 5)[2],
  "CD4+ T cell" = "lightblue",
  "CD8+ T cell" = "darkblue",
  "NK" = "purple",
  "B cell" = "lightgreen",
  "Plasma cell" = "darkgreen"
)
so.all$manual.cluster <- factor(so.all$manual.cluster, levels = names(AML.colors))

DoHeatmap(so.all, features = pbmc.markers$gene, group.by = "manual.cluster", group.colors = AML.colors) +
  scale_fill_gradientn(colours = BuenColors::jdb_palette(name = "brewer_yes"))

saveRDS(file = "./data/AML_coevolution/objects/20240416_AML_so.rds", so.all)
