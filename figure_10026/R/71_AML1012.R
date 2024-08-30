library(ArchR)
library(Rsamtools)
library(gplots)

source("./R/20200708_filter_doublets_archr_fix.R")

addArchRGenome("hg38")

AML.input.files <- createArrowFiles(
  inputFiles = c(
    "./data/10026/AML1012_12/fragments.tsv.gz",
    "./data/10026/AML1012_34/fragments.tsv.gz",
    "./data/10026/Pool148_17/fragments.tsv.gz"
  ),
  sampleNames = c("AML1012_12", "AML1012_34", "AML1012_A"),
  minTSS = 4, minFrags = 1000, addTileMat = T, addGeneScoreMat = T, threads = 1
)

doubScores <- addDoubletScores(input = AML.input.files, k = 10, knnMethod = "UMAP", LSIMethod = 1, threads = 1)

# create ArchR Project
AML.1012 <- ArchRProject(ArrowFiles = AML.input.files, outputDirectory = "./data/10026/AML.1012.all", copyArrows = F)

# filter doublets
AML.1012 <- filterDoublets(AML.1012)

# basic clustering to identify CLL cells
AML.1012 <- addIterativeLSI(
  ArchRProj = AML.1012, useMatrix = "TileMatrix", name = "IterativeLSI",
  iterations = 4, clusterParams = list(resolution = 1, sampleCells = 10000, n.start = 10),
  varFeatures = 25000, dimsToUse = 1:30
)
AML.1012 <- addClusters(input = AML.1012, reducedDims = "IterativeLSI", method = "Seurat", name = "Clusters", resolution = 1)
set.seed(1987)
AML.1012 <- addUMAP(ArchRProj = AML.1012, reducedDims = "IterativeLSI", name = "UMAP", nNeighbors = 30, minDist = 0.5, metric = "cosine", force = T)
AML.1012 <- addImputeWeights(AML.1012)

saveArchRProject(AML.1012)

### subset cells with available mitochondrial barcoding
combined.mutation.frequencies <- readRDS("./data/10026/mtDNA/20240516_AML1012_combined_mutation_frequencies.rds")
AML.1012.mito <- subsetArchRProject(
  ArchRProj = AML.1012,
  cells = AML.1012$cellNames[which(AML.1012$cellNames %in% colnames(combined.mutation.frequencies))],
  outputDirectory = "./data/10026/AML.1012.all.mito/", threads = 12
)
AML.1012.mito <- addClusters(input = AML.1012.mito, reducedDims = "IterativeLSI", method = "Seurat", name = "Clusters", resolution = 0.3, force = T)
AML.1012.mito <- addUMAP(ArchRProj = AML.1012.mito, reducedDims = "IterativeLSI", name = "UMAP", nNeighbors = 30, minDist = 0.5, metric = "cosine", force = T)
AML.1012.mito <- addImputeWeights(AML.1012.mito)
saveArchRProject(AML.1012.mito)

cells.mat <- read.csv2(file = "./data/10026/20210611_cells_by_sample.csv", sep = " ") %>% as.data.frame()
rownames(cells.mat) <- cells.mat$barcode
cells.mat <- cells.mat[intersect(cells.mat$barcode, AML.1012.mito$cellNames), ]
AML.1012.mito$Sample.hash <- cells.mat[AML.1012.mito$cellNames, "sample"]
AML.1012.mito$Sample.hash[which(AML.1012.mito$Sample == "AML1012_A")] <- "AML1012_A"
saveArchRProject(AML.1012.mito)

# plot mitochondrial barcodes
p <- plotEmbedding(AML.1012.mito)
plotPDF(p, name = "AML.1012.mito.umap", ArchRProj = AML.1012.mito, width = 4, height = 4, addDOC = F)
for (mutation in rownames(combined.mutation.frequencies)) {
  AML.1012.mito$vaf <- unlist(combined.mutation.frequencies[mutation, AML.1012.mito$cellNames])
  AML.1012.mito$vaf[which(AML.1012.mito$vaf > 0.1)] <- 0.1
  p <- plotEmbedding(AML.1012.mito, name = "vaf", pal = c("grey", "red"), plotAs = "points", na.rm = T) +
    ggtitle(mutation) + theme(plot.title = element_text(hjust = 0.5))
  plotPDF(p, name = paste0("AML.1012.mito.", mutation), ArchRProj = AML.1012.mito, addDOC = F, width = 4, height = 4)
}
