library(ArchR)
library(Rsamtools)
library(gplots)

source("./R/20200708_filter_doublets_archr_fix.R")

addArchRGenome("hg38")

AML.input.files <- createArrowFiles(
  inputFiles = c(
    "./data/10026/Pool148_13/fragments.tsv.gz",
    "./data/10026/Pool148_14/fragments.tsv.gz",
    "./data/10026/Pool148_15/fragments.tsv.gz",
    "./data/10026/Pool148_16/fragments.tsv.gz"
  ),
  sampleNames = c("AML1011_A", "AML1011_B", "AML1011_1", "AML1011_2"),
  minTSS = 4, minFrags = 1000, addTileMat = T, addGeneScoreMat = T, threads = 1
)

doubScores <- addDoubletScores(input = AML.input.files, k = 10, knnMethod = "UMAP", LSIMethod = 1, threads = 1)

# create ArchR Project
AML.1011 <- ArchRProject(ArrowFiles = AML.input.files, outputDirectory = "./AML.1011", copyArrows = F)

# filter doublets
AML.1011 <- filterDoublets(AML.1011)

# basic clustering to identify CLL cells
AML.1011 <- addIterativeLSI(
  ArchRProj = AML.1011, useMatrix = "TileMatrix", name = "IterativeLSI",
  iterations = 4, clusterParams = list(resolution = 1, sampleCells = 10000, n.start = 10),
  varFeatures = 25000, dimsToUse = 1:30
)
AML.1011 <- addClusters(input = AML.1011, reducedDims = "IterativeLSI", method = "Seurat", name = "Clusters", resolution = 1)
set.seed(1987)
AML.1011 <- addUMAP(ArchRProj = AML.1011, reducedDims = "IterativeLSI", name = "UMAP", nNeighbors = 30, minDist = 0.5, metric = "cosine", force = T)
AML.1011 <- addImputeWeights(AML.1011)

saveArchRProject(AML.1011)

### subset cells with available mitochondrial barcoding
combined.mutation.frequencies <- readRDS("./data/10026/mtDNA/20240516_AML1011_combined_mutation_frequencies.rds")
AML.1011.mito <- subsetArchRProject(
  ArchRProj = AML.1011,
  cells = AML.1011$cellNames[which(AML.1011$cellNames %in% colnames(combined.mutation.frequencies))],
  outputDirectory = "./data/10026/AML.1011.mito/", threads = 12
)
AML.1011.mito <- addClusters(input = AML.1011.mito, reducedDims = "IterativeLSI", method = "Seurat", name = "Clusters", resolution = 0.3, force = T)
AML.1011.mito <- addUMAP(ArchRProj = AML.1011.mito, reducedDims = "IterativeLSI", name = "UMAP", nNeighbors = 30, minDist = 0.5, metric = "cosine", force = T)
AML.1011.mito <- addImputeWeights(AML.1011.mito)
saveArchRProject(AML.1011.mito)

# plot mitochondrial barcodes
p <- plotEmbedding(AML.1011.mito)
plotPDF(p, name = "AML.1011.mito.umap", ArchRProj = AML.1011.mito, width = 4, height = 4, addDOC = F)
for (mutation in rownames(combined.mutation.frequencies)) {
  AML.1011.mito$vaf <- unlist(combined.mutation.frequencies[mutation, AML.1011.mito$cellNames])
  AML.1011.mito$vaf[which(AML.1011.mito$vaf > 0.1)] <- 0.1
  p <- plotEmbedding(AML.1011.mito, name = "vaf", pal = c("grey", "red"), plotAs = "points", na.rm = T) +
    ggtitle(mutation) + theme(plot.title = element_text(hjust = 0.5))
  plotPDF(p, name = paste0("AML.1011.mito.", mutation), ArchRProj = AML.1011.mito, addDOC = F, width = 4, height = 4)
}
