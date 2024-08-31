library(ArchR)
library(ggplot2)
library(parallel)
library(Seurat)
library(Signac)

# create ArchR object

addArchRGenome("hg38")

input.files <- createArrowFiles(
  inputFiles = c(
    "./data/Pool148_13/fragments.tsv.gz",
    "./data/Pool148_17/fragments.tsv.gz"
  ),
  sampleNames = c("AML1011_A", "AML1012_A"),
  minTSS = 4, minFrags = 1000, addTileMat = T, addGeneScoreMat = T
)
doubScores <- addDoubletScores(input = input.files, k = 10, knnMethod = "UMAP", LSIMethod = 1, threads = 1)

AML.mix <- ArchRProject(ArrowFiles = input.files, outputDirectory = "./data/mixing/AML.mix", copyArrows = F)
AML.mix <- filterDoublets(AML.mix)

AML.mix <- addIterativeLSI(
  ArchRProj = AML.mix, useMatrix = "TileMatrix", name = "IterativeLSI",
  iterations = 4, clusterParams = list(resolution = 1, sampleCells = 10000, n.start = 10),
  varFeatures = 25000, dimsToUse = 1:30
)
AML.mix <- addClusters(input = AML.mix, reducedDims = "IterativeLSI", method = "Seurat", name = "Clusters", resolution = 1)
AML.mix <- addUMAP(ArchRProj = AML.mix, reducedDims = "IterativeLSI", name = "UMAP", nNeighbors = 30, minDist = 1, metric = "cosine")
AML.mix <- addImputeWeights(AML.mix)
AML.mix$manual.cluster <- "AML"
AML.mix$manual.cluster[which(AML.mix$Clusters %in% c("C9", "C10"))] <- "B cell"
AML.mix$manual.cluster[which(AML.mix$Clusters %in% c("C16", "C17", "C18", "C19", "C20"))] <- "T cell"
AML.mix$manual.cluster[which(AML.mix$Clusters %in% c("C21", "C22", "C23"))] <- "NK"

AML.mix$manual.cluster[which(AML.mix$manual.cluster == "AML" & AML.mix$Sample == "AML1011_A")] <- "AML1011"
AML.mix$manual.cluster[which(AML.mix$manual.cluster == "AML" & AML.mix$Sample == "AML1012_A")] <- "AML1012"

# 1, 2 - Ery
# 4, 6, 7, 8 - Mono
# 11, 13, 14, 15 - HSC
# 3, 5 - GMP

saveArchRProject(AML.mix)
