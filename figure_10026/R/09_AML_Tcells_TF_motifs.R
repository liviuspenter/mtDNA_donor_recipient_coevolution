# transcription factor motifs in AML T cells

library(ArchR)
library(gplots)
library(Seurat)

source('./figure_10026/R/Tcell.colors.R')

AML.Tcell = loadArchRProject('./data/10026/AML.Tcell/')

MotifMatrix_tmp = getMatrixFromProject(ArchRProj = AML.Tcell,
                                       useMatrix = "MotifMatrix",
                                       useSeqnames = NULL,
                                       verbose = TRUE,
                                       binarize = FALSE,
                                       threads = 12)

zscore_matrix     = as.matrix(assays(MotifMatrix_tmp)$z)          %>% as.data.frame()
deviations_matrix = as.matrix(assays(MotifMatrix_tmp)$deviations) %>% as.data.frame()

VarDeviation = getVarDeviations(AML.Tcell, name = "MotifMatrix", plot = FALSE) %>% as.data.frame()
Top_TF_target = head(VarDeviation$name, n=75)

zscore_part = zscore_matrix[Top_TF_target,]
zscore_part[zscore_part>2]  = 2
zscore_part[-2>zscore_part] = -2

rownames(zscore_part) = stringr::str_split_fixed(rownames(zscore_part), pattern = '_', n=2)[,1]

blueYellow = c("1"="#352A86","2"="#343DAE","3"="#0262E0","4"="#1389D2","5"="#2DB7A3","6"="#A5BE6A","7"="#F8BA43","8"="#F6DA23","9"="#F8FA0D")
col_list   = colorRampPalette(blueYellow)

cells.order = c()
side.colors = c()
for (cluster in unique(names(Tcell.colors))) {
  cells.in.cluster = AML.Tcell$cellNames[which(AML.Tcell$manual_clusters == cluster)]
  if (length(cells.in.cluster) > 200) {
    cells.in.cluster = cells.in.cluster[sample(length(cells.in.cluster), 200)]
  }
  cells.order = c(cells.order, cells.in.cluster)
  side.colors = c(side.colors, rep(Tcell.colors[cluster], length(cells.in.cluster)))
}

svglite::svglite('./figure_10026/figures/heatmaps/20210603_Tcell_TF.svg', width = 5, height = 5)
heatmap.2(as.matrix(zscore_part[,cells.order]), trace = 'none', Colv = F, dendrogram = 'none', ColSideColors = side.colors, col = col_list, 
          labCol = F, key = F, lwid = c(0.1, 5), lhei = c(0.1,5), margins = c(1,5), cexRow = 0.5)
dev.off()