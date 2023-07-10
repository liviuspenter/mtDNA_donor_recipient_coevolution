setwd('/Users/shaka87/dfci/asap_seq/')

library(BuenColors)
library(ComplexHeatmap)
library(data.table)
library(parallel)
library(stringr)
library(svglite)
library(ArchR)


AML.1010.mito = loadArchRProject('/Users/shaka87/ArchRProjects/AML.asap/AML.1010.mito/')
AML.1010.mito = addGroupCoverages(AML.1010.mito, threads = 12)
AML.1010.mito = addReproduciblePeakSet(AML.1010.mito, groupBy = 'Clusters', threads = 12)
AML.1010.mito = addPeakMatrix(AML.1010.mito, threads = 12)
AML.1010.mito = addMotifAnnotations(AML.1010.mito, motifSet = 'cisbp')
AML.1010.mito <- addBgdPeaks(AML.1010.mito)
AML.1010.mito = addDeviationsMatrix(AML.1010.mito, threads = 12)

AML.1012.mito = loadArchRProject('/Users/shaka87/ArchRProjects/AML.asap/AML.1012.mito/')
AML.1012.mito = addGroupCoverages(AML.1012.mito, threads = 12)
AML.1012.mito = addReproduciblePeakSet(AML.1012.mito, groupBy = 'Clusters', threads = 12)
AML.1012.mito = addPeakMatrix(AML.1012.mito, threads = 12)
AML.1012.mito = addMotifAnnotations(AML.1012.mito, motifSet = 'cisbp')
AML.1012.mito <- addBgdPeaks(AML.1012.mito)
AML.1012.mito = addDeviationsMatrix(AML.1012.mito, threads = 12)

AML.1026.mito = loadArchRProject('/Users/shaka87/ArchRProjects/AML.asap/AML.1026.mito/')
AML.1026.mito = addGroupCoverages(AML.1026.mito, threads = 12)
AML.1026.mito = addReproduciblePeakSet(AML.1026.mito, groupBy = 'Clusters', threads = 12)
AML.1026.mito = addPeakMatrix(AML.1026.mito, threads = 12)
AML.1026.mito = addMotifAnnotations(AML.1026.mito, motifSet = 'cisbp')
AML.1026.mito <- addBgdPeaks(AML.1026.mito)
AML.1026.mito = addDeviationsMatrix(AML.1026.mito, threads = 12)

# Differential chromatin analysis
MotifMatrix_tmp.1010 = getMatrixFromProject(AML.1010.mito, useMatrix = "MotifMatrix", useSeqnames = NULL, binarize = FALSE, threads = 12)
MotifMatrix_tmp.1012 = getMatrixFromProject(AML.1012.mito, useMatrix = "MotifMatrix", useSeqnames = NULL, binarize = FALSE, threads = 12)
MotifMatrix_tmp.1026 = getMatrixFromProject(AML.1026.mito, useMatrix = "MotifMatrix", useSeqnames = NULL, binarize = FALSE, threads = 12)

zscore_matrix.1010 = as.matrix(assays(MotifMatrix_tmp.1010)$z) %>% as.data.frame()
zscore_matrix.1012 = as.matrix(assays(MotifMatrix_tmp.1012)$z) %>% as.data.frame()
zscore_matrix.1026 = as.matrix(assays(MotifMatrix_tmp.1026)$z) %>% as.data.frame()

VarDeviation.1010 = getVarDeviations(AML.1010.mito[which(AML.1010.mito$manual.clusters %in% c('HSCT', 'GMP', 'Mono', 'Ery'))], 
                                     name = "MotifMatrix", plot = FALSE) %>% as.data.frame()
VarDeviation.1012 = getVarDeviations(AML.1012.mito[which(AML.1012.mito$manual.clusters %in% c('HSCT', 'GMP', 'Mono', 'Ery'))], 
                                     name = "MotifMatrix", plot = FALSE) %>% as.data.frame()
VarDeviation.1026 = getVarDeviations(AML.1026.mito[which(AML.1026.mito$manual.clusters %in% c('HSCT', 'GMP', 'Mono', 'Ery'))], 
                                     name = "MotifMatrix", plot = FALSE) %>% as.data.frame()

Top_TF_target = intersect(head(VarDeviation.1010$name, n=50), head(VarDeviation.1012$name, n=50))
Top_TF_target = intersect(Top_TF_target, head(VarDeviation.1026$name, n=50))

HSC.cells.1010 = AML.1010.mito$cellNames[which(AML.1010.mito$manual.clusters == 'HSCT')]
HSC.cells.1010 = HSC.cells.1010[sample(length(HSC.cells.1010), size = 200)]
GMP.cells.1010 = AML.1010.mito$cellNames[which(AML.1010.mito$manual.clusters == 'GMP')]
GMP.cells.1010 = GMP.cells.1010[sample(length(GMP.cells.1010), size = 200)]
Mono.cells.1010 = AML.1010.mito$cellNames[which(AML.1010.mito$manual.clusters == 'Mono')]
Mono.cells.1010 = Mono.cells.1010[sample(length(Mono.cells.1010), size = 200)]
Ery.cells.1010 = AML.1010.mito$cellNames[which(AML.1010.mito$manual.clusters == 'erythroid')]
Ery.cells.1010 = Ery.cells.1010[sample(length(Ery.cells.1010), size = 200)]

zscore_part.1010 = zscore_matrix.1010[Top_TF_target,c(HSC.cells.1010, GMP.cells.1010, Mono.cells.1010, Ery.cells.1010)]
zscore_part.1010[zscore_part.1010>2]  = 2
zscore_part.1010[-2>zscore_part.1010] = -2

HSC.cells.1012 = AML.1012.mito$cellNames[which(AML.1012.mito$manual.clusters == 'HSCT')]
HSC.cells.1012 = HSC.cells.1012[sample(length(HSC.cells.1012), size = 200)]
GMP.cells.1012 = AML.1012.mito$cellNames[which(AML.1012.mito$manual.clusters == 'GMP')]
GMP.cells.1012 = GMP.cells.1012[sample(length(GMP.cells.1012), size = 200)]
Mono.cells.1012 = AML.1012.mito$cellNames[which(AML.1012.mito$manual.clusters == 'Mono')]
Mono.cells.1012 = Mono.cells.1012[sample(length(Mono.cells.1012), size = 200)]
Ery.cells.1012 = AML.1012.mito$cellNames[which(AML.1012.mito$manual.clusters == 'erythroid')]
Ery.cells.1012 = Ery.cells.1012[sample(length(Ery.cells.1012), size = 200)]

zscore_part.1012 = zscore_matrix.1012[Top_TF_target,c(HSC.cells.1012, GMP.cells.1012, Mono.cells.1012, Ery.cells.1012)]
zscore_part.1012[zscore_part.1012>2]  = 2
zscore_part.1012[-2>zscore_part.1012] = -2

HSC.cells.1026 = AML.1026.mito$cellNames[which(AML.1026.mito$manual.clusters == 'HSCT')]
HSC.cells.1026 = HSC.cells.1026[sample(length(HSC.cells.1026), size = 200)]
GMP.cells.1026 = AML.1026.mito$cellNames[which(AML.1026.mito$manual.clusters == 'GMP')]
GMP.cells.1026 = GMP.cells.1026[sample(length(GMP.cells.1026), size = 200)]
Mono.cells.1026 = AML.1026.mito$cellNames[which(AML.1026.mito$manual.clusters == 'Mono')]
Mono.cells.1026 = Mono.cells.1026[sample(length(Mono.cells.1026), size = 200)]
Ery.cells.1026 = AML.1026.mito$cellNames[which(AML.1026.mito$manual.clusters == 'erythroid')]
Ery.cells.1026 = Ery.cells.1026[sample(length(Ery.cells.1026), size = 200)]

zscore_part.1026 = zscore_matrix.1026[Top_TF_target,c(HSC.cells.1026, GMP.cells.1026, Mono.cells.1026, Ery.cells.1026)]
zscore_part.1026[zscore_part.1026>2]  = 2
zscore_part.1026[-2>zscore_part.1026] = -2

zscore_part = cbind(zscore_part.1010, zscore_part.1026)
zscore_part = cbind(zscore_part,zscore_part.1012)

TF.order = c('CTCF_177', 'RUNX1_733', 'RUNX2_732', 'RUNX3_731', 'CBFB_801', 'CEBPB_140', 'CEBPA_155', 'CEBPD_152', 'PITX2_504', 
             'ATF4_122', 'HLF_112', 'CEBPE_107', 'CEBPG_128', 'SPIC_344', 'SPIB_336', 'SPI1_322', 'BCL11A_194', 'BCL11B_825',
             'GATA1_383', 'GATA2_388', 'GATA3_384', 'GATA4_386', 'GATA5_385', 'GATA6_387', 'MECOM_169')

ha = HeatmapAnnotation(Cluster = rep(c(rep('HSC', 200), rep('GMP', 200), rep('Mono', 200), rep('erythroid', 200)), 3), 
                       col = list(Cluster = c('HSC' = RColorBrewer::brewer.pal(n=5, name='Reds')[5], 
                               'GMP' = RColorBrewer::brewer.pal(n=5, name='Reds')[4], 
                               'Mono' = RColorBrewer::brewer.pal(n=5, name='Reds')[3],
                               'erythroid' = RColorBrewer::brewer.pal(n=5, name='Reds')[2])), 
                       simple_anno_size = unit(5, 'pt'), border = T)
blueYellow = c("1"="#352A86","2"="#343DAE","3"="#0262E0","4"="#1389D2","5"="#2DB7A3","6"="#A5BE6A","7"="#F8BA43","8"="#F6DA23","9"="#F8FA0D")
col_fun = circlize::colorRamp2(breaks = seq(-2,2,4/8), colors = blueYellow)
svglite::svglite('./figures/heatmaps/20220302_AML_TF.svg', width = 5, height = 3)
Heatmap(as.matrix(zscore_part[TF.order,]), col = col_fun, 
        cluster_rows = F, cluster_columns = F, show_row_names = T, show_column_names = F, show_row_dend = F, 
        top_annotation = ha, column_split = c(rep('AML1010', 800), rep('AML1026', 800), rep('AML1012', 800)), 
        use_raster = T, raster_quality = 10, border = T, row_names_gp = gpar(fontsize=8), row_names_side = 'left', 
        row_labels = stringr::str_split_fixed(TF.order, pattern = '_', n=2)[,1])
dev.off()

TSB.so = readRDS('./data/objects/20210601_TSB.so.rds')

TSB.data = Seurat::GetAssayData(TSB.so, slot = 'scale.data')
colnames(TSB.data) = paste0(colnames(TSB.data), '-1')
col_fun = circlize::colorRamp2(breaks = seq(-2,2,4/8), colors = BuenColors::jdb_palette(name = 'brewer_yes', n=9))
svglite::svglite('./figures/heatmaps/20220302_AML_ASAP.svg', width = 5, height = 1.5)
Heatmap(as.matrix(TSB.data[c('CD117', 'CD38', 'CD33', 'CD14', 'CD11c', 'CD16', 'CD39'),colnames(zscore_part)]), col = col_fun, 
        cluster_rows = F, cluster_columns = F, show_row_names = T, show_column_names = F, show_row_dend = F, 
        top_annotation = ha, column_split = c(rep('AML1010', 800), rep('AML1012', 800), rep('AML1026', 800)), 
        use_raster = T, raster_quality = 10, 
        border = T, row_names_gp = gpar(fontsize=8), row_names_side = 'left')
dev.off()



trajectory.erythroid = c('C5', 'C4','C7')
trajectory.mono = c('C5', 'C6','C2', 'C3', 'C1')

AML.1012.mito = addTrajectory(AML.1012.mito, name = 'trajectory.erythroid', groupBy = 'Clusters', trajectory = trajectory.erythroid, embedding = 'UMAP', force = T)
AML.1012.mito = addTrajectory(AML.1012.mito, name = 'trajectory.mono', groupBy = 'Clusters', trajectory = trajectory.mono, embedding = 'UMAP', force = T)

p=plotTrajectory(AML.1012.mito, trajectory = 'trajectory.erythroid', colorBy = 'cellColData', name = 'trajectory.erythroid', rastr = T, plotAs = 'points')
plotPDF(p[[1]], name = 'AML1012.trajectory.erythroid', ArchRProj = AML.1012.mito, addDOC = F)
p=plotTrajectory(AML.1012.mito, trajectory = 'trajectory.mono', colorBy = 'cellColData', name = 'trajectory.mono')
plotPDF(p[[1]], name = 'AML1012.trajectory.mono', ArchRProj = AML.1012.mito, addDOC = F)

trajGSM.ery <- getTrajectory(ArchRProj = AML.1012.mito, name = "trajectory.erythroid", useMatrix = "GeneScoreMatrix", log2Norm = TRUE)
trajMM.ery <- getTrajectory(ArchRProj = AML.1012.mito, name = "trajectory.erythroid", useMatrix = "MotifMatrix", log2Norm = FALSE)
cor_traj.ery = correlateTrajectories(trajGSM.ery, trajMM.ery)

trajGSM.ery.2 = trajGSM.ery[cor_traj.ery[[1]]$name1, ]
trajMM.ery.2 = trajMM.ery[cor_traj.ery[[1]]$name2, ]

match(rownames(trajGSM.ery.2), rownames(trajMM.ery.2))
ht1 <- plotTrajectoryHeatmap(trajGSM.ery.2,  pal = paletteContinuous(set = "horizonExtra"),  varCutOff = 0)
ht2 <- plotTrajectoryHeatmap(trajMM.ery.2,  pal = paletteContinuous(set = "horizonExtra"),  varCutOff = 0)


trajCombined <- trajGSM.ery.2
trajCombined <- t(apply(assay(trajGSM.ery.2), 1, scale)) + t(apply(assay(trajMM.ery.2), 1, scale))


p2 <- plotTrajectoryHeatmap(trajGSM.ery,  pal = paletteContinuous(set = "horizonExtra"))
p3 <- plotTrajectoryHeatmap(trajMM.ery,  pal = paletteContinuous(set = "horizonExtra"))