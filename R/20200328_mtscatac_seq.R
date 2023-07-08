library(ggplot2)
library(ggrepel)
library(gplots)
library(Signac)
library(Seurat)

source('/Users/shaka87/dfci/scripts/20200122_colors_chorddiagram.R')
color_palette = colorRampPalette(c('white','red'))(100)

get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

get_coverage_reads_in_peaks_qc <- function(peak_matrix_h5, singlecell_csv, fragments_tsv, mgatk_experiment) {
  counts = Read10X_h5(filename = peak_matrix_h5)
  metadata = read.csv(singlecell_csv, header = TRUE, row.names = 1)
  atac = CreateSeuratObject(counts = counts, assay = 'peaks', project = 'Ibrutinib', 
                            min.cells = 1, meta.data = metadata)
  atac = SetFragments(object = atac, file = fragments_tsv)
  atac$pct_reads_in_peaks = atac$peak_region_fragments / atac$passed_filters * 100
  qc = data.frame(coverage = mgatk_experiment@colData@listData[["depth"]], 
                  reads_in_peaks = atac$pct_reads_in_peaks[which(colnames(atac) %in% colnames(mgatk_experiment))])
  qc$density = get_density(qc$coverage, qc$reads_in_peaks, n=100)
  return (qc)
}

plot_coverage_reads_in_peaks_qc <- function(qc, coverage_limits, reads_limits, pngfile) {
  p=ggplot(data=qc, aes(x=coverage,y=reads_in_peaks,color=density)) + geom_point(size=0.5) +
    geom_hline(yintercept = reads_limits, linetype='dashed') +
    geom_vline(xintercept = coverage_limits, linetype='dashed') +
    scale_x_log10('mtDNA Coverage',limits = c(1,200), breaks = c(1,10,20,50,100)) + 
    scale_y_continuous('% Reads in Peaks') +
    scale_color_gradientn(colours = c('blue','lightblue','yellow','red')) +
    theme_bw() +
    theme(panel.grid  = element_blank(),
          panel.border = element_blank(),
          legend.position = 'None',
          axis.line = element_line(color='black'),
          axis.ticks = element_line(color='black'),
          axis.text = element_text('Arial',size=14,color='black'),
          axis.title = element_text('Arial',size=18,color='black'))
  ggsave(filename = pngfile, width = 4, height = 4,dpi = 350, plot = p)
}

plot_concordance_vmr <- function(variants.qc, concordance_limits, vmr_limits, pngfile, n_cell_limit=5,
                                 color1='grey',color2='blue') {
  p=ggplot() +
    geom_point(data=variants.qc[which(variants.qc$n_cells > n_cell_limit & (log10(variants.qc$vmr) < vmr_limits
                                                                            | variants.qc$concordance < concordance_limits)),], 
               aes(x=concordance,y=log10(vmr)),size=0.5,color=color1) +
    geom_point(data=variants.qc[which(variants.qc$n_cells > n_cell_limit & log10(variants.qc$vmr) > vmr_limits
                                      & variants.qc$concordance > concordance_limits),], 
               aes(x=concordance,y=log10(vmr)),size=0.5,color=color2) +
    geom_hline(yintercept = vmr_limits, linetype='dashed') +
    geom_vline(xintercept = concordance_limits, linetype='dashed') +
    scale_x_continuous('Strand concordance',limits = c(-0.5,1)) + 
    scale_y_continuous('log10 VMR', limits = c(-5,0)) +
    theme_bw() +
    theme(panel.grid  = element_blank(),
          panel.border = element_blank(),
          legend.position = 'None',
          axis.line = element_line(color='black'),
          axis.ticks = element_line(color='black'),
          axis.text = element_text('Arial',size=14,color='black'),
          axis.title = element_text('Arial',size=18,color='black'))
  ggsave(pngfile, width = 4, height = 4, dpi = 350, plot = p)
}

relevant_mutations <- function(mutations_mgatk, variants.qc, vmr_limit, concordance_limit, n_cells_limit=5) {
  top.mutations = which(log10(variants.qc$vmr) > vmr_limit & variants.qc$concordance > concordance_limit &
                          variants.qc$n_cells > n_cells_limit)
  frequencies = as.matrix(mutations_mgatk@assays@data@listData[["allele_frequency"]][top.mutations,])
  mutations = mutations_mgatk@elementMetadata@listData[["variant"]][top.mutations]
  rownames(frequencies) = mutations
  colnames(frequencies) = mutations_mgatk@colData@rownames
  return (frequencies)
}

cluster_relevant_mutations <- function(frequencies, mutations_mgatk, dims = 1:10, k_param = 20, resolution = 0.8) {
  frequencies.seurat = t(frequencies)
  rownames(frequencies.seurat) = mutations_mgatk@colData@rownames
  frequencies.seurat = CreateSeuratObject(t(frequencies.seurat))
  frequencies.seurat = ScaleData(frequencies.seurat, do.scale = F)
  frequencies.seurat = FindVariableFeatures(frequencies.seurat)
  frequencies.seurat = RunPCA(frequencies.seurat)
  frequencies.seurat = FindNeighbors(frequencies.seurat, k.param = k_param, dims = dims, annoy.metric = 'cosine')
  frequencies.seurat = FindClusters(frequencies.seurat, resolution = resolution)
  
  cell.order = c()
  cells.cluster = c()
  for (cluster in names(summary(frequencies.seurat$seurat_clusters))) {
    cells.in.cluster = colnames(frequencies.seurat)[which(frequencies.seurat$seurat_clusters == cluster)]
    for (cell in cells.in.cluster) {
      cell.order = c(cell.order, which(colnames(frequencies) == cell))
    }
    cells.cluster = c(cells.cluster, rep(cluster,length(cells.in.cluster)))
  }
  frequencies = frequencies[,cell.order]
  cluster.colors = choose_colour_clone(as.integer(cells.cluster))
  return (list(frequencies, cells.cluster, cluster.colors))
}

cluster_relevant_mutations <- function(frequencies, dims = 1:10, k_param = 20, resolution = 0.8) {
  frequencies.seurat = t(frequencies)
  frequencies.seurat = CreateSeuratObject(t(frequencies.seurat))
  frequencies.seurat = ScaleData(frequencies.seurat, do.scale = F)
  frequencies.seurat = FindVariableFeatures(frequencies.seurat)
  frequencies.seurat = RunPCA(frequencies.seurat)
  frequencies.seurat = FindNeighbors(frequencies.seurat, k.param = k_param, dims = dims, annoy.metric = 'cosine')
  frequencies.seurat = FindClusters(frequencies.seurat, resolution = resolution)
  
  cell.order = c()
  cells.cluster = c()
  for (cluster in names(summary(frequencies.seurat$seurat_clusters))) {
    cells.in.cluster = colnames(frequencies.seurat)[which(frequencies.seurat$seurat_clusters == cluster)]
    for (cell in cells.in.cluster) {
      cell.order = c(cell.order, which(colnames(frequencies) == cell))
    }
    cells.cluster = c(cells.cluster, rep(cluster,length(cells.in.cluster)))
  }
  frequencies = frequencies[,cell.order]
  cluster.colors = choose_colour_clone(as.integer(cells.cluster))
  return (list(frequencies, cells.cluster, cluster.colors))
}
