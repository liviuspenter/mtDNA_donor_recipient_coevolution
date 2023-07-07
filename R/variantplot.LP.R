
# tweak VariantPlot function
VariantPlot.LP <- function(variants, min.cells = 2, concordance.threshold = 0.65, vmr.threshold = 0.01, pt.size=0.5) {
  high.conf <- variants[variants$n_cells_conf_detected >= min.cells, ]
  high.conf$pos = 'non.significant'
  high.conf$pos[high.conf$vmr > vmr.threshold & high.conf$strand_correlation > concordance.threshold] = 'significant'
  high.conf$pos[high.conf$vmr < vmr.threshold & high.conf$strand_correlation > concordance.threshold] = 'germline'
  
  p <- ggplot(data = high.conf, aes(x=strand_correlation, y=log10(vmr), color=pos)) +
    geom_point(size=pt.size) +
    labs(x = "Strand concordance", y = "Variance-mean ratio") +
    geom_vline(xintercept = concordance.threshold, color = "black", linetype = 2) +
    geom_hline(yintercept = log10(vmr.threshold), color = "black", linetype = 2) +
    scale_color_manual(values = c('non.significant' = "grey", 'significant' = "firebrick", 'germline' = "blue")) +
    scale_y_continuous('log10(VMR)') +
    theme_classic() +
    theme(legend.position = "none",
          axis.text = element_text('Arial', size=10, color='black'),
          axis.title = element_text('Arial', size=10, color='black'))
  return(p)
}