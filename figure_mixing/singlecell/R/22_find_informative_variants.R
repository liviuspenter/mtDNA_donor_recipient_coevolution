library(ArchR)
library(ggplot2)
library(parallel)
library(Seurat)
library(Signac)

source("./R/variantplot.LP.R")

AML.mix <- loadArchRProject("./data/mixing/AML.mix/")

# extract mtDNA variants from "germline" controls
mito.data.1 <- ReadMGATK("./data/10026/Pool148_13.mgatk/final/")
mito.data.2 <- ReadMGATK("./data/10026/Pool148_17.mgatk/final/")
variant.df.1 <- IdentifyVariants(mito.data.1$counts, refallele = mito.data.1$refallele)
write.table(variant.df.1, quote = F, sep = "\t", row.names = F, file = "./data/mixing/singlecell/artificial_mixing_AML/Pool148_13.variants.csv")
variant.df.2 <- IdentifyVariants(mito.data.2$counts, refallele = mito.data.2$refallele)
write.table(variant.df.2, quote = F, sep = "\t", row.names = F, file = "./data/mixing/singlecell/artificial_mixing_AML/Pool148_17.variants.csv")

VariantPlot.LP(variant.df.1)
ggsave("./figure_mixing/singlecell/figures/20240603_germline1.svg", width = 1.5, height = 1.5)
VariantPlot.LP(variant.df.2)
ggsave("./figure_mixing/singlecell/figures/20240603_germline2.svg", width = 1.5, height = 1.5)

variants.1 <- variant.df.1$variant[which(variant.df.1$vmr < 0.01 & variant.df.1$strand_correlation > 0.65)]
variants.2 <- variant.df.2$variant[which(variant.df.2$vmr < 0.01 & variant.df.2$strand_correlation > 0.65)]
informative.variants.1 <- setdiff(variants.1, variants.2)
informative.variants.2 <- setdiff(variants.2, variants.1)

informative.variants <- data.frame(variant = c(informative.variants.1, informative.variants.2))
informative.variants$individual <- ifelse(informative.variants$variant %in% variants.1, "AML1011", "AML1012")

write.table(informative.variants, quote = F, sep = "\t", row.names = F, file = "./data/mixing/singlecell/artificial_mixing_AML/informative_variants.csv")

# test for false-positive calls in opposite germline and unannotated cells
VAFs.1 <- AlleleFreq(mito.data.1$counts, variants = informative.variants$variant, assay = NULL)
colnames(VAFs.1) <- paste0("AML1011_A#", colnames(VAFs.1))
VAFs.2 <- AlleleFreq(mito.data.2$counts, variants = informative.variants$variant, assay = NULL)
colnames(VAFs.2) <- paste0("AML1012_A#", colnames(VAFs.2))

df <- data.frame(
  bc = c(colnames(VAFs.1), colnames(VAFs.2)),
  variants.1 = c(
    colMeans(VAFs.1[informative.variants.1, ]),
    colMeans(VAFs.2[informative.variants.1, ])
  ),
  variants.2 = c(
    colMeans(VAFs.1[informative.variants.2, ]),
    colMeans(VAFs.2[informative.variants.2, ])
  ),
  sample = c(rep("AML1011_A", ncol(VAFs.1)), rep("AML1012_A", ncol(VAFs.2)))
)

df <- df[intersect(df$bc, AML.mix$cellNames), ]

df$coverage <- c(
  mito.data.1$depth[stringr::str_split_fixed(df$bc[which(df$sample == "AML1011_A")], pattern = "#", n = 2)[, 2], ],
  mito.data.2$depth[stringr::str_split_fixed(df$bc[which(df$sample == "AML1012_A")], pattern = "#", n = 2)[, 2], ]
)
df$individual <- "none"
df[which(df$variants.1 > 0.8 & df$variants.2 < 0.2), "individual"] <- "AML1011_A"
df[which(df$variants.2 > 0.8 & df$variants.1 < 0.2), "individual"] <- "AML1012_A"

# get information for UMAPs from ArchR object
df <- cbind(df, getEmbedding(AML.mix)[df$bc, ])
df$manual.cluster <- AML.mix[df$bc]$manual.cluster

# plot stuff

ggplot(df, aes(x = 100 * variants.1, y = 100 * variants.2)) +
  geom_hline(yintercept = c(20, 80)) +
  geom_vline(xintercept = c(20, 80)) +
  ggrastr::rasterize(geom_point(aes(color = individual), size = 0.5), dpi = 600) +
  scale_x_continuous("% informative variants 1") +
  scale_y_continuous("% informative variants 2") +
  scale_color_manual(values = c("none" = "black", "AML1011_A" = "orange", "AML1012_A" = "purple")) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black")
  )
ggsave("./figure_mixing/singlecell/figures/20240603_VAF_plot.svg", width = 2, height = 2)

p <- ggplot(df, aes(x = individual, y = coverage)) +
  ggrastr::rasterize(geom_jitter(size = 0.5, color = "grey"), dpi = 600) +
  geom_boxplot(outlier.color = NA, fill = NA, aes(color = individual)) +
  scale_x_discrete(labels = c("AML1011", "AML1012", "unannotated")) +
  scale_y_log10("mtDNA coverage") +
  scale_color_manual(values = c("none" = "black", "AML1011_A" = "orange", "AML1012_A" = "purple")) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black"),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )
ggsave("./figure_mixing/singlecell/figures/20240603_mtDNA_coverage.svg", width = 1.5, height = 2, plot = p)

ggplot() +
  geom_point(
    data = df, size = 0.5,
    aes(x = `IterativeLSI#UMAP_Dimension_1`, y = `IterativeLSI#UMAP_Dimension_2`, color = manual.cluster)
  ) +
  scale_color_manual(values = c("NK" = "purple", "T cell" = "blue", "B cell" = "yellow", "AML1011" = "orange", "AML1012" = "purple")) +
  theme_classic() +
  theme(legend.position = "none") +
  NoAxes()
ggsave("./figure_mixing/singlecell/figures/20240603_UMAP.png", width = 4, height = 4, dpi = 600)

ggplot() +
  geom_point(
    data = df, size = 0.5,
    aes(x = `IterativeLSI#UMAP_Dimension_1`, y = `IterativeLSI#UMAP_Dimension_2`, color = manual.cluster), alpha = 0.05
  ) +
  geom_point(
    data = df[which(df$individual == "none"), ], size = 1,
    aes(x = `IterativeLSI#UMAP_Dimension_1`, y = `IterativeLSI#UMAP_Dimension_2`), color = "black"
  ) +
  scale_color_manual(values = c("NK" = "purple", "T cell" = "blue", "B cell" = "yellow", "AML1011" = "orange", "AML1012" = "purple")) +
  theme_classic() +
  theme(legend.position = "none") +
  NoAxes()
ggsave("./figure_mixing/singlecell/figures/20240603_UMAP_not_annotated.png", width = 4, height = 4, dpi = 600)
