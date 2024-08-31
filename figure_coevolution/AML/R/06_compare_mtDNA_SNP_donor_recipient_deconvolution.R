library(ComplexHeatmap)
library(dplyr)
library(ggplot2)
library(ggExtra)
library(Seurat)

"%ni%" <- function(x, y) !("%in%"(x, y))

source("./R/filehandling.functions.tapestri.R")

cluster.colors <- c(
  "HSC" = RColorBrewer::brewer.pal(name = "YlOrRd", n = 5)[5],
  "Progenitor" = RColorBrewer::brewer.pal(name = "YlOrRd", n = 5)[4],
  "Mono" = RColorBrewer::brewer.pal(name = "YlOrRd", n = 5)[3],
  "Erythroid" = RColorBrewer::brewer.pal(name = "YlOrRd", n = 5)[2],
  "CD4" = "lightblue",
  "CD8" = "darkblue",
  "Plasma" = "darkgreen"
)

# read data using helper functions
AML1026.1.vafs <- read_vafs("./data/coevolution/1026_1/AF.csv", "AML1026.1")
AML1026.3.vafs <- read_vafs("./data/coevolution/1026_3/AF.csv", "AML1026.3")

AML1026.1.depths <- read_vafs("./data/coevolution/1026_1/DP.csv", "AML1026.1")
AML1026.3.depths <- read_vafs("./data/coevolution/1026_3/DP.csv", "AML1026.3")

AML1026.1.vafs.combined <- merge_with_mgatk(AML1026.1.vafs, "./data/coevolution/vafs/20220211_AML1026.1_mtDNA_vafs.csv", "AML1026.1", truncate_barcodes = T)
AML1026.3.vafs.combined <- merge_with_mgatk(AML1026.3.vafs, "./data/coevolution/vafs/20220211_AML1026.3_mtDNA_vafs.csv", "AML1026.3", truncate_barcodes = T)

AML1026.vafs <- rbind(
  AML1026.1.vafs[, intersect(colnames(AML1026.1.vafs), colnames(AML1026.3.vafs))],
  AML1026.3.vafs[, intersect(colnames(AML1026.1.vafs), colnames(AML1026.3.vafs))]
)
AML1026.vafs.combined <- rbind(
  AML1026.1.vafs.combined[, intersect(colnames(AML1026.1.vafs.combined), colnames(AML1026.3.vafs.combined))],
  AML1026.3.vafs.combined[, intersect(colnames(AML1026.1.vafs.combined), colnames(AML1026.3.vafs.combined))]
)

AML1026.depths.combined <- rbind(
  AML1026.1.depths[, intersect(colnames(AML1026.1.depths), colnames(AML1026.3.depths))],
  AML1026.3.depths[, intersect(colnames(AML1026.1.depths), colnames(AML1026.3.depths))]
)

## combined deconvolution variants
# recipient
variants.1 <- c("NF1:chr17:29559932:C/A", "TP53:chr17:7579801:G/C", "TP53:chr17:7578115:T/C", "16304T>C", "4336T>C")
# donor
variants.2 <- c(
  "SF3B1:chr2:198267770:G/GAA",
  "NPM1:chr5:170837457:A/G", "KDM6A:chrX:44938563:G/A", "DNMT3A:chr2:25463483:G/A",
  "BRAF:chr7:140449071:C/G", "FLT3:chr13:28592546:T/C", "NF1:chr17:29483195:G/C", "16294C>T", "16296C>T"
)

## separate deconvolution variants
# recipient
variants.mtDNA.1 <- c("16304T>C", "4336T>C")
# donor
variants.mtDNA.2 <- c("16294C>T", "16296C>T")
# recipient
variants.somatic.1 <- c("NF1:chr17:29559932:C/A", "TP53:chr17:7579801:G/C", "TP53:chr17:7578115:T/C")
# donor
variants.somatic.2 <- c(
  "SF3B1:chr2:198267770:G/GAA",
  "NPM1:chr5:170837457:A/G", "KDM6A:chrX:44938563:G/A", "DNMT3A:chr2:25463483:G/A",
  "BRAF:chr7:140449071:C/G", "FLT3:chr13:28592546:T/C", "NF1:chr17:29483195:G/C"
)

DP1 <- as.data.frame(data.table::fread("./data/coevolution/1026_1/DP.csv"))
DP1$Barcode <- paste0("AML1026.1#", DP1$Barcode)
DP1[DP1 == 0] <- NA
DP3 <- as.data.frame(data.table::fread("./data/coevolution/1026_3/DP.csv"))
DP3$Barcode <- paste0("AML1026.3#", DP3$Barcode)
DP3[DP3 == 0] <- NA
DP <- rbind(DP1[c("Sample", "Barcode", colnames(AML1026.vafs))], DP3[c("Sample", "Barcode", colnames(AML1026.vafs))])

VARIANT_CUTOFF <- 10
variant.df <- data.frame(
  variants.1 = rowMeans(AML1026.vafs.combined[, variants.1]),
  variants.2 = rowMeans(AML1026.vafs.combined[, variants.2]),
  variants.mtDNA.1 = rowMeans(AML1026.vafs.combined[, variants.mtDNA.1]),
  variants.mtDNA.2 = rowMeans(AML1026.vafs.combined[, variants.mtDNA.2]),
  variants.somatic.1 = rowMeans(AML1026.vafs.combined[, variants.somatic.1]),
  variants.somatic.2 = rowMeans(AML1026.vafs.combined[, variants.somatic.2]),
  barcode = rownames(AML1026.vafs)
)

variant.df$annotation <- "none"
variant.df$annotation[which(variant.df$variants.1 > VARIANT_CUTOFF & variant.df$variants.2 > VARIANT_CUTOFF)] <- "doublet"
variant.df$annotation[which(variant.df$variants.1 > VARIANT_CUTOFF & variant.df$variants.2 < VARIANT_CUTOFF)] <- "recipient"
variant.df$annotation[which(variant.df$variants.1 < VARIANT_CUTOFF & variant.df$variants.2 > VARIANT_CUTOFF)] <- "donor"

variant.df$annotation.somatic <- "none"
variant.df$annotation.somatic[which(variant.df$variants.somatic.1 > VARIANT_CUTOFF & variant.df$variants.somatic.2 > VARIANT_CUTOFF)] <- "doublet"
variant.df$annotation.somatic[which(variant.df$variants.somatic.1 > VARIANT_CUTOFF & variant.df$variants.somatic.2 < VARIANT_CUTOFF)] <- "recipient"
variant.df$annotation.somatic[which(variant.df$variants.somatic.1 < VARIANT_CUTOFF & variant.df$variants.somatic.2 > VARIANT_CUTOFF)] <- "donor"

variant.df$annotation.mtDNA <- "none"
variant.df$annotation.mtDNA[which(variant.df$variants.mtDNA.1 > 20 & variant.df$variants.mtDNA.2 > 20)] <- "doublet"
variant.df$annotation.mtDNA[which(variant.df$variants.mtDNA.1 > 80 & variant.df$variants.mtDNA.2 < 20)] <- "recipient"
variant.df$annotation.mtDNA[which(variant.df$variants.mtDNA.1 < 20 & variant.df$variants.mtDNA.2 > 80)] <- "donor"

# visualize mtDNA-based donor-recipient deconvolution
ggplot(variant.df, aes(x = variants.mtDNA.1, y = variants.mtDNA.2)) +
  ggrastr::rasterize(geom_point(aes(color = annotation.mtDNA), size = 0.5), dpi = 600) +
  scale_x_continuous("% mtDNA recipient") +
  scale_y_continuous("% mtDNA donor") +
  scale_color_manual(values = c("none" = "grey", "doublet" = "black", recipient = "orange", donor = "purple")) +
  geom_hline(yintercept = c(20, 80)) +
  geom_vline(xintercept = c(20, 80)) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black")
  )
ggsave("./figure_coevolution/AML/figures/20230428_donor_recipient_mtDNA.svg", width = 1.5, height = 1.5)

# visualize SNP-based donor-recipient deconvolution
ggplot(variant.df, aes(x = variants.somatic.1, y = variants.somatic.2)) +
  ggrastr::rasterize(geom_point(aes(color = annotation.somatic), size = 0.5), dpi = 600) +
  scale_x_continuous("% VAF recipient") +
  scale_y_continuous("% VAF donor") +
  scale_color_manual(values = c("none" = "grey", "doublet" = "black", recipient = "orange", donor = "purple")) +
  geom_hline(yintercept = c(10)) +
  geom_vline(xintercept = c(10)) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black")
  )
ggsave("./figure_coevolution/AML/figures/20230428_donor_recipient_somatic_mutation.svg", width = 1.5, height = 1.5)

recipient.cells <- rownames(variant.df)[which(variant.df$annotation.mtDNA == "recipient")]
recipient.cells <- recipient.cells[sample(length(recipient.cells), 100)]
donor.cells <- rownames(variant.df)[which(variant.df$annotation.mtDNA == "donor")]
donor.cells <- donor.cells[sample(length(donor.cells), 100)]
none.cells <- rownames(variant.df)[which(variant.df$annotation.mtDNA == "none")]
doublet.cells <- rownames(variant.df)[which(variant.df$annotation.mtDNA == "doublet")]

# visualize heatmap with mtDNA- and SNP-based donor-recipient deconvolution
col_fun <- circlize::colorRamp2(breaks = seq(0, 100, 100 / 8), colors = BuenColors::jdb_palette(name = "solar_rojos"))
ha <- columnAnnotation(
  annotation = variant.df[
    c(recipient.cells, donor.cells, none.cells, doublet.cells),
    "annotation.mtDNA"
  ],
  col = list("annotation" = c(
    "none" = "grey", "doublet" = "black",
    "recipient" = "orange", "donor" = "purple"
  )),
  simple_anno_size = unit(5, "pt"),
  border = T
)
svglite::svglite("./figure_coevolution/AML/figures/20230428_donor_recipient_heatmap.svg", width = 4, height = 2)
Heatmap(
  t(AML1026.vafs.combined[
    c(recipient.cells, donor.cells, none.cells, doublet.cells),
    c(variants.mtDNA.1, variants.mtDNA.2, variants.somatic.1, variants.somatic.2)
  ]),
  show_row_dend = F, show_column_dend = F, cluster_rows = F, cluster_columns = F,
  show_row_names = F, show_column_names = F,
  column_split = factor(variant.df[c(recipient.cells, donor.cells, none.cells, doublet.cells), "annotation.mtDNA"],
    levels = c("recipient", "donor", "doublet", "none")
  ),
  col = col_fun, border = T, row_split = factor(c(rep("mtDNA", 4), rep("somatic", 10)),
    levels = c("mtDNA", "somatic")
  ),
  top_annotation = ha, use_raster = T, raster_quality = 10
)
dev.off()

df <- data.frame()
for (mtDNA.annotation in unique(variant.df$annotation.mtDNA)) {
  for (somatic.annotation in unique(variant.df$annotation.somatic)) {
    df <- rbind(df, data.frame(
      mtDNA.annotation = mtDNA.annotation,
      somatic.annotation = somatic.annotation,
      cells = length(which(variant.df$annotation.mtDNA == mtDNA.annotation &
        variant.df$annotation.somatic == somatic.annotation))
    ))
  }
}

###
# track NRAS and SF3B1 mutation across recipient-derived cell types
###

so.13 <- readRDS("./data/coevolution/objects/20220112_AML1026.rds")
so.13 <- RenameCells(so.13, new.names = paste0("AML", colnames(so.13)))

df <- data.frame(
  bc = rownames(AML1026.vafs.combined),
  NRAS.depth = AML1026.depths.combined$`NRAS:chr1:115258745:C/G`,
  NRAS = AML1026.vafs.combined$`NRAS:chr1:115258745:C/G`,
  SF3B1.depth = AML1026.depths.combined$`SF3B1:chr2:198266512:C/T`,
  SF3B1 = AML1026.vafs.combined$`SF3B1:chr2:198266512:C/T`
)
rownames(df) <- df$bc
df$annotation.mtDNA <- variant.df[rownames(df), "annotation.mtDNA"]
df$manual.cluster <- so.13$manual.cluster[df$bc]
df$sample <- so.13$orig.ident[df$bc]

stats <- df %>%
  filter(df$annotation.mtDNA %in% c("recipient")) %>%
  group_by(sample, manual.cluster) %>%
  summarize(
    NRAS.mut = length(which(NRAS != 0)),
    NRAS.wt = length(which(NRAS == 0)),
    SF3B1.mut = length(which(SF3B1 != 0)),
    SF3B1.wt = length(which(SF3B1 == 0))
  )
stats$NRAS.freq <- stats$NRAS.mut / (stats$NRAS.mut + stats$NRAS.wt)
stats$SF3B1.freq <- stats$SF3B1.mut / (stats$SF3B1.mut + stats$SF3B1.wt)

stats <- stats %>% tidyr::pivot_longer(names_to = "gene", values_to = "freq", cols = c("SF3B1.freq", "NRAS.freq"))
stats$sample.gene <- factor(paste0(stats$sample, ".", stats$gene),
  levels = c(
    "1026.1.NRAS.freq", "1026.3.NRAS.freq",
    "1026.1.SF3B1.freq", "1026.3.SF3B1.freq"
  )
)

# visualize regardless of sequencing depth
ggplot(
  stats[which(stats$manual.cluster %ni% c("CD4", "CD8", "Plasma")), ],
  aes(x = sample.gene, y = 100 * freq, color = manual.cluster)
) +
  geom_point() +
  geom_line(aes(group = manual.cluster)) +
  scale_x_discrete(labels = c("Screening", "Treatment", "Screening", "Treatment")) +
  scale_y_continuous("% recipient positive", limits = c(0, 100)) +
  scale_color_manual(values = cluster.colors) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black"),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )
ggsave("./figure_coevolution/AML/figures/20230428_SF3B1_NRAS.svg", width = 1.3, height = 1.5)

# filter cells with sequencing depth >10x
stats.filtered <- df %>%
  filter(df$annotation.mtDNA %in% c("recipient")) %>%
  group_by(sample, manual.cluster) %>%
  summarize(
    NRAS.mut = length(which(NRAS != 0 & NRAS.depth >= 10)),
    NRAS.wt = length(which(NRAS == 0 & NRAS.depth >= 10)),
    SF3B1.mut = length(which(SF3B1 != 0 & SF3B1.depth >= 10)),
    SF3B1.wt = length(which(SF3B1 == 0 & SF3B1.depth >= 10))
  )
stats.filtered$NRAS.freq <- stats.filtered$NRAS.mut / (stats.filtered$NRAS.mut + stats.filtered$NRAS.wt)
stats.filtered$SF3B1.freq <- stats.filtered$SF3B1.mut / (stats.filtered$SF3B1.mut + stats.filtered$SF3B1.wt)

stats.filtered <- stats.filtered %>% tidyr::pivot_longer(names_to = "gene", values_to = "freq", cols = c("SF3B1.freq", "NRAS.freq"))
stats.filtered$sample.gene <- factor(paste0(stats.filtered$sample, ".", stats.filtered$gene),
  levels = c(
    "1026.1.NRAS.freq", "1026.3.NRAS.freq",
    "1026.1.SF3B1.freq", "1026.3.SF3B1.freq"
  )
)

# visualize for cells with sequencing depth >10x
ggplot(
  stats.filtered[which(stats.filtered$manual.cluster %ni% c("CD4", "CD8", "Plasma")), ],
  aes(x = sample.gene, y = 100 * freq, color = manual.cluster)
) +
  geom_point() +
  geom_line(aes(group = manual.cluster)) +
  scale_x_discrete(labels = c("Screening", "Treatment", "Screening", "Treatment")) +
  scale_y_continuous("% recipient positive", limits = c(0, 100)) +
  scale_color_manual(values = cluster.colors) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black"),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )
ggsave("./figure_coevolution/AML/figures/20230428_SF3B1_NRAS_filtered.svg", width = 1.3, height = 1.5)
