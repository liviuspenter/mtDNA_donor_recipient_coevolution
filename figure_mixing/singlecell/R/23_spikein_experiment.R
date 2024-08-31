library(ArchR)
library(dplyr)
library(ggplot2)
library(parallel)
library(Seurat)
library(Signac)

AML.mix <- loadArchRProject("./data/mixing/singlecell/AML.mix/")

informative.variants <- as.data.frame(read.csv2(file = "./data/mixing/singlecell/artificial_mixing_AML/informative_variants.csv", sep = "\t"))

# read out mixing experiment
df <- data.frame()
for (exp in c(1, 5, 10, 50, 100, 500, 1000)) {
  message(exp)
  # load mgatk output
  mito.data <- ReadMGATK(dir = paste0("./data/mixing/singlecell/artificial_mixing_AML/mtscATAC-seq/mixing_", exp, ".mgatk/"))

  # read out informative mtDNA mutations
  VAFs <- AlleleFreq(mito.data$counts, variants = informative.variants$variant, assay = NULL)

  df <- rbind(df, data.frame(
    bc = colnames(VAFs),
    VAFs.1 = colMeans(VAFs[informative.variants$variant[which(informative.variants$individual == "AML1011")], ]),
    VAFs.2 = colMeans(VAFs[informative.variants$variant[which(informative.variants$individual == "AML1012")], ]),
    spike.in = exp
  ))
}
write.table(df, quote = F, sep = "\t", row.names = F, file = "./data/mixing/singlecell/artificial_mixing_AML/mtscATAC-seq/results.csv")
df <- read.table(file = "./data/mixing/singlecell/artificial_mixing_AML/mtscATAC-seq/results.csv", header = T)

# read ground-truth data
df$sample.orig <- "none"
for (exp in c(1, 5, 10, 50, 100, 500, 1000)) {
  barcodes.1 <- as.data.frame(data.table::fread(file = paste0("./data/mixing/singlecell/artificial_mixing_AML/barcodes_mtscATACseq/AML1011_barcodes.", exp), header = F))
  barcodes.2 <- as.data.frame(data.table::fread(file = paste0("./data/mixing/singlecell/artificial_mixing_AML/barcodes_mtscATACseq/AML1012_barcodes.", exp), header = F))
  df[which(df$spike.in == exp & df$bc %in% barcodes.1$V1), "sample.orig"] <- "AML1011"
  df[which(df$spike.in == exp & df$bc %in% barcodes.2$V1), "sample.orig"] <- "AML1012"
}

# annotate based on mtDNA
df$sample.annotation <- "none"
df[which(df$VAFs.1 > 0.8 & df$VAFs.2 < 0.2), "sample.annotation"] <- "AML1011"
df[which(df$VAFs.1 < 0.8 & df$VAFs.2 > 0.8), "sample.annotation"] <- "AML1012"

# assess annotation
df$result <- apply(df, MARGIN = 1, FUN = function(x) {
  if (x["sample.annotation"] == "none") {
    "none"
  } else if (x["sample.orig"] == x["sample.annotation"]) {
    "correct"
  } else if (x["sample.orig"] != x["sample.annotation"]) {
    "mismatch"
  }
})

results <- df %>%
  group_by(spike.in) %>%
  summarize(
    none.1 = length(which(sample.orig == "AML1011" & result == "none")),
    none.2 = length(which(sample.orig == "AML1012" & result == "none")),
    correct.1 = length(which(sample.orig == "AML1011" & result == "correct")),
    correct.2 = length(which(sample.orig == "AML1012" & result == "correct"))
  ) %>%
  as.data.frame()

ggplot() +
  geom_abline(slope = 1, color = "lightgrey") +
  geom_point(data = results, aes(x = spike.in, y = correct.1), size = 0.5, color = "orange") +
  geom_point(data = results, aes(x = spike.in, y = correct.2), size = 0.5, color = "purple") +
  scale_x_log10("cells spiked-in") +
  scale_y_log10("correct annotation") +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black")
  )
ggsave("./figure_mixing/singlecell/figures/20240603_spikein_mtscATAC.svg", width = 1.5, height = 1.5)

ggplot(results, aes(x = spike.in, y = none.1)) +
  geom_abline(slope = 1) +
  geom_point(size = 0.5) +
  scale_x_log10("cells spiked-in", limits = c(1, 1000)) +
  scale_y_log10("missing annotation", limits = c(1, 1000)) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black")
  )
ggsave("./figure_mixing/singlecell/figures/20240603_spikein_mtscATAC_unannotated.svg", width = 1.5, height = 1.5)
