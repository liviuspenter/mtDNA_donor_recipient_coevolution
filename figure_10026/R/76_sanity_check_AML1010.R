library(ArchR)

AML1010.mito <- loadArchRProject("./data/10026/AML.1010.mito/")
combined.frequences <- readRDS("./data/10026/mtDNA/20210429_AML1010_combined_mutation_frequencies.rds")

AML1010.1 <- ReadMGATK("./data/10026/AML1010_1.mgatk/")

AML1010.variants <- data.table::fread("./data/10026/mtDNA/20210429_AML1010_high_confidence_variants.csv")

df.all <- data.frame()

GC.variants <- AML1010.variants$variant[grepl("G>C", AML1010.variants$variant)]
GA.variants <- AML1010.variants$variant[grepl("G>A", AML1010.variants$variant)]

for (variant in GC.variants) {
  message(variant)
  pos <- gsub(variant, pattern = "G>C", replacement = "")
  df <- data.frame(cb = AML1010.mito$cellNames[which(AML1010.mito$Sample == "AML1010_1")])
  df$cb.simple <- stringr::str_split_fixed(df$cb, pattern = "#", n = 2)[, 2]
  df$coverage <- AML1010.1$counts[paste0("C-", pos, "-fwd"), df$cb.simple] +
    AML1010.1$counts[paste0("C-", pos, "-rev"), df$cb.simple]
  df$heteroplasmy <- as.character(combined.frequences[variant, df$cb])
  df$variant <- variant
  df$type <- "GC"
  df.all <- rbind(df.all, df)
}

for (variant in GA.variants) {
  message(variant)
  pos <- gsub(variant, pattern = "G>A", replacement = "")
  df <- data.frame(cb = AML1010.mito$cellNames[which(AML1010.mito$Sample == "AML1010_1")])
  df$cb.simple <- stringr::str_split_fixed(df$cb, pattern = "#", n = 2)[, 2]
  df$coverage <- AML1010.1$counts[paste0("A-", pos, "-fwd"), df$cb.simple] +
    AML1010.1$counts[paste0("A-", pos, "-rev"), df$cb.simple]
  df$heteroplasmy <- as.character(combined.frequences[variant, df$cb])
  df$variant <- variant
  df$type <- "GA"
  df.all <- rbind(df.all, df)
}

p <- ggplot() +
  ggrastr::rasterize(geom_point(data = df.all, aes(x = 100 * as.numeric(heteroplasmy), y = coverage, color = type), size = 0.5), dpi = 600) +
  scale_color_manual(values = c("GC" = "firebrick", "GA" = "grey")) +
  scale_x_log10("% single cell heteroplasmy of variant") +
  scale_y_log10("single cell coverage of variant") +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black")
  )
ggsave("./figure_10026/figures/qc/20240625_AML1010.1_heteroplasmy_coverage.svg", width = 3, height = 3, plot = p)
ggsave("./figure_10026/figures/qc/20240625_AML1010.1_heteroplasmy_coverage.png", width = 3, height = 3, plot = p)
