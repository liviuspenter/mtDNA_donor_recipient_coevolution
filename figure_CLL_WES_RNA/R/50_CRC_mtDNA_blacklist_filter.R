# identify mtDNA variants that will be considered artifacts
# criteria:

library(dplyr)
library(ggplot2)
library(tidyr)

source("./figure_mixing/bulk/R/00_variant_calling_atac.R")

### compare RNA with ATAC-seq data from Vago dataset
matched.files <- readxl::read_excel("./data/Vago/20231128_paired_RNA_ATAC_seq.xlsx")
matched.files$sample <- seq(1, nrow(matched.files))

### read all data
mgatk.data.ATAC <- 0
for (library in c(
  paste0("./data/Vago/PRJNA810181/", matched.files$RNA, ".mgatk"),
  paste0("./data/Vago/PRJNA810182/", matched.files$ATAC, ".mgatk")
)) {
  boo <- readRDS(paste0(library, "/final/mgatk.rds"))
  if (class(mgatk.data.ATAC) == "numeric") {
    mgatk.data.ATAC <- boo
  } else {
    mgatk.data.ATAC <- cbind(mgatk.data.ATAC, boo)
  }
}

### read variants
results.df <- data.frame()
for (library in c(
  paste0("./data/Vago/PRJNA810181/", matched.files$RNA, ".mgatk"),
  paste0("./data/Vago/PRJNA810182/", matched.files$ATAC, ".mgatk")
)) {
  message(library)
  results <- variant_calling_bulk(
    sample.bulk = readRDS(paste0(library, "/final/mgatk.rds")),
    coverage.position = 10, strand.coordination = 0.5, total.coverage = 100,
    sample.name = s
  )
  results <- results[which(!is.na(results$heteroplasmy)), ]

  if (nrow(results) > 0) {
    results$sample <- stringr::str_split_fixed(library, pattern = "/", n = 5)[, 5]
    results.df <- rbind(results.df, results)
  }
}
results.df$run <- stringr::str_split_fixed(results.df$sample, pattern = "\\.", n = 2)[, 1]

matched.files <- matched.files %>%
  pivot_longer(cols = c("RNA", "ATAC")) %>%
  as.data.frame()
rownames(matched.files) <- matched.files$value
results.df$sample <- matched.files[results.df$run, "sample"]
results.df$type <- matched.files[results.df$run, "name"]

boo <-
  results.df %>% pivot_wider(
    id_cols = c("sample", "variant"),
    names_from = c("type"),
    values_from = c("coverage", "heteroplasmy")
  )
boo$match <- ifelse(!is.na(boo$heteroplasmy_RNA) & !is.na(boo$heteroplasmy_ATAC), "yes", "no")
# exclude variants not covered in RNA data
boo <- boo[-which(is.na(boo$heteroplasmy_RNA) & !is.na(boo$heteroplasmy_ATAC)), ]

stats <-
  boo %>%
  dplyr::count(variant, match) %>%
  pivot_wider(names_from = c("match"), values_from = "n")

blacklist <- stats[which(stats$no > 2 & is.na(stats$yes)), ]

blacklist$position <- rank(-blacklist$no, ties.method = "random")

ggplot(blacklist, aes(x = position, y = no)) +
  geom_col() +
  scale_x_continuous("mtDNA variant") +
  scale_y_continuous("samples") +
  theme_classic() +
  theme(
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black")
  )
ggsave("./figure_CLL_WES_RNA/figures/blacklist/20240103_blacklist_raw.svg", width = 2, height = 1.5)

write.table(blacklist, file = "./data/Vago/20240103_blacklist_raw_Gambacorta.csv", quote = F, sep = "\t")

### read all data
mgatk.data.ATAC <- 0
for (library in list.files("./data/CRC/", pattern = "*.mgatk")) {
  boo <- readRDS(paste0("./data/CRC/", library, "/final/mgatk.rds"))
  if (class(mgatk.data.ATAC) == "numeric") {
    mgatk.data.ATAC <- boo
  } else {
    mgatk.data.ATAC <- cbind(mgatk.data.ATAC, boo)
  }
}

### read variants
results.df <- data.frame()
for (library in list.files("./data/CRC/", pattern = "*.mgatk")) {
  message(library)
  results <- variant_calling_bulk(
    sample.bulk = readRDS(paste0("./data/CRC/", library, "/final/mgatk.rds")),
    coverage.position = 10, strand.coordination = 0.5, total.coverage = 100,
    sample.name = s
  )
  results <- results[which(!is.na(results$heteroplasmy)), ]

  if (nrow(results) > 0) {
    results$sample <- library
    results.df <- rbind(results.df, results)
  }
}
results.df$run <- stringr::str_split_fixed(results.df$sample, pattern = "\\.", n = 2)[, 1]
results.df$patient <- stringr::str_split_fixed(results.df$run, pattern = "-", n = 3)[, 2]
results.df$timepoint <- stringr::str_split_fixed(results.df$run, pattern = "-", n = 3)[, 3]

# filtering statistics
length(unique(results.df$variant[which(results.df$patient != "0008")]))

length(unique(results.df$variant[which(results.df$patient != "0008" &
  results.df$heteroplasmy > 0.005 &
  results.df$heteroplasmy < 0.95)]))

df <- sort(table(results.df$variant[which(results.df$heteroplasmy > 0.005 &
  results.df$heteroplasmy < 0.95)]), decreasing = T) %>% as.data.frame()
colnames(df) <- c("variant", "n")
df$freq <- df$n / length(unique(results.df$sample))
df$blacklist <- ifelse(df$variant %in% blacklist$variant & df$freq > 0.2, "yes", "no")

df$position <- rank(-df$n, ties.method = "random")

ggplot(df, aes(x = position, y = 100 * freq, fill = blacklist)) +
  geom_col() +
  geom_hline(yintercept = 20, color = "firebrick") +
  scale_x_continuous("mtDNA variant") +
  scale_y_continuous("% samples") +
  scale_fill_manual(values = c("yes" = "black", "no" = "grey")) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black")
  )
ggsave("./figure_CLL_WES_RNA/figures/blacklist/20240103_blacklist_CRC.svg", width = 2, height = 1.5)

write.table(df, file = "./data/CRC/processed/20231128_blacklist.csv", quote = F, sep = "\t")
