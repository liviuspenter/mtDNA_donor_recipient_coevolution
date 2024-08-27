library(dplyr)

source("./figure_mixing/bulk/R/00_variant_calling_atac.R")

# blacklist
mtDNA.blacklist <- read.table(file = "./data/CRC/processed/20231128_blacklist.csv", sep = "\t") %>%
  as.data.frame()
exclude.variants <- mtDNA.blacklist$variant[which(mtDNA.blacklist$blacklist == "yes")]

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

vaf <- results.df %>%
  tidyr::pivot_wider(names_from = variant, values_from = heteroplasmy, id_cols = sample) %>%
  as.data.frame()

###
### gather information per patient across different identifiers
###

meta <- readxl::read_excel("./data/CRC/20231127_CRC_mtDNA_meta.xlsx")

for (p in unique(meta$Patient)) {
  files <- paste0("./data/CRC/CRC-", meta$Sample[which(meta$Patient == p)], ".mgatk")

  results.df <- data.frame()
  for (library in files) {
    message(library)
    results <- variant_calling_bulk(
      sample.bulk = readRDS(paste0(library, "/final/mgatk.rds")),
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

  variant.df.all <- data.frame()

  # make a list of variants
  variants <- unique(results.df$variant[which(results.df$heteroplasmy > 0.005 &
    results.df$heteroplasmy < 0.95)])

  # exclude blacklisted variants
  variants <- variants[which(!variants %in% exclude.variants)]

  variant.df <- data.frame(
    variant = variants,
    pos = as.numeric(gsub("\\D+", "", variants)),
    REF = stringr::str_split_fixed(gsub("\\d", "", variants), pattern = ">", n = 2)[, 1],
    ALT = stringr::str_split_fixed(gsub("\\d", "", variants), pattern = ">", n = 2)[, 2]
  )

  colnames(mgatk.data.ATAC) <- gsub(pattern = "Tumor-", replacement = "", colnames(mgatk.data.ATAC))
  colnames(mgatk.data.ATAC) <- gsub(pattern = "CLL-CRC-", replacement = "", colnames(mgatk.data.ATAC))

  for (s in meta$Sample[which(meta$Patient == p)]) {
    timepoint <- meta$Timepoint[which(meta$Sample == s)]
    # perform de-novo variant calling from mgatk data
    boo <- cbind(
      variant.df[, c(1, 2, 3, 4)],
      assays(mgatk.data.ATAC)[["coverage"]][
        variant.df$pos,
        s
      ] %>% as.matrix() %>% as.data.frame()
    )
    colnames(boo)[5] <- "coverage"

    boo$REF.reads <-
      apply(boo, MARGIN = 1, FUN = function(x) {
        refs <- paste0(x["REF"], c("_counts_fw", "_counts_rev"))
        assays(mgatk.data.ATAC)[[refs[1]]][as.numeric(x["pos"]), s] +
          assays(mgatk.data.ATAC)[[refs[2]]][as.numeric(x["pos"]), s]
      })
    boo$ALT.reads <-
      apply(variant.df, MARGIN = 1, FUN = function(x) {
        refs <- paste0(x["ALT"], c("_counts_fw", "_counts_rev"))
        assays(mgatk.data.ATAC)[[refs[1]]][as.numeric(x["pos"]), s] +
          assays(mgatk.data.ATAC)[[refs[2]]][as.numeric(x["pos"]), s]
      })

    boo$heteroplasmy <- boo$ALT.reads / boo$coverage
    colnames(boo)[c(5, 6, 7, 8)] <- paste0(timepoint, ".", colnames(boo)[c(5, 6, 7, 8)])
    variant.df <- cbind(variant.df, boo[, c(5, 6, 7, 8)])
  }

  last.sample <- length(meta$Sample[which(meta$Patient == p)])

  # perform statistical testing
  variant.df$p.value <- apply(variant.df, MARGIN = 1, FUN = function(x) {
    fisher.test(matrix(
      data = as.numeric(c(
        x["1.REF.reads"], x[paste0(last.sample, ".REF.reads")],
        x["1.ALT.reads"], x[paste0(last.sample, ".ALT.reads")]
      )),
      nrow = 2
    ))$p.value
  })

  variant.df$change <- apply(variant.df, MARGIN = 1, FUN = function(x) {
    as.numeric(x["1.heteroplasmy"]) / as.numeric(x[paste0(last.sample, ".heteroplasmy")])
  })

  variant.df$direction <- ifelse(variant.df$change > 1, "decrease", "increase")

  variant.df$significant <- ifelse(variant.df$p.value < 0.05 &
    (variant.df$change > 1.3 | variant.df$change < 0.7),
  "significant", "non.significant"
  )

  colnames(variant.df) <- gsub(colnames(variant.df), pattern = "-", replacement = "\\.")
  # save data
  write.table(variant.df,
    file = paste0("./data/CRC/processed_together/", p, "_de_novo.csv"),
    sep = "\t", quote = F
  )
}
