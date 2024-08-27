library(cowplot)
library(dplyr)
library(ggplot2)
library(ggsignif)
library(scales)
library(tidyr)
library(SummarizedExperiment)

meta <- readxl::read_excel("./data/CRC/20231127_CRC_mtDNA_meta.xlsx")

### read all data
mgatk.data.RNA <- 0
for (library in list.files("./data/CRC/", pattern = "*.mgatk")) {
  boo <- readRDS(paste0("./data/CRC/", library, "/final/mgatk.rds"))
  if (class(mgatk.data.RNA) == "numeric") {
    mgatk.data.RNA <- boo
  } else {
    mgatk.data.RNA <- cbind(mgatk.data.RNA, boo)
  }
}
colnames(mgatk.data.RNA) <- gsub(pattern = "Tumor-", replacement = "", colnames(mgatk.data.RNA))
colnames(mgatk.data.RNA) <- gsub(pattern = "CLL-CRC-", replacement = "", colnames(mgatk.data.RNA))

# provide lower limit of detection for variant outside of indicated patient
detection_limit <- function(variant, patient) {
  pos <- as.numeric(gsub("\\D+", "", variant))
  REF <- stringr::str_split_fixed(gsub("\\d", "", variant), pattern = ">", n = 2)[, 1]
  ALT <- stringr::str_split_fixed(gsub("\\d", "", variant), pattern = ">", n = 2)[, 2]

  samples <- meta$Sample[which(meta$Patient != patient)]

  refs <- paste0(REF, c("_counts_fw", "_counts_rev"))
  alts <- paste0(ALT, c("_counts_fw", "_counts_rev"))

  REF.reads <-
    assays(mgatk.data.RNA)[[refs[1]]][as.numeric(pos), samples] +
    assays(mgatk.data.RNA)[[refs[2]]][as.numeric(pos), samples]

  ALT.reads <-
    assays(mgatk.data.RNA)[[alts[1]]][as.numeric(pos), samples] +
    assays(mgatk.data.RNA)[[alts[2]]][as.numeric(pos), samples]

  sum(ALT.reads) / sum(REF.reads)
}

coverage_pos <- function(variant, patient) {
  pos <- as.numeric(gsub("\\D+", "", variant))
  REF <- stringr::str_split_fixed(gsub("\\d", "", variant), pattern = ">", n = 2)[, 1]
  ALT <- stringr::str_split_fixed(gsub("\\d", "", variant), pattern = ">", n = 2)[, 2]

  mean(assays(mgatk.data.RNA)[["coverage"]][as.numeric(pos), meta$Sample])
}

### detection limit aggregated
detection.limit.all <- data.frame()
for (p in unique(meta$Patient)) {
  if (p == "8") {
    next()
  }
  mtDNA.df <- read.table(paste0("./data/CRC/processed_together/", p, "_de_novo.csv"),
    check.names = F
  ) %>%
    as.data.frame()
  rownames(mtDNA.df) <- mtDNA.df$variant

  variants <- mtDNA.df$variant

  if (isEmpty(variants)) {
    next()
  }

  detection.limit.df <- data.frame(
    variant = variants,
    detection.limit = sapply(variants, FUN = function(x) {
      detection_limit(x, p)
    }),
    cov = sapply(variants, FUN = function(x) {
      coverage_pos(x, p)
    }),
    patient = p
  )

  if (nrow(detection.limit.df) > 0) {
    detection.limit.all <- rbind(detection.limit.all, detection.limit.df)
  }
}

ggplot(detection.limit.all, aes(x = cov, y = 100 * detection.limit)) +
  geom_smooth(method = "lm", color = "firebrick") +
  geom_point(size = 0.5) +
  scale_x_continuous("mean coverage at position") +
  scale_y_log10("% heteroplasmy", breaks = c(0.01, 0.1, 1, 10, 100), labels = c(0.01, 0.1, 1, 10, 100)) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title = element_text("Arial", size = 10, color = "black"),
    axis.text = element_text("Arial", size = 10, color = "black")
  )
ggsave("./figure_CLL_WES_RNA/figures/CRC/detection_limit/20240527_detection_limit_coverage.svg", width = 2, height = 2)

detection.limit.all$pos <- as.numeric(gsub("\\D+", "", detection.limit.all$variant))

GenePos.tib <- tibble(
  Names = c(
    "ATP6", "ATP8", "COX1", "COX2", "COX3", "CYTB", "ND1", "ND2", "ND3",
    "ND4", "ND4L", "ND5", "ND6", "RNR1", "RNR2"
  ),
  start = c(8527, 8366, 5904, 7586, 9207, 14747, 3307, 4470, 10059, 10760, 10470, 12337, 14149, 648, 1671),
  end = c(9207, 8572, 7445, 8269, 9990, 15887, 4262, 5511, 10404, 12137, 10766, 14148, 14673, 1601, 3229)
)
GenePos.tib$Names <- factor(GenePos.tib$Names, levels = GenePos.tib$Names)
y.max <- 100 * max(detection.limit.all$detection.limit)
GenePos.tib <- GenePos.tib %>%
  arrange(start) %>%
  mutate(mid = round((end - start) / 2 + start, 0), ycoord = rep(c(
    5 * y.max,
    3 * y.max,
    1 * y.max
  ), length.out = 15))

ggplot(detection.limit.all, aes(x = pos, y = 100 * detection.limit)) +
  geom_line(color = "#89CFF0") +
  geom_point(size = 0.5) +
  scale_x_continuous("position on chrM") +
  scale_y_log10("% heteroplasmy", breaks = c(0.01, 0.1, 1, 10, 100), labels = c(0.01, 0.1, 1, 10, 100)) +
  geom_segment(data = GenePos.tib, aes(x = start, y = ycoord, xend = end, yend = ycoord)) +
  geom_text(data = GenePos.tib, aes(x = mid, y = 1.4 * ycoord, label = Names), size = 2) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title = element_text("Arial", size = 10, color = "black"),
    axis.text = element_text("Arial", size = 10, color = "black")
  )
ggsave("./figure_CLL_WES_RNA/figures/CRC/detection_limit/20240528_detection_limit_position.svg", width = 3, height = 2)

cor.test(detection.limit.all$cov, detection.limit.all$detection.limit)
summary(detection.limit.all$detection.limit)

### for each case plot mtDNA mutations with most extreme changes
for (p in unique(meta$Patient)) {
  if (p == "8") {
    next()
  }
  mtDNA.df <- read.table(paste0("./data/CRC/processed_together/", p, "_de_novo.csv"),
    check.names = F
  ) %>%
    as.data.frame()
  rownames(mtDNA.df) <- mtDNA.df$variant

  variants <- mtDNA.df$variant[which(mtDNA.df$significant == "significant")]
  # exclude variants that are not detectable at either end
  variants <- variants[which(mtDNA.df[variants, "change"] != 0 & mtDNA.df[variants, "change"] != Inf)]

  # variants with change >10 or <0.1
  variants <- variants[which(mtDNA.df[variants, "change"] > 10 | mtDNA.df[variants, "change"] < 0.1)]

  detection.limit.df <- data.frame(
    variant = variants,
    detection.limit = sapply(variants, FUN = function(x) {
      detection_limit(x, p)
    })
  )

  # filter variants with detection limit >0.1%
  variants <- variants[which(detection.limit.df[variants, "detection.limit"] < 0.001)]
  if (isEmpty(variants)) {
    next()
  }

  detection.limit.df <- detection.limit.df[which(detection.limit.df$detection.limit < 0.001), ]

  timepoints <- colnames(mtDNA.df)[grepl("heteroplasmy", colnames(mtDNA.df))]
  mtDNA.plot.data <- mtDNA.df %>%
    dplyr::select(c("variant", timepoints)) %>%
    dplyr::filter(variant %in% variants) %>%
    pivot_longer(cols = timepoints)
  mtDNA.plot.data$day <- sapply(mtDNA.plot.data$name, function(x) {
    meta$Day[which(paste0(meta$Timepoint, ".heteroplasmy") == x & meta$Patient == p)]
  })
  mtDNA.plot.data$value.log <- mtDNA.plot.data$value
  mtDNA.plot.data$value.log[which(mtDNA.plot.data$value < 0.0001)] <- 0.0001

  color.scale <- as.character(BuenColors::jdb_palette(name = "corona", n = length(variants)))
  names(color.scale) <- variants
  color.scale <- c(color.scale, c(
    "pretreatment" = "blue", "posttreatment" = "red",
    "posttreatment2" = "firebrick"
  ))

  mtDNA.plot <- ggplot() +
    geom_hline(
      yintercept = 100 * as.numeric(detection.limit.df[variants, "detection.limit"]),
      color = color.scale[variants], linetype = "dashed"
    ) +
    geom_line(data = mtDNA.plot.data, aes(
      x = day, y = 100 * value.log,
      group = variant, color = as.character(variant)
    )) +
    geom_point(data = mtDNA.plot.data, aes(
      x = day, y = 100 * value.log,
      group = variant, color = as.character(variant)
    )) +
    geom_point(data = meta[which(meta$Patient == p), ], aes(x = Day, y = 250, color = Description), size = 2, shape = 18) +
    scale_y_log10("% heteroplasmy",
      limits = c(0.01, 250), breaks = c(0.01, 0.1, 1, 10, 100),
      labels = c("0.01", "0.1", "1", "10", "100")
    ) +
    scale_color_manual(values = color.scale) +
    ggtitle(p) +
    theme_classic() +
    theme(
      legend.position = "none",
      axis.title = element_text("Arial", size = 10, color = "black"),
      axis.text = element_text("Arial", size = 10, color = "black"),
      axis.title.x = element_blank(),
      plot.title = element_text("Arial", size = 10, color = "black", hjust = 0.5)
    )
  ggsave(paste0("./figure_CLL_WES_RNA/figures/CRC/detection_limit/20240103_", p, ".svg"),
    width = 2, height = 2, plot = mtDNA.plot
  )
}
