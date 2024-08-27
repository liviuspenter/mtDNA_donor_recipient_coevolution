library(dplyr)
library(ggplot2)
library(ggsignif)
library(scales)
library(tidyr)

meta <- readxl::read_excel("./data/CLL_RIC/20231229_CLL_RIC_mtDNA_meta.xlsx") %>%
  as.data.frame()
meta2 <- readxl::read_excel("./data/CLL_RIC/20231229_CLL_RIC_mtDNA_meta.xlsx", sheet = 2) %>%
  as.data.frame()
rownames(meta2) <- meta2$Patient
meta$Day2 <- 0
meta[which(meta$Timepoint == 2), "Day2"] <- meta2[as.character(meta$Patient[which(meta$Timepoint == 2)]), "Days_between_samples"]

for (p in unique(meta$Patient)) {
  variant.df <- read.table(paste0("./data/CLL_RIC/processed_together/", p, "_de_novo.csv"),
    check.names = F
  ) %>%
    as.data.frame()

  timepoints <- colnames(variant.df)[grepl("heteroplasmy", colnames(variant.df))]
  boo <- variant.df %>%
    select(c("variant", timepoints)) %>%
    as.data.frame()

  boo[is.na(boo)] <- 0

  rownames(boo) <- boo$variant
  boo <- boo[, -1]

  variant.df$variance <- rowVars(as.matrix(boo))

  boo <- variant.df %>%
    select(c("variant", timepoints, variance, significant, direction)) %>%
    pivot_longer(cols = timepoints)

  boo$day <- sapply(boo$name, function(x) {
    meta$Day2[which(paste0(meta$Timepoint, ".heteroplasmy") == x & meta$Patient == p)]
  })

  boo$value.log <- boo$value
  boo$value.log[which(boo$value == 0)] <- 0.0001
  ggplot() +
    geom_line(data = boo[which(boo$significant == "non.significant"), ], aes(
      x = day, y = 100 * value.log,
      group = variant
    ), color = "lightgrey") +
    geom_line(data = boo[which(boo$significant == "significant"), ], aes(
      x = day, y = 100 * value.log,
      group = variant, color = direction
    )) +
    geom_point(data = meta[which(meta$Patient == p), ], aes(x = Day, y = 250, color = Description), size = 2) +
    scale_y_log10("% heteroplasmy",
      limits = c(0.001, 260), breaks = c(0.01, 0.1, 1, 10, 100),
      labels = c("ND", "0.1", "1", "10", "100")
    ) +
    scale_color_manual(values = c(
      "pretreatment" = "blue", "posttreatment" = "red",
      "posttreatment2" = "firebrick",
      "increase" = "purple", "decrease" = "darkgreen",
      BuenColors::jdb_palette(name = "corona")
    )) +
    ggtitle(label = p) +
    theme_classic() +
    theme(
      legend.position = "none",
      axis.title = element_text("Arial", size = 10, color = "black"),
      axis.text = element_text("Arial", size = 10, color = "black"),
      axis.title.x = element_blank(),
      plot.title = element_text("Arial", size = 10, color = "black", hjust = 0.5)
    )
  ggsave(paste0("./figure_CLL_WES_RNA/figures/CLL_RIC/individual_plots_together/", p, "_log.svg"),
    width = 2, height = 2
  )
}

# CCF
CCF.data.all <- readxl::read_excel("./data/CLL_RIC/MAF/20231229_CLL_RIC_VAFs.xlsx") %>%
  as.data.frame()
CCF.data.all$patient <- stringr::str_split_fixed(CCF.data.all$Tumor_Sample_Barcode, pattern = "-", n = 3)[, 1]

for (p in unique(meta$Patient)) {
  CCF.data <- CCF.data.all[which(CCF.data.all$patient == p), ]

  CCF.data$identifier <- paste0(
    CCF.data$Chromosome, ":", CCF.data$Start_position, "_",
    CCF.data$Reference_Allele, ">", CCF.data$Tumor_Seq_Allele1
  )
  CCF.data$VAF <- CCF.data$validation_tumor_alt_count_wex / (CCF.data$validation_tumor_ref_count_wex + CCF.data$validation_tumor_alt_count_wex)
  CCF.data$Day <- sapply(CCF.data$Tumor_Sample_Barcode, function(x) {
    meta$Day2[which(meta$Sample == x)]
  })

  sample.first <- meta$Sample[which(meta$Patient == p)][1]
  sample.last <- last(meta$Sample[which(meta$Patient == p)])

  CCF.data$direction <- sapply(CCF.data$identifier, FUN = function(x) {
    VAF.first <- CCF.data$VAF[which(CCF.data$identifier == x & CCF.data$Tumor_Sample_Barcode == sample.first)]
    VAF.last <- CCF.data$VAF[which(CCF.data$identifier == x & CCF.data$Tumor_Sample_Barcode == sample.last)]
    ifelse(VAF.first > VAF.last, "decrease", "increase")
  })

  CCF.data$p.value <- sapply(CCF.data$identifier, FUN = function(x) {
    REF.first <- CCF.data$validation_tumor_ref_count_wex[which(CCF.data$identifier == x & CCF.data$Tumor_Sample_Barcode == sample.first)]
    REF.last <- CCF.data$validation_tumor_ref_count_wex[which(CCF.data$identifier == x & CCF.data$Tumor_Sample_Barcode == sample.last)]
    ALT.first <- CCF.data$validation_tumor_alt_count_wex[which(CCF.data$identifier == x & CCF.data$Tumor_Sample_Barcode == sample.first)]
    ALT.last <- CCF.data$validation_tumor_alt_count_wex[which(CCF.data$identifier == x & CCF.data$Tumor_Sample_Barcode == sample.last)]

    REF.first <- ifelse(isEmpty(REF.first), 0, REF.first)
    REF.last <- ifelse(isEmpty(REF.last), 0, REF.last)
    ALT.first <- ifelse(isEmpty(ALT.first), 0, ALT.first)
    ALT.last <- ifelse(isEmpty(ALT.last), 0, ALT.last)

    fisher.test(matrix(c(REF.first, REF.last, ALT.first, ALT.last), nrow = 2))$p.value
  })

  CCF.data$significant <- ifelse(CCF.data$p.value < 0.05, "significant", "non.significant")

  ggplot() +
    geom_line(
      data = CCF.data[which(CCF.data$significant == "non.significant"), ],
      aes(x = Day, y = VAF, group = identifier), color = "lightgrey"
    ) +
    geom_line(
      data = CCF.data[which(CCF.data$significant == "significant" & !isEmpty(CCF.data$direction)), ],
      aes(x = Day, y = VAF, group = identifier, color = as.character(direction))
    ) +
    scale_y_continuous(limits = c(0, 1)) +
    scale_color_manual(values = c("decrease" = "darkgreen", "increase" = "purple")) +
    theme_classic() +
    theme(
      legend.position = "none",
      axis.title = element_text("Arial", size = 10, color = "black"),
      axis.text = element_text("Arial", size = 10, color = "black"),
      axis.title.x = element_blank(),
      plot.title = element_text("Arial", size = 10, color = "black", hjust = 0.5)
    )
  ggsave(paste0("./figure_CLL_WES_RNA/figures/CLL_RIC//individual_plots_together/", p, "_ccf.svg"),
    width = 2, height = 2
  )
}
