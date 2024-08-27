library(dplyr)
library(ggplot2)
library(ggsignif)
library(scales)
library(tidyr)

meta <- readxl::read_excel("./data/CRC/20231127_CRC_mtDNA_meta.xlsx")

for (p in unique(meta$Patient)) {
  variant.df <- read.table(paste0("./data/CRC/processed_together/", p, "_de_novo.csv"),
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
    meta$Day[which(paste0(meta$Timepoint, ".heteroplasmy") == x & meta$Patient == p)]
  })

  boo$value.log <- boo$value
  boo$value.log[which(boo$value == 0)] <- 0.00001
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
      limits = c(0.001, 260), breaks = c(0.001, 0.01, 0.1, 1, 10, 100),
      labels = c("ND", "0.01", "0.1", "1", "10", "100")
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
  ggsave(paste0("./figure_CLL_WES_RNA/figures/CRC/individual_plots_together/", p, "_log.svg"),
    width = 2, height = 2
  )
}

# CCF
for (p in unique(meta$Patient)) {
  CCF.data <- data.table::fread(paste0("./data/CRC/MAF/edited/", p, ".mut_ccfs.txt"))

  CCF.data <- CCF.data[which(CCF.data$Sample_ID %in% meta$Sample[which(meta$Patient == p)]), ]

  CCF.data$identifier <- paste0(
    CCF.data$Chromosome, ":", CCF.data$Start_position, "_",
    CCF.data$Reference_Allele, ">", CCF.data$Tumor_Seq_Allele
  )
  CCF.data$VAF <- CCF.data$t_alt_count / (CCF.data$t_ref_count + CCF.data$t_alt_count)
  CCF.data$Day <- sapply(CCF.data$Sample_ID, function(x) {
    meta$Day[which(meta$Sample == x)]
  })

  sample.first <- meta$Sample[which(meta$Patient == p)][1]
  sample.last <- last(meta$Sample[which(meta$Patient == p)])

  CCF.data$direction <- sapply(CCF.data$identifier, FUN = function(x) {
    VAF.first <- CCF.data$VAF[which(CCF.data$identifier == x & CCF.data$Sample_ID == sample.first)]
    VAF.last <- CCF.data$VAF[which(CCF.data$identifier == x & CCF.data$Sample_ID == sample.last)]
    ifelse(VAF.first > VAF.last, "decrease", "increase")
  })

  CCF.data$p.value <- sapply(CCF.data$identifier, FUN = function(x) {
    REF.first <- CCF.data$t_ref_count[which(CCF.data$identifier == x & CCF.data$Sample_ID == sample.first)]
    REF.last <- CCF.data$t_ref_count[which(CCF.data$identifier == x & CCF.data$Sample_ID == sample.last)]
    ALT.first <- CCF.data$t_alt_count[which(CCF.data$identifier == x & CCF.data$Sample_ID == sample.first)]
    ALT.last <- CCF.data$t_alt_count[which(CCF.data$identifier == x & CCF.data$Sample_ID == sample.last)]

    fisher.test(matrix(c(REF.first, REF.last, ALT.first, ALT.last), nrow = 2))$p.value
  })

  CCF.data$significant <- ifelse(CCF.data$p.value < 0.05, "significant", "non.significant")

  ggplot() +
    geom_line(
      data = CCF.data[which(CCF.data$significant == "non.significant"), ],
      aes(x = Day, y = VAF, group = identifier), color = "lightgrey"
    ) +
    geom_line(
      data = CCF.data[which(CCF.data$significant == "significant"), ],
      aes(x = Day, y = VAF, group = identifier, color = direction)
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
  ggsave(paste0("./figure_CLL_WES_RNA//figures/CRC/individual_plots_together/", p, "_ccf.svg"),
    width = 2, height = 2
  )
}
