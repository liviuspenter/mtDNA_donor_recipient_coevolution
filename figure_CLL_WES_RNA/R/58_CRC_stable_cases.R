library(cowplot)
library(dplyr)
library(ggplot2)
library(ggsignif)
library(scales)
library(tidyr)

meta <- readxl::read_excel("./data/CRC/20231127_CRC_mtDNA_meta.xlsx")

stable.cases <- c(2, 14, 15, 17, 18, 36, 40, 41)

mtDNA.dynamics <- data.frame()
CCF.dynamics <- data.frame()
for (p in unique(meta$Patient)) {
  if (p == "8") {
    next()
  }

  mtDNA.df <- read.table(paste0("./data/CRC/processed_together/", p, "_de_novo.csv"),
    check.names = F
  ) %>%
    as.data.frame()

  # remove variants that cannot be compared
  mtDNA.df <- mtDNA.df[which(!mtDNA.df$change %in% c(0, Inf)), ]

  # calculate mtDNA coefficient of variation
  columns <- colnames(mtDNA.df)[grepl("heteroplasmy", colnames(mtDNA.df))]
  mtDNA.df$sd <- apply(mtDNA.df[, columns], MARGIN = 1, FUN = sd)
  mtDNA.df$mean <- apply(mtDNA.df[, columns], MARGIN = 1, FUN = mean)
  mtDNA.df$variance <- apply(mtDNA.df[, columns], MARGIN = 1, FUN = var)
  mtDNA.df$cv <- 100 * (mtDNA.df$sd / mtDNA.df$mean)
  mtDNA.df$vmr <- mtDNA.df$variance / mtDNA.df$mean

  # calculate CCF coefficient of variation
  CCF.df <- data.table::fread(paste0("./data/CRC/MAF/edited/", p, ".mut_ccfs.txt"))
  CCF.df <- CCF.df[which(CCF.df$Sample_ID %in% meta$Sample[which(meta$Patient == p)]), ]

  # pull cluster ccf mean values
  cluster.data <- CCF.df %>%
    group_by(Sample_ID, Cluster_Assignment) %>%
    summarize(clust_ccf_mean = unique(clust_ccf_mean)) %>%
    pivot_wider(names_from = Sample_ID, values_from = clust_ccf_mean)

  cluster.data$sd <- apply(cluster.data[, seq(2, ncol(cluster.data))], MARGIN = 1, FUN = sd)
  cluster.data$mean <- apply(cluster.data[, seq(2, ncol(cluster.data))], MARGIN = 1, FUN = mean)
  cluster.data$variance <- apply(cluster.data[, seq(2, ncol(cluster.data))], MARGIN = 1, FUN = var)
  cluster.data$cv <- 100 * (cluster.data$sd / cluster.data$mean)
  cluster.data$vmr <- cluster.data$variance / cluster.data$mean

  mtDNA.df$patient <- p
  mtDNA.dynamics <- rbind(mtDNA.dynamics, mtDNA.df[, c("patient", "variant", "sd", "mean", "cv", "variance", "vmr")])

  cluster.data$patient <- p
  CCF.dynamics <- rbind(CCF.dynamics, cluster.data[, c("patient", "Cluster_Assignment", "sd", "mean", "cv", "variance", "vmr")])
}

mtDNA.dynamics$category <- ifelse(mtDNA.dynamics$patient %in% stable.cases, "stable", "non.stable")
CCF.dynamics$category <- ifelse(CCF.dynamics$patient %in% stable.cases, "stable", "non.stable")

ggplot(CCF.dynamics, aes(x = cv)) +
  geom_histogram(aes(fill = category)) +
  scale_fill_manual(values = c("non.stable" = "#9e15c5", "stable" = "#00fc02")) +
  scale_x_continuous("coefficient of variance") +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title = element_text("Arial", size = 10, color = "black"),
    axis.text = element_text("Arial", size = 10, color = "black")
  )
ggsave("./figure_CLL_WES_RNA/figures/CRC/20240102_CV_CCF.svg", width = 2, height = 1.5)

t.test(CCF.dynamics$cv[which(CCF.dynamics$category == "stable")], CCF.dynamics$cv[which(CCF.dynamics$category == "non.stable")])

ggplot(mtDNA.dynamics, aes(x = cv)) +
  geom_histogram(aes(fill = category)) +
  scale_fill_manual(values = c("non.stable" = "#9e15c5", "stable" = "#00fc02")) +
  scale_x_continuous("coefficient of variance") +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title = element_text("Arial", size = 10, color = "black"),
    axis.text = element_text("Arial", size = 10, color = "black")
  )
ggsave("./figure_CLL_WES_RNA/figures/CRC/20240102_CV_mtDNA.svg", width = 2, height = 1.5)

t.test(
  mtDNA.dynamics$cv[which(mtDNA.dynamics$category == "stable")],
  mtDNA.dynamics$cv[which(mtDNA.dynamics$category == "non.stable")]
)
