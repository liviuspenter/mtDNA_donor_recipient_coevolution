library(cowplot)
library(dplyr)
library(ggplot2)
library(ggsignif)
library(scales)
library(tidyr)

meta <- readxl::read_excel("./data/CRC/20231127_CRC_mtDNA_meta.xlsx")

dynamics.df <- data.frame()
cluster.stats <- data.frame()
for (p in unique(meta$Patient)) {
  if (p == "8") {
    next()
  }
  sample.first <- meta$Sample[which(meta$Patient == p)][1]
  sample.last <- last(meta$Sample[which(meta$Patient == p)])

  mtDNA.df <- read.table(paste0("./data/CRC/processed_together/", p, "_de_novo.csv"),
    check.names = F
  ) %>%
    as.data.frame()

  # remove variants that cannot be compared
  mtDNA.df <- mtDNA.df[which(!mtDNA.df$change %in% c(0, Inf)), ]

  CCF.df <- data.table::fread(paste0("./data/CRC/MAF/edited/", p, ".mut_ccfs.txt"))
  CCF.df <- CCF.df[which(CCF.df$Sample_ID %in% meta$Sample[which(meta$Patient == p)]), ]

  # pull cluster ccf mean values
  cluster.data <- CCF.df %>%
    group_by(Sample_ID, Cluster_Assignment) %>%
    summarize(clust_ccf_mean = unique(clust_ccf_mean)) %>%
    pivot_wider(names_from = Sample_ID, values_from = clust_ccf_mean)

  # calculate dynamics of clusters
  cluster.data$dynamics <- unlist(cluster.data[, sample.first] / cluster.data[, sample.last])

  # match mtDNA mutations
  mtDNA.df$cluster.match <- sapply(mtDNA.df$change, FUN = function(x) {
    which.min(abs(cluster.data$dynamics - x))
  })
  mtDNA.df$cluster.match.change <- unlist(cluster.data[mtDNA.df$cluster.match, "dynamics"])

  dynamics.df <- rbind(dynamics.df, data.frame(
    patient = p,
    variant = mtDNA.df$variant,
    cluster.match = mtDNA.df$cluster.match,
    mtDNA.change = mtDNA.df$change,
    CCF.change = mtDNA.df$cluster.match.change
  ))

  # plot CCF and matched mtDNA changes
  CCF.plot.data <- CCF.df %>%
    group_by(Sample_ID, Cluster_Assignment) %>%
    summarize(clust_ccf_mean = unique(clust_ccf_mean))
  CCF.plot.data$day <- sapply(CCF.plot.data$Sample_ID, function(x) {
    meta$Day[which(meta$Sample == x)]
  })

  timepoints <- colnames(mtDNA.df)[grepl("heteroplasmy", colnames(mtDNA.df))]
  mtDNA.plot.data <- mtDNA.df %>%
    select(c("variant", timepoints, cluster.match)) %>%
    pivot_longer(cols = timepoints)
  mtDNA.plot.data$day <- sapply(mtDNA.plot.data$name, function(x) {
    meta$Day[which(paste0(meta$Timepoint, ".heteroplasmy") == x & meta$Patient == p)]
  })
  mtDNA.plot.data$value.log <- mtDNA.plot.data$value
  mtDNA.plot.data$value.log[which(mtDNA.plot.data$value == 0)] <- 0.00001

  color.scale <- as.character(BuenColors::jdb_palette(name = "corona"))
  names(color.scale) <- as.character(seq(1, length(color.scale)))
  color.scale <- c(color.scale, c(
    "pretreatment" = "blue", "posttreatment" = "red",
    "posttreatment2" = "firebrick"
  ))

  CCF.plot <- ggplot() +
    geom_line(data = CCF.plot.data, aes(
      x = day, y = 100 * clust_ccf_mean,
      group = Cluster_Assignment, color = as.character(Cluster_Assignment)
    )) +
    geom_point(data = CCF.plot.data, aes(
      x = day, y = 100 * clust_ccf_mean,
      group = Cluster_Assignment, color = as.character(Cluster_Assignment)
    ), size = 0.5) +
    geom_point(data = meta[which(meta$Patient == p), ], aes(x = Day, y = 110, color = Description), size = 1) +
    scale_y_continuous("% CCF", breaks = c(0, 25, 50, 75, 100)) +
    scale_color_manual(values = color.scale) +
    theme_classic() +
    theme(
      legend.position = "none",
      axis.title = element_text("Arial", size = 10, color = "black"),
      axis.text = element_text("Arial", size = 10, color = "black"),
      axis.title.x = element_blank(),
      plot.title = element_text("Arial", size = 10, color = "black", hjust = 0.5)
    )

  mtDNA.plot <- ggplot() +
    geom_line(data = mtDNA.plot.data, aes(
      x = day, y = 100 * value.log,
      group = variant, color = as.character(cluster.match)
    )) +
    geom_point(data = mtDNA.plot.data, aes(
      x = day, y = 100 * value.log,
      group = variant, color = as.character(cluster.match)
    ), size = 0.5) +
    geom_point(data = meta[which(meta$Patient == p), ], aes(x = Day, y = 250, color = Description), size = 1) +
    scale_y_log10("% heteroplasmy",
      limits = c(0.001, 250), breaks = c(0.001, 0.01, 0.1, 1, 10, 100),
      labels = c("ND", "0.01", "0.1", "1", "10", "100")
    ) +
    scale_color_manual(values = color.scale) +
    theme_classic() +
    theme(
      legend.position = "none",
      axis.title = element_text("Arial", size = 10, color = "black"),
      axis.text = element_text("Arial", size = 10, color = "black"),
      axis.title.x = element_blank(),
      plot.title = element_text("Arial", size = 10, color = "black", hjust = 0.5)
    )

  # now add the title
  title <- ggdraw() +
    draw_label(
      p,
      size = 10,
      x = 0.5,
      hjust = 0.5
    )
  plot_row <- plot_grid(plotlist = list(CCF.plot, mtDNA.plot))
  integrated.plot <- plot_grid(title, plot_row, ncol = 1, rel_heights = c(0.05, 1))
  ggsave(paste0("./figure_CLL_WES_RNA/figures/CRC/integrated/20231130_", p, ".svg"),
    width = 3, height = 1.5, plot = integrated.plot
  )

  # how many CCF clusters have a match at mtDNA level
  cluster.stats <-
    rbind(
      cluster.stats,
      data.frame(
        patient = p,
        CCF.clusters = length(cluster.data$Cluster_Assignment),
        mtDNA.mutations = length(mtDNA.df$variant),
        CCF.clusters.matched = length(unique(mtDNA.df$cluster.match)),
        CCF.clusters.not.matched =
          length(setdiff(cluster.data$Cluster_Assignment, unique(mtDNA.df$cluster.match)))
      )
    )
}

# correlation of mtDNA and somatic nuclear changes
# cor.test(dynamics.df$mtDNA.change, dynamics.df$CCF.change)
# r = 0.84
# p < 2.2x10^-16
ggplot(dynamics.df, aes(x = CCF.change, y = mtDNA.change)) +
  geom_smooth(method = "lm") +
  geom_point(size = 1) +
  scale_x_log10("CCF change",
    breaks = c(0.01, 0.1, 1, 10, 100),
    labels = c("0.01", "0.1", "1", "10", "100")
  ) +
  scale_y_log10("mtDNA change",
    breaks = c(0.01, 0.1, 1, 10, 100),
    labels = c("0.01", "0.1", "1", "10", "100")
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title = element_text("Arial", size = 10, color = "black"),
    axis.text = element_text("Arial", size = 10, color = "black")
  )
ggsave("./figure_CLL_WES_RNA/figures/CRC/20231130_mtDNA_somatic_mutations_change.svg", width = 1.5, height = 1.5)
