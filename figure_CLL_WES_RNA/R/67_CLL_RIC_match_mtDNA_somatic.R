library(cowplot)
library(dplyr)
library(ggplot2)
library(ggsignif)
library(scales)
library(tidyr)

meta <- readxl::read_excel("./data/CLL_RIC/20231229_CLL_RIC_mtDNA_meta.xlsx") %>%
  as.data.frame()
CCF.data <- readxl::read_excel("./data/CLL_RIC/MAF/20231229_CLL_RIC_CCF.xlsx") %>%
  as.data.frame()
meta2 <- readxl::read_excel("./data/CLL_RIC/20231229_CLL_RIC_mtDNA_meta.xlsx", sheet = 2) %>%
  as.data.frame()
rownames(meta2) <- meta2$Patient
meta$Day2 <- 0
meta[which(meta$Timepoint == 2), "Day2"] <- meta2[as.character(meta$Patient[which(meta$Timepoint == 2)]), "Days_between_samples"]


dynamics.df <- data.frame()
cluster.stats <- data.frame()
for (p in unique(meta$Patient)) {
  sample.first <- meta$Sample[which(meta$Patient == p)][1]
  sample.last <- last(meta$Sample[which(meta$Patient == p)])

  mtDNA.df <- read.table(paste0("./data/CLL_RIC//processed_together/", p, "_de_novo.csv"),
    check.names = F
  ) %>%
    as.data.frame()

  # remove variants that cannot be compared
  mtDNA.df <- mtDNA.df[which(!mtDNA.df$change %in% c(0, Inf)), ]

  cluster.data <- CCF.data[which(CCF.data$Patient == p), ]

  # calculate dynamics of clusters
  cluster.data$dynamics <- cluster.data$`Pre-HSCT CCF` / cluster.data$`Post-HSCT CCF`

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
  CCF.plot.data <- cluster.data %>%
    select("Cluster", "Pre-HSCT CCF", "Post-HSCT CCF") %>%
    pivot_longer(values_to = "ccf", cols = c("Pre-HSCT CCF", "Post-HSCT CCF"))
  CCF.plot.data$name <- ifelse(CCF.plot.data$name == "Pre-HSCT CCF", "pre-HSCT", "post-HSCT")

  CCF.plot.data$day <- sapply(CCF.plot.data$name, function(x) {
    meta$Day2[which(meta$Description == x & meta$Patient == p)]
  })

  timepoints <- colnames(mtDNA.df)[grepl("heteroplasmy", colnames(mtDNA.df))]
  mtDNA.plot.data <- mtDNA.df %>%
    select(c("variant", timepoints, cluster.match)) %>%
    pivot_longer(cols = timepoints)
  mtDNA.plot.data$day <- sapply(mtDNA.plot.data$name, function(x) {
    meta$Day2[which(paste0(meta$Timepoint, ".heteroplasmy") == x & meta$Patient == p)]
  })
  mtDNA.plot.data$value.log <- mtDNA.plot.data$value
  mtDNA.plot.data$value.log[which(mtDNA.plot.data$value == 0)] <- 0.0001

  color.scale <- as.character(BuenColors::jdb_palette(name = "corona"))
  names(color.scale) <- as.character(seq(1, length(color.scale)))

  CCF.plot <- ggplot() +
    geom_line(data = CCF.plot.data, aes(
      x = day, y = 100 * ccf,
      group = Cluster, color = as.character(Cluster)
    )) +
    geom_point(data = CCF.plot.data, aes(
      x = day, y = 100 * ccf,
      group = Cluster, color = as.character(Cluster)
    )) +
    scale_y_continuous("% CCF") +
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
    )) +
    scale_y_log10("% heteroplasmy",
      limits = c(0.01, 100), breaks = c(0.01, 0.1, 1, 10, 100),
      labels = c("ND", "0.1", "1", "10", "100")
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
  ggsave(paste0("./figure_CLL_WES_RNA/figures/CLL_RIC/integrated/20231229_", p, ".svg"),
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
# r = 0.92
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
ggsave("./figure_CLL_WES_RNA/figures/CLL_RIC/20231229_mtDNA_somatic_mutations_change.svg", width = 2, height = 2)
