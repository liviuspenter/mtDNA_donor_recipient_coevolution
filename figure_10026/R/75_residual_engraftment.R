library(ArchR)
library(dplyr)
library(gplots)
library(parallel)
library(Seurat)
library(ComplexHeatmap)

AML.colors <- c(
  "HSC" = RColorBrewer::brewer.pal(n = 5, name = "Reds")[5],
  "GMP" = RColorBrewer::brewer.pal(n = 5, name = "Reds")[4],
  "Mono" = RColorBrewer::brewer.pal(n = 5, name = "Reds")[3],
  "erythroid" = RColorBrewer::brewer.pal(n = 5, name = "Reds")[2]
)

AML.1010.mito <- loadArchRProject("./data/10026/AML.1010.mito/")
AML.1011.mito <- loadArchRProject("./data/10026/AML.1011.mito/")
AML.1012.mito <- loadArchRProject("./data/10026/AML.1012.mito/")
AML.1026.mito <- loadArchRProject("./data/10026/AML.1026.mito/")

AML.1010.chimerism <- data.table::fread("./data/10026/mtDNA/20200302_AML1010_donor_recipient.csv", header = T, dec = ",") %>% as.data.frame()
rownames(AML.1010.chimerism) <- AML.1010.chimerism$V1
AML.1010.chimerism[AML.1010.mito$cellNames, "celltype"] <- AML.1010.mito$manual.clusters
AML.1010.chimerism[AML.1010.mito$cellNames, "sample"] <- AML.1010.mito$Sample.hash

AML.1011.chimerism <- data.table::fread("./data/10026/mtDNA/20240522_AML1011_donor_recipient.csv", header = T, dec = ",") %>% as.data.frame()
AML.1011.chimerism$sample <- stringr::str_split_fixed(AML.1011.chimerism$V1, pattern = "#", n = 2)[, 1]

AML.1026.chimerism <- data.table::fread("./data/10026/mtDNA/20200302_AML1026_donor_recipient.csv", header = T, dec = ",") %>% as.data.frame()
rownames(AML.1026.chimerism) <- AML.1026.chimerism$V1
AML.1026.chimerism[AML.1026.mito$cellNames, "celltype"] <- AML.1026.mito$manual.clusters
AML.1026.chimerism[AML.1026.mito$cellNames, "sample"] <- AML.1026.mito$Sample.hash

df <- bind_rows(AML.1010.chimerism, AML.1011.chimerism, AML.1026.chimerism)
df$bc <- stringr::str_split_fixed(df$V1, pattern = "#", n = 2)[, 2]
df$patient <- stringr::str_split_fixed(df$sample, pattern = "_", n = 2)[, 1]

df <- df[which(df$sample != "none"), ]

df$celltype[which(df$celltype == "HSCT")] <- "HSC"
df$celltype[which(df$celltype == "Ery")] <- "erythroid"


stats <- df %>%
  filter(!sample %in% c("AML1011_A", "AML1011_B")) %>%
  filter(celltype %in% c("HSC", "GMP", "Mono", "erythroid")) %>%
  group_by(patient) %>%
  summarize(
    donor = length(which(individual == "donor")),
    recipient = length(which(individual == "recipient"))
  )
stats$residual.donor <- stats$donor / (stats$donor + stats$recipient)


stats <- df %>%
  filter(!sample %in% c("AML1011_A", "AML1011_B")) %>%
  filter(celltype %in% c("HSC", "GMP", "Mono", "erythroid")) %>%
  group_by(patient, celltype) %>%
  summarize(
    donor = length(which(individual == "donor")),
    recipient = length(which(individual == "recipient"))
  )
stats$residual.donor <- stats$donor / (stats$donor + stats$recipient)

ggplot(stats, aes(x = celltype, y = 100 * residual.donor, group = patient, fill = celltype)) +
  geom_col(position = "dodge", color = "black") +
  scale_x_discrete(limits = c("HSC", "GMP", "Mono", "erythroid")) +
  scale_y_continuous("% residual donor") +
  scale_fill_manual(values = AML.colors) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title = element_text("Arial", size = 10, color = "black"),
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )
ggsave("./figure_10026/figures/plots/20240604_residual_engraftment.svg", width = 1.3, height = 2)

## get T cell chimerism
Tcell.stats <- df %>%
  filter(!sample %in% c("AML1011_A", "AML1011_B")) %>%
  filter(celltype %in% c("CD4", "CD8", "T cell")) %>%
  group_by(patient, sample) %>%
  summarize(
    donor.cells = length(which(individual == "donor")),
    recipient.cells = length(which(individual == "recipient")),
    chimerism = length(which(individual == "donor")) / length(which(individual %in% c("donor", "recipient")))
  )

Tcell.stats$cells.sample <- Tcell.stats$donor.cells + Tcell.stats$recipient.cells

Tcell.stats$CI.1 <- apply(Tcell.stats, 1, function(x) prop.test(as.numeric(x["donor.cells"]), as.numeric(x["cells.sample"]))[["conf.int"]][1])
Tcell.stats$CI.2 <- apply(Tcell.stats, 1, function(x) prop.test(as.numeric(x["donor.cells"]), as.numeric(x["cells.sample"]))[["conf.int"]][2])

ggplot(Tcell.stats[which(Tcell.stats$patient == "AML1011"), ], aes(x = sample, y = 100 * chimerism)) +
  geom_point(size = 0.5) +
  geom_line(aes(group = 1)) +
  geom_errorbar(aes(ymin = 100 * CI.1, ymax = 100 * CI.2), width = 0.5) +
  scale_x_discrete(labels = c("Screening", "C2")) +
  scale_y_continuous("% donor T cell chimerism", limits = c(0, 100)) +
  theme_classic() +
  theme(
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black"),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )
ggsave("./figure_10026/figures/plots/20240604_chimerism_1011.svg", width = 1.2, height = 2.5)
