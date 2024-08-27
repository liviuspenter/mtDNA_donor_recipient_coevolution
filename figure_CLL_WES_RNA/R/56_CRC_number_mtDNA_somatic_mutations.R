library(dplyr)
library(ggplot2)
library(ggsignif)
library(scales)
library(tidyr)

stable.patients <- c(15, 36, 14, 17, 18, 2, 40, 41)
natural.patients <- c(12, 3, 7, 1, 9, 19, 4)
dynamic.patients <- c(11, 20, 21, 39, 6, 10, 37, 13, 5, 16, 38)

category.colors <- c("stable" = "darkgreen", "natural.progression" = "magenta", "therapeutic.bottleneck" = "purple")

meta <- readxl::read_excel("./data/CRC/20231127_CRC_mtDNA_meta.xlsx")

# patients per category
meta$category <- "stable"
meta$category[which(meta$Patient %in% natural.patients)] <- "natural.progression"
meta$category[which(meta$Patient %in% dynamic.patients)] <- "therapeutic.bottleneck"
meta$category <- factor(meta$category, levels = c("stable", "natural.progression", "therapeutic.bottleneck"))

ggplot(data = meta %>% group_by(Patient) %>% filter(Day == max(Day)), aes(x = category, y = Day)) +
  geom_boxplot(outlier.size = 0, color = "black") +
  # stat_summary(geom='crossbar', fun = median, fun.min = median, fun.max = median, width=0.5, color='black') +
  scale_color_manual(values = category.colors) +
  scale_x_discrete(labels = c("stability", "natural\nprogression", "therapeutic\nbottleneck")) +
  scale_y_continuous("follow-up time") +
  geom_jitter(width = 0.1, size = 1, aes(color = category)) +
  geom_signif(
    comparisons = list(
      c("stable", "therapeutic.bottleneck"),
      c("natural.progression", "therapeutic.bottleneck")
    ),
    step_increase = 0.1, test = "t.test", color = "black", textsize = 3
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black"),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    axis.title.x = element_blank()
  )
ggsave("./figure_CLL_WES_RNA/figures/CRC/20240122_follow_up_cohort.svg", width = 1.5, height = 2)


mutations.df <- data.frame(
  patient = as.numeric(),
  mtDNA.mutations = as.numeric(),
  somatic.mutations = as.numeric()
)
variant.df.all <- data.frame()

for (p in unique(meta$Patient)) {
  variant.df <- read.table(paste0("./data/CRC/processed_together/", p, "_de_novo.csv"),
    check.names = F
  ) %>%
    as.data.frame()

  mutations.df <- rbind(mutations.df, data.frame(
    patient = p,
    mtDNA.mutations = length(unique(variant.df$variant)),
    mtDNA.mutations.stable = length(which(variant.df$significant == "non.significant")),
    mtDNA.mutations.dynamic = length(which(variant.df$significant == "significant")),
    somatic.mutations = NA
  ))

  variant.df$patient <- p
  variant.df.all <- rbind(variant.df.all, variant.df[, c(
    "variant", "pos", "REF", "ALT", "1.coverage", "1.REF.reads",
    "1.ALT.reads", "1.heteroplasmy", "p.value", "change", "direction",
    "significant", "patient"
  )])
}

# statistics
length(unique(variant.df.all$variant))
boo <- variant.df.all %>%
  group_by(patient) %>%
  tally() %>%
  filter(patient != 8)
summary(boo)

for (p in unique(meta$Patient)) {
  CCF.data <- data.table::fread(paste0("./data/CRC/MAF/edited/", p, ".mut_ccfs.txt"))
  CCF.data$identifier <- paste0(
    CCF.data$Chromosome, ":", CCF.data$Start_position, "_",
    CCF.data$Reference_Allele, ">", CCF.data$Tumor_Seq_Allele
  )
  mutations.df[which(mutations.df$patient == p), "somatic.mutations"] <- length(unique(CCF.data$identifier))

  baseline.sample <- meta$Sample[which(meta$Patient == p & meta$Day == 0)][1]
  followup.sample <- meta$Sample[which(meta$Patient == p & meta$Day == max(meta$Day[which(meta$Patient == p)]))]

  # count dynamic mutations
  dynamic.mutations <- 0
  for (mutation in unique(CCF.data$identifier)) {
    baseline.ref <- CCF.data$t_ref_count[which(CCF.data$identifier == mutation & CCF.data$Sample_ID == baseline.sample)]
    followup.ref <- CCF.data$t_ref_count[which(CCF.data$identifier == mutation & CCF.data$Sample_ID == followup.sample)]
    baseline.alt <- CCF.data$t_alt_count[which(CCF.data$identifier == mutation & CCF.data$Sample_ID == baseline.sample)]
    followup.alt <- CCF.data$t_alt_count[which(CCF.data$identifier == mutation & CCF.data$Sample_ID == followup.sample)]

    if (fisher.test(matrix(c(baseline.alt, baseline.ref, followup.alt, followup.ref), ncol = 2))$p.value < 0.05) {
      dynamic.mutations <- dynamic.mutations + 1
    }
  }
  mutations.df[which(mutations.df$patient == p), "dynamic.somatic.mutations"] <- dynamic.mutations

  # baseline somatic mutations with VAF >5%
  mutations.df[which(mutations.df$patient == p), "baseline.somatic.mutations"] <-
    length(unique(CCF.data$identifier[which(CCF.data$Sample_ID == baseline.sample &
      (CCF.data$t_alt_count / (CCF.data$t_alt_count + CCF.data$t_ref_count) > 0.05))]))

  mutations.df[which(mutations.df$patient == p), "follow.up"] <- max(meta$Day[which(meta$Patient == p)])
}

ggplot(mutations.df, aes(x = somatic.mutations, y = mtDNA.mutations)) +
  geom_point(size = 1) +
  geom_smooth(method = "lm") +
  scale_x_continuous("somatic nuclear mutations", limits = c(0, 200)) +
  scale_y_continuous("mtDNA mutations", limits = c(0, 30)) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title = element_text("Arial", size = 10, color = "black"),
    axis.text = element_text("Arial", size = 10, color = "black")
  )
ggsave("./figure_CLL_WES_RNA/figures/CRC/20231130_mtDNA_somatic_mutations_correlations.svg", width = 2, height = 2)

# correlation
# r=0.45
# p=0.02
cor.test(mutations.df$somatic.mutations, mutations.df$mtDNA.mutations)

# statistics of mutations per category
mutations.df$category <- "stable"
mutations.df$category[which(mutations.df$patient %in% natural.patients)] <- "natural.progression"
mutations.df$category[which(mutations.df$patient %in% dynamic.patients)] <- "therapeutic.bottleneck"
mutations.df$category <- factor(mutations.df$category, levels = c("stable", "natural.progression", "therapeutic.bottleneck"))

# dynamic mutations
ggplot(mutations.df, aes(x = category, y = mtDNA.mutations.dynamic, color = category)) +
  # stat_summary(geom='crossbar', fun = median, fun.min = median, fun.max = median, width=0.5, color='black') +
  geom_boxplot(outlier.size = 0, color = "black") +
  scale_color_manual(values = category.colors) +
  scale_x_discrete(labels = c("stability", "natural\nprogression", "therapeutic\nbottleneck")) +
  scale_y_continuous("dynamic mtDNA mutations") +
  geom_jitter(width = 0.1, size = 1) +
  geom_signif(
    comparisons = list(
      c("stable", "natural.progression"),
      c("stable", "therapeutic.bottleneck")
    ),
    step_increase = 0.1, test = "t.test", color = "black", textsize = 3
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black"),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    axis.title.x = element_blank()
  )
ggsave("./figure_CLL_WES_RNA/figures/CRC/20240119_dynamic_mtDNA_mutations.svg", width = 1.5, height = 2)

boo <- mutations.df %>%
  group_by(category, patient) %>%
  tally(mtDNA.mutations.dynamic)
summary(boo$n[which(boo$category == "stable")])

ggplot(mutations.df, aes(x = category, y = mtDNA.mutations, color = category)) +
  geom_boxplot(outlier.size = 0, color = "black") +
  # stat_summary(geom='crossbar', fun = median, fun.min = median, fun.max = median, width=0.5, color='black') +
  scale_color_manual(values = category.colors) +
  scale_x_discrete(labels = c("stability", "natural\nprogression", "therapeutic\nbottleneck")) +
  scale_y_continuous("mtDNA mutations") +
  geom_jitter(width = 0.1, size = 1) +
  geom_signif(
    comparisons = list(
      c("stable", "therapeutic.bottleneck"),
      c("natural.progression", "therapeutic.bottleneck")
    ),
    step_increase = 0.1, test = "t.test", color = "black", textsize = 3
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black"),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    axis.title.x = element_blank()
  )
ggsave("./figure_CLL_WES_RNA/figures/CRC/20240119_mtDNA_mutations.svg", width = 1.5, height = 2)

ggplot(mutations.df, aes(x = category, y = somatic.mutations, color = category)) +
  geom_boxplot(outlier.size = 0, color = "black") +
  # stat_summary(geom='crossbar', fun = median, fun.min = median, fun.max = median, width=0.5, color='black') +
  scale_color_manual(values = category.colors) +
  scale_x_discrete(labels = c("stability", "natural\nprogression", "therapeutic\nbottleneck")) +
  scale_y_continuous("somatic mutations") +
  geom_jitter(width = 0.1, size = 1) +
  geom_signif(
    comparisons = list(
      c("stable", "therapeutic.bottleneck"),
      c("natural.progression", "therapeutic.bottleneck")
    ),
    step_increase = 0.1, test = "t.test", color = "black", textsize = 3
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black"),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    axis.title.x = element_blank()
  )
ggsave("./figure_CLL_WES_RNA/figures/CRC/20240119_somatic_mutations.svg", width = 1.5, height = 2)

# dynamic somatic mutations
ggplot(mutations.df, aes(x = category, y = dynamic.somatic.mutations, color = category)) +
  # stat_summary(geom='crossbar', fun = median, fun.min = median, fun.max = median, width=0.5, color='black') +
  geom_boxplot(outlier.size = 0, color = "black") +
  scale_color_manual(values = category.colors) +
  scale_x_discrete(labels = c("stability", "natural\nprogression", "therapeutic\nbottleneck")) +
  scale_y_continuous("dynamic somatic mutations") +
  geom_jitter(width = 0.1, size = 1) +
  geom_signif(
    comparisons = list(
      c("stable", "natural.progression"),
      c("stable", "therapeutic.bottleneck"),
      c("natural.progression", "therapeutic.bottleneck")
    ),
    step_increase = 0.1, test = "t.test", color = "black", textsize = 3
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black"),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    axis.title.x = element_blank()
  )
ggsave("./figure_CLL_WES_RNA/figures/CRC/20240119_dynamic_somatic_mutations.svg", width = 1.5, height = 2)

# also significant for stable versus clonal evolution (natural progression + therapeutic bottleneck)
# p = 0.04
t.test(
  mutations.df$mtDNA.mutations[which(mutations.df$category == "stable")],
  mutations.df$mtDNA.mutations[which(mutations.df$category %in% c("natural.progression", "therapeutic.bottleneck"))]
)

ggplot(mutations.df, aes(x = follow.up, y = dynamic.somatic.mutations, color = category)) +
  geom_smooth(method = "lm", aes(group = 1), color = "black") +
  geom_point(size = 1) +
  scale_x_continuous("days follow-up", breaks = c(0, 2000, 4000), limits = c(0, 4000)) +
  scale_y_continuous("dynamic somatic mutations") +
  scale_color_manual(values = category.colors) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black")
  )
ggsave("./figure_CLL_WES_RNA/figures/CRC/20240122_follow_up_dynamic_somatic_mutations.svg", width = 1.5, height = 1.5)

ggplot(mutations.df, aes(x = follow.up, y = mtDNA.mutations.dynamic, color = category)) +
  geom_smooth(method = "lm", aes(group = 1), color = "black") +
  geom_point(size = 1) +
  scale_x_continuous("days follow-up", breaks = c(0, 2000, 4000), limits = c(0, 4000)) +
  scale_y_continuous("dynamic mtDNA mutations") +
  scale_color_manual(values = category.colors) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black")
  )
ggsave("./figure_CLL_WES_RNA/figures/CRC/20240122_follow_up_dynamic_mtDNA_mutations.svg", width = 1.5, height = 1.5)

ggplot(mutations.df, aes(x = dynamic.somatic.mutations, y = mtDNA.mutations.dynamic, color = category)) +
  geom_smooth(method = "lm", aes(group = 1), color = "black") +
  geom_point(size = 1) +
  scale_x_continuous("follow-up time") +
  scale_y_continuous("dynamic mtDNA mutations") +
  scale_color_manual(values = category.colors) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black")
  )
ggsave("./figure_CLL_WES_RNA/figures/CRC/20240122_follow_up_dynamic_mtDNA_mutations.svg", width = 1.5, height = 1.5)


# correlation followup and dynamic somatic mutations: r = 0.55; p = 0.03
cor.test(mutations.df$follow.up, mutations.df$dynamic.somatic.mutations)

# correlation followup and dynamic mtDNA mutations: r = 0.49; p = 0.01
cor.test(mutations.df$follow.up, mutations.df$mtDNA.mutations.dynamic)

# correlation dynamic mtDNA and somatic mutations: r = 0.45; p = 0.02
cor.test(mutations.df$dynamic.somatic.mutations, mutations.df$mtDNA.mutations.dynamic)







variant.df.all$identifier <- paste0(variant.df.all$patient, ".", variant.df.all$variant)
variant.df.all$category <- paste0(variant.df.all$direction, ".", variant.df.all$significant)
variant.df.all$category2 <- "stable"
variant.df.all$category2[which(variant.df.all$patient %in% natural.patients)] <- "natural.progression"
variant.df.all$category2[which(variant.df.all$patient %in% dynamic.patients)] <- "therapeutic.bottleneck"
variant.df.all$category2 <- factor(variant.df.all$category2, levels = c("stable", "natural.progression", "therapeutic.bottleneck"))

# number of baseline mutations with heteroplasmy > 0.5%
stats <- variant.df.all %>%
  group_by(category2, patient) %>%
  summarize(length(which(`1.heteroplasmy` > 0.005))) %>%
  as.data.frame()
colnames(stats) <- c("category2", "patient", "n")

ggplot(stats, aes(x = category2, y = n)) +
  geom_boxplot(outlier.size = 0) +
  geom_jitter(width = 0.1, size = 1, aes(color = category2)) +
  geom_signif(
    comparisons = list(
      c("stable", "therapeutic.bottleneck"),
      c("natural.progression", "therapeutic.bottleneck")
    ),
    step_increase = 0.2, test = "t.test", textsize = 3
  ) +
  scale_color_manual(values = category.colors) +
  scale_x_discrete(labels = c("stability", "natural\nprogression", "therapeutic\nbottleneck")) +
  scale_y_continuous("baseline mtDNA mutations", limits = c(0, 15)) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black"),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    axis.title.x = element_blank()
  )
ggsave("./figure_CLL_WES_RNA/figures/CRC/20240119_mtDNA_mutations_baseline.svg", width = 1.5, height = 2)
