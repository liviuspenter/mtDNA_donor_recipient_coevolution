library(dplyr)
library(ggplot2)
library(ggsignif)
library(scales)
library(tidyr)

source("./figure_mixing/bulk/R/00_variant_calling_atac.R")

### code from Caleb Lareau

# Simple reverse complement function
reverse_complement <- function(s) {
  chartr("ATGC", "TACG", s)
}

get_enrich_mutation_df <- function(called_variants, ref_all) {
  # Process 3 digit signature based on letters
  colnames(ref_all) <- c("pos", "ref")
  ref_all$ref <- toupper(ref_all$ref)
  l <- as.character(ref_all$ref)

  # Gs happen to be at the first and last position
  ref_all$three <- paste0(c("G", l[-length(l)]), l, c(l[-1], "G"))

  # Remove Ns
  ref_all <- ref_all[!grepl("N", ref_all$three), ]

  # Make every possible mutation
  ref_all_long <- rbind(ref_all, ref_all, ref_all, ref_all)
  ref_all_long$alt <- rep(c("A", "C", "G", "T"), each = dim(ref_all)[1])
  ref_all_long <- ref_all_long[ref_all_long$ref != ref_all_long$alt, ]

  # add some meta data
  ref_all_long$variant <- paste0(as.character(ref_all_long$pos), ref_all_long$ref, ">", ref_all_long$alt)
  ref_all_long$change <- paste0(ref_all_long$ref, ref_all_long$alt)
  ref_all_long$change_rc <- reverse_complement(paste0(ref_all_long$ref, ref_all_long$alt))

  # A/G rich strand is "heavy" -- https://en.wikipedia.org/wiki/Heavy_strand
  table(ref_all$ref) # so the reference strand is light (more C/T)
  ref_all_long$strand <- ifelse(ref_all_long$ref %in% c("C", "T"), "L", "H")

  # Change to C/T as ref allele
  ref_all_long$rc3 <- reverse_complement(ref_all_long$three)
  ref_all_long$three_plot <- ifelse(ref_all_long$strand == "L", ref_all_long$three, ref_all_long$rc3)
  ref_all_long$group_change <- ifelse(ref_all_long$strand == "L", ref_all_long$change, ref_all_long$change_rc)

  # Annotate with called variants
  ref_all_long$called <- ref_all_long$variant %in% called_variants

  # Compute changes in expected/observed
  total <- dim(ref_all_long)[1]
  total_called <- sum(ref_all_long$called)
  prop_df <- ref_all_long %>%
    group_by(three_plot, group_change, strand) %>%
    summarize(observed_prop_called = sum(called) / total_called, expected_prop = n() / total, n = n()) %>%
    mutate(fc_called = observed_prop_called / expected_prop)

  prop_df$change_plot <- paste0(prop_df$group_change, "_", prop_df$three_plot)
  prop_df
}
###

ref_all <- data.table::fread("./data/10026/AML1010_1.mgatk/chrM_refAllele.txt")

AML1007.variants <- data.table::fread("./data/10026/mtDNA/20240516_AML1007_high_confidence_variants.csv")
AML1010.variants <- data.table::fread("./data/10026/mtDNA/20210429_AML1010_high_confidence_variants.csv")
AML1011.variants <- data.table::fread("./data/10026/mtDNA/20240516_AML1011_high_confidence_variants.csv")
AML1012.variants <- data.table::fread("./data/10026/mtDNA/20210429_AML1012_high_confidence_variants.csv")
AML1026.variants <- data.table::fread("./data/10026/mtDNA/20210429_AML1026_high_confidence_variants.csv")
IST1.variants <- data.table::fread("./data/IST/mtDNA/20220110_IST1_2_high_confidence_variants.csv")
IST2.variants <- data.table::fread("./data/IST/mtDNA/20220110_IST3_high_confidence_variants.csv")
IST3.variants <- data.table::fread("./data/IST/mtDNA/20220110_IST4_high_confidence_variants.csv")
IST4.variants <- data.table::fread("./data/IST/mtDNA/20220110_IST5_high_confidence_variants.csv")
MART1.variants <- data.table::fread("./data/mixing/MART1/mtDNA/20211030_MART1_high_confidence_variants.csv")


mutational.signature <- get_enrich_mutation_df(AML1007.variants$variant, ref_all)
mutational.signature$fc_called[which(mutational.signature$strand == "L")] <- -mutational.signature$fc_called[which(mutational.signature$strand == "L")]
ggplot(mutational.signature, aes(x = change_plot, fill = strand, color = strand, y = fc_called)) +
  # geom_boxplot(size=0.1, outlier.size = 0.2, outlier.colour = NULL, outlier.stroke = 0.1) +
  geom_col(color = "black") +
  scale_fill_manual(values = c("firebrick", "dodgerblue3")) +
  scale_color_manual(values = c("firebrick", "dodgerblue3")) +
  theme(legend.position = "bottom") +
  scale_y_continuous("Substitution Rate\n(Expected / Observed)", expand = c(0, 0), limits = c(-35, 35)) +
  geom_hline(yintercept = 0, linetype = 2, color = "black", size = 0.25) +
  geom_vline(xintercept = c(16, 32, 48, 64, 80), linetype = 2, size = 0.25) +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(colour = "black", size = 0.25),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title = element_text("Arial", size = 10, colour = "black"),
    axis.text = element_text("Arial", size = 10, colour = "black"),
    axis.ticks = element_line(size = 0.25, colour = "black")
  )
ggsave("./figure_10026/figures/qc/20240625_AML1007_mutational_profile.svg", width = 3.5, height = 1.5)


mutational.signature <- get_enrich_mutation_df(AML1010.variants$variant, ref_all)
mutational.signature$fc_called[which(mutational.signature$strand == "L")] <- -mutational.signature$fc_called[which(mutational.signature$strand == "L")]
ggplot(mutational.signature, aes(x = change_plot, fill = strand, color = strand, y = fc_called)) +
  # geom_boxplot(size=0.1, outlier.size = 0.2, outlier.colour = NULL, outlier.stroke = 0.1) +
  geom_col(color = "black") +
  scale_fill_manual(values = c("firebrick", "dodgerblue3")) +
  scale_color_manual(values = c("firebrick", "dodgerblue3")) +
  theme(legend.position = "bottom") +
  scale_y_continuous("Substitution Rate\n(Expected / Observed)", expand = c(0, 0), limits = c(-30, 30)) +
  geom_hline(yintercept = 0, linetype = 2, color = "black", size = 0.25) +
  geom_vline(xintercept = c(16, 32, 48, 64, 80), linetype = 2, size = 0.25) +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(colour = "black", size = 0.25),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title = element_text("Arial", size = 10, colour = "black"),
    axis.text = element_text("Arial", size = 10, colour = "black"),
    axis.ticks = element_line(size = 0.25, colour = "black")
  )
ggsave("./figure_10026/figures/qc/20240625_AML1010_mutational_profile.svg", width = 3.5, height = 1.5)


mutational.signature <- get_enrich_mutation_df(AML1011.variants$variant, ref_all)
mutational.signature$fc_called[which(mutational.signature$strand == "L")] <- -mutational.signature$fc_called[which(mutational.signature$strand == "L")]
ggplot(mutational.signature, aes(x = change_plot, fill = strand, color = strand, y = fc_called)) +
  # geom_boxplot(size=0.1, outlier.size = 0.2, outlier.colour = NULL, outlier.stroke = 0.1) +
  geom_col(color = "black") +
  scale_fill_manual(values = c("firebrick", "dodgerblue3")) +
  scale_color_manual(values = c("firebrick", "dodgerblue3")) +
  theme(legend.position = "bottom") +
  scale_y_continuous("Substitution Rate\n(Expected / Observed)", expand = c(0, 0), limits = c(-25, 25)) +
  geom_hline(yintercept = 0, linetype = 2, color = "black", size = 0.25) +
  geom_vline(xintercept = c(16, 32, 48, 64, 80), linetype = 2, size = 0.25) +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(colour = "black", size = 0.25),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title = element_text("Arial", size = 10, colour = "black"),
    axis.text = element_text("Arial", size = 10, colour = "black"),
    axis.ticks = element_line(size = 0.25, colour = "black")
  )
ggsave("./figure_10026/figures/qc/20240625_AML1011_mutational_profile.svg", width = 3.5, height = 1.5)


mutational.signature <- get_enrich_mutation_df(AML1012.variants$variant, ref_all)
mutational.signature$fc_called[which(mutational.signature$strand == "L")] <- -mutational.signature$fc_called[which(mutational.signature$strand == "L")]
ggplot(mutational.signature, aes(x = change_plot, fill = strand, color = strand, y = fc_called)) +
  # geom_boxplot(size=0.1, outlier.size = 0.2, outlier.colour = NULL, outlier.stroke = 0.1) +
  geom_col(color = "black") +
  scale_fill_manual(values = c("firebrick", "dodgerblue3")) +
  scale_color_manual(values = c("firebrick", "dodgerblue3")) +
  theme(legend.position = "bottom") +
  scale_y_continuous("Substitution Rate\n(Expected / Observed)", expand = c(0, 0), limits = c(-30, 30)) +
  geom_hline(yintercept = 0, linetype = 2, color = "black", size = 0.25) +
  geom_vline(xintercept = c(16, 32, 48, 64, 80), linetype = 2, size = 0.25) +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(colour = "black", size = 0.25),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title = element_text("Arial", size = 10, colour = "black"),
    axis.text = element_text("Arial", size = 10, colour = "black"),
    axis.ticks = element_line(size = 0.25, colour = "black")
  )
ggsave("./figure_10026/figures/qc/20240625_AML1012_mutational_profile.svg", width = 3.5, height = 1.5)


mutational.signature <- get_enrich_mutation_df(AML1026.variants$variant, ref_all)
mutational.signature$fc_called[which(mutational.signature$strand == "L")] <- -mutational.signature$fc_called[which(mutational.signature$strand == "L")]
ggplot(mutational.signature, aes(x = change_plot, fill = strand, color = strand, y = fc_called)) +
  # geom_boxplot(size=0.1, outlier.size = 0.2, outlier.colour = NULL, outlier.stroke = 0.1) +
  geom_col(color = "black") +
  scale_fill_manual(values = c("firebrick", "dodgerblue3")) +
  scale_color_manual(values = c("firebrick", "dodgerblue3")) +
  theme(legend.position = "bottom") +
  scale_y_continuous("Substitution Rate\n(Expected / Observed)", expand = c(0, 0), limits = c(-35, 35)) +
  geom_hline(yintercept = 0, linetype = 2, color = "black", size = 0.25) +
  geom_vline(xintercept = c(16, 32, 48, 64, 80), linetype = 2, size = 0.25) +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(colour = "black", size = 0.25),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title = element_text("Arial", size = 10, colour = "black"),
    axis.text = element_text("Arial", size = 10, colour = "black"),
    axis.ticks = element_line(size = 0.25, colour = "black")
  )
ggsave("./figure_10026/figures/qc/20240625_AML1026_mutational_profile.svg", width = 3.5, height = 1.5)




mutational.signature <- get_enrich_mutation_df(MART1.variants$variant, ref_all)
mutational.signature$fc_called[which(mutational.signature$strand == "L")] <- -mutational.signature$fc_called[which(mutational.signature$strand == "L")]
ggplot(mutational.signature, aes(x = change_plot, fill = strand, color = strand, y = fc_called)) +
  # geom_boxplot(size=0.1, outlier.size = 0.2, outlier.colour = NULL, outlier.stroke = 0.1) +
  geom_col(color = "black") +
  scale_fill_manual(values = c("firebrick", "dodgerblue3")) +
  scale_color_manual(values = c("firebrick", "dodgerblue3")) +
  theme(legend.position = "bottom") +
  scale_y_continuous("Substitution Rate\n(Expected / Observed)", expand = c(0, 0), limits = c(-25, 25)) +
  geom_hline(yintercept = 0, linetype = 2, color = "black", size = 0.25) +
  geom_vline(xintercept = c(16, 32, 48, 64, 80), linetype = 2, size = 0.25) +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(colour = "black", size = 0.25),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title = element_text("Arial", size = 10, colour = "black"),
    axis.text = element_text("Arial", size = 10, colour = "black"),
    axis.ticks = element_line(size = 0.25, colour = "black")
  )
ggsave("./figure_10026/figures/qc/20240625_MART1_mutational_profile.svg", width = 3.5, height = 1.5)



mutational.signature <- get_enrich_mutation_df(IST1.variants$variant, ref_all)
mutational.signature$fc_called[which(mutational.signature$strand == "L")] <- -mutational.signature$fc_called[which(mutational.signature$strand == "L")]
ggplot(mutational.signature, aes(x = change_plot, fill = strand, color = strand, y = fc_called)) +
  # geom_boxplot(size=0.1, outlier.size = 0.2, outlier.colour = NULL, outlier.stroke = 0.1) +
  geom_col(color = "black") +
  scale_fill_manual(values = c("firebrick", "dodgerblue3")) +
  scale_color_manual(values = c("firebrick", "dodgerblue3")) +
  theme(legend.position = "bottom") +
  scale_y_continuous("Substitution Rate\n(Expected / Observed)", expand = c(0, 0), limits = c(-25, 25)) +
  geom_hline(yintercept = 0, linetype = 2, color = "black", size = 0.25) +
  geom_vline(xintercept = c(16, 32, 48, 64, 80), linetype = 2, size = 0.25) +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(colour = "black", size = 0.25),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title = element_text("Arial", size = 10, colour = "black"),
    axis.text = element_text("Arial", size = 10, colour = "black"),
    axis.ticks = element_line(size = 0.25, colour = "black")
  )
ggsave("./figure_10026/figures/qc/20240625_IST1_mutational_profile.svg", width = 3.5, height = 1.5)

mutational.signature <- get_enrich_mutation_df(IST2.variants$variant, ref_all)
mutational.signature$fc_called[which(mutational.signature$strand == "L")] <- -mutational.signature$fc_called[which(mutational.signature$strand == "L")]
ggplot(mutational.signature, aes(x = change_plot, fill = strand, color = strand, y = fc_called)) +
  # geom_boxplot(size=0.1, outlier.size = 0.2, outlier.colour = NULL, outlier.stroke = 0.1) +
  geom_col(color = "black") +
  scale_fill_manual(values = c("firebrick", "dodgerblue3")) +
  scale_color_manual(values = c("firebrick", "dodgerblue3")) +
  theme(legend.position = "bottom") +
  scale_y_continuous("Substitution Rate\n(Expected / Observed)", expand = c(0, 0), limits = c(-30, 30)) +
  geom_hline(yintercept = 0, linetype = 2, color = "black", size = 0.25) +
  geom_vline(xintercept = c(16, 32, 48, 64, 80), linetype = 2, size = 0.25) +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(colour = "black", size = 0.25),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title = element_text("Arial", size = 10, colour = "black"),
    axis.text = element_text("Arial", size = 10, colour = "black"),
    axis.ticks = element_line(size = 0.25, colour = "black")
  )
ggsave("./figure_10026/figures/qc/20240625_IST2_mutational_profile.svg", width = 3.5, height = 1.5)


mutational.signature <- get_enrich_mutation_df(IST3.variants$variant, ref_all)
mutational.signature$fc_called[which(mutational.signature$strand == "L")] <- -mutational.signature$fc_called[which(mutational.signature$strand == "L")]
ggplot(mutational.signature, aes(x = change_plot, fill = strand, color = strand, y = fc_called)) +
  # geom_boxplot(size=0.1, outlier.size = 0.2, outlier.colour = NULL, outlier.stroke = 0.1) +
  geom_col(color = "black") +
  scale_fill_manual(values = c("firebrick", "dodgerblue3")) +
  scale_color_manual(values = c("firebrick", "dodgerblue3")) +
  theme(legend.position = "bottom") +
  scale_y_continuous("Substitution Rate\n(Expected / Observed)", expand = c(0, 0), limits = c(-30, 30)) +
  geom_hline(yintercept = 0, linetype = 2, color = "black", size = 0.25) +
  geom_vline(xintercept = c(16, 32, 48, 64, 80), linetype = 2, size = 0.25) +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(colour = "black", size = 0.25),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title = element_text("Arial", size = 10, colour = "black"),
    axis.text = element_text("Arial", size = 10, colour = "black"),
    axis.ticks = element_line(size = 0.25, colour = "black")
  )
ggsave("./figure_10026/figures/qc/20240625_IST3_mutational_profile.svg", width = 3.5, height = 1.5)


mutational.signature <- get_enrich_mutation_df(IST4.variants$variant, ref_all)
mutational.signature$fc_called[which(mutational.signature$strand == "L")] <- -mutational.signature$fc_called[which(mutational.signature$strand == "L")]
ggplot(mutational.signature, aes(x = change_plot, fill = strand, color = strand, y = fc_called)) +
  # geom_boxplot(size=0.1, outlier.size = 0.2, outlier.colour = NULL, outlier.stroke = 0.1) +
  geom_col(color = "black") +
  scale_fill_manual(values = c("firebrick", "dodgerblue3")) +
  scale_color_manual(values = c("firebrick", "dodgerblue3")) +
  theme(legend.position = "bottom") +
  scale_y_continuous("Substitution Rate\n(Expected / Observed)", expand = c(0, 0), limits = c(-40, 40)) +
  geom_hline(yintercept = 0, linetype = 2, color = "black", size = 0.25) +
  geom_vline(xintercept = c(16, 32, 48, 64, 80), linetype = 2, size = 0.25) +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(colour = "black", size = 0.25),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title = element_text("Arial", size = 10, colour = "black"),
    axis.text = element_text("Arial", size = 10, colour = "black"),
    axis.ticks = element_line(size = 0.25, colour = "black")
  )
ggsave("./figure_10026/figures/qc/20240625_IST4_mutational_profile.svg", width = 3.5, height = 1.5)




## look into AML1011
library(ArchR)
library(Signac)
AML1011.mito <- loadArchRProject("./data/10026//AML.1011.mito/")
AML1011.A.mgatk <- ReadMGATK("./data/10026/Pool148_13.mgatk/final/")
AML1011.B.mgatk <- ReadMGATK("./data/10026/Pool148_14.mgatk/final/")
AML1011.1.mgatk <- ReadMGATK("./data/10026/Pool148_15.mgatk/final/")
AML1011.2.mgatk <- ReadMGATK("./data/10026/Pool148_16.mgatk/final/")

colnames(AML1011.A.mgatk$counts) <- paste0("AML1011_A#", colnames(AML1011.A.mgatk$counts))
rownames(AML1011.A.mgatk$depth) <- paste0("AML1011_A#", rownames(AML1011.A.mgatk$depth))
colnames(AML1011.B.mgatk$counts) <- paste0("AML1011_B#", colnames(AML1011.B.mgatk$counts))
rownames(AML1011.B.mgatk$depth) <- paste0("AML1011_B#", rownames(AML1011.B.mgatk$depth))

colnames(AML1011.1.mgatk$counts) <- paste0("AML1011_1#", colnames(AML1011.1.mgatk$counts))
rownames(AML1011.1.mgatk$depth) <- paste0("AML1011_1#", rownames(AML1011.1.mgatk$depth))
colnames(AML1011.2.mgatk$counts) <- paste0("AML1011_2#", colnames(AML1011.2.mgatk$counts))
rownames(AML1011.2.mgatk$depth) <- paste0("AML1011_2#", rownames(AML1011.2.mgatk$depth))

AML1011.pre.mgatk <- AML1011.A.mgatk
AML1011.pre.mgatk$counts <- cbind(AML1011.A.mgatk$counts, AML1011.B.mgatk$counts)
AML1011.pre.mgatk$depth <- rbind(AML1011.A.mgatk$depth, AML1011.B.mgatk$depth)
AML1011.pre.mgatk$counts <- AML1011.pre.mgatk$counts[, intersect(colnames(AML1011.pre.mgatk$counts), AML1011.mito$cellNames)]
AML1011.pre.mgatk$depth <- subset(AML1011.pre.mgatk$depth, rownames(AML1011.pre.mgatk$depth) %in% colnames(AML1011.pre.mgatk$counts))
AML1011.pre.mutations <- IdentifyVariants(object = AML1011.pre.mgatk$counts, refallele = AML1011.pre.mgatk$refallele)
AML1011.pre.mutations.high.conf <- subset(AML1011.pre.mutations, subset = n_cells_conf_detected >= 5 & strand_correlation >= 0.65 & vmr > 0.01)
write.table(AML1011.pre.mutations.high.conf, file = "./data/10026/mtDNA/20240625_AML1011_pre_high_confidence_variants.csv")

AML1011.post.mgatk <- AML1011.1.mgatk
AML1011.post.mgatk$counts <- cbind(AML1011.1.mgatk$counts, AML1011.2.mgatk$counts)
AML1011.post.mgatk$depth <- rbind(AML1011.1.mgatk$depth, AML1011.2.mgatk$depth)
AML1011.post.mgatk$counts <- AML1011.post.mgatk$counts[, intersect(colnames(AML1011.post.mgatk$counts), AML1011.mito$cellNames)]
AML1011.post.mgatk$depth <- subset(AML1011.post.mgatk$depth, rownames(AML1011.post.mgatk$depth) %in% colnames(AML1011.post.mgatk$counts))
AML1011.post.mutations <- IdentifyVariants(object = AML1011.post.mgatk$counts, refallele = AML1011.post.mgatk$refallele)
AML1011.post.mutations.high.conf <- subset(AML1011.post.mutations, subset = n_cells_conf_detected >= 5 & strand_correlation >= 0.65 & vmr > 0.01)
write.table(AML1011.post.mutations.high.conf, file = "./data/10026/mtDNA/20240625_AML1011_post_high_confidence_variants.csv")

AML1011.pre.mutations.high.conf <- data.table::fread("./data/10026/mtDNA/20240625_AML1011_pre_high_confidence_variants.csv")
AML1011.post.mutations.high.conf <- data.table::fread("./data/10026/mtDNA/20240625_AML1011_post_high_confidence_variants.csv")

mutational.signature <- get_enrich_mutation_df(AML1011.pre.mutations.high.conf$variant, ref_all)
mutational.signature$fc_called[which(mutational.signature$strand == "L")] <- -mutational.signature$fc_called[which(mutational.signature$strand == "L")]
ggplot(mutational.signature, aes(x = change_plot, fill = strand, color = strand, y = fc_called)) +
  # geom_boxplot(size=0.1, outlier.size = 0.2, outlier.colour = NULL, outlier.stroke = 0.1) +
  geom_col(color = "black") +
  scale_fill_manual(values = c("firebrick", "dodgerblue3")) +
  scale_color_manual(values = c("firebrick", "dodgerblue3")) +
  theme(legend.position = "bottom") +
  scale_y_continuous("Substitution Rate\n(Expected / Observed)", expand = c(0, 0), limits = c(-45, 45)) +
  geom_hline(yintercept = 0, linetype = 2, color = "black", size = 0.25) +
  geom_vline(xintercept = c(16, 32, 48, 64, 80), linetype = 2, size = 0.25) +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(colour = "black", size = 0.25),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title = element_text("Arial", size = 10, colour = "black"),
    axis.text = element_text("Arial", size = 10, colour = "black"),
    axis.ticks = element_line(size = 0.25, colour = "black")
  )
ggsave("./figure_10026/figures/qc/20240625_AML1011_pre_mutational_profile.svg", width = 3.5, height = 1.5)

mutational.signature <- get_enrich_mutation_df(AML1011.post.mutations.high.conf$variant, ref_all)
mutational.signature$fc_called[which(mutational.signature$strand == "L")] <- -mutational.signature$fc_called[which(mutational.signature$strand == "L")]
ggplot(mutational.signature, aes(x = change_plot, fill = strand, color = strand, y = fc_called)) +
  # geom_boxplot(size=0.1, outlier.size = 0.2, outlier.colour = NULL, outlier.stroke = 0.1) +
  geom_col(color = "black") +
  scale_fill_manual(values = c("firebrick", "dodgerblue3")) +
  scale_color_manual(values = c("firebrick", "dodgerblue3")) +
  theme(legend.position = "bottom") +
  scale_y_continuous("Substitution Rate\n(Expected / Observed)", expand = c(0, 0), limits = c(-30, 30)) +
  geom_hline(yintercept = 0, linetype = 2, color = "black", size = 0.25) +
  geom_vline(xintercept = c(16, 32, 48, 64, 80), linetype = 2, size = 0.25) +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(colour = "black", size = 0.25),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title = element_text("Arial", size = 10, colour = "black"),
    axis.text = element_text("Arial", size = 10, colour = "black"),
    axis.ticks = element_line(size = 0.25, colour = "black")
  )
ggsave("./figure_10026/figures/qc/20240625_AML1011_post_mutational_profile.svg", width = 3.5, height = 1.5)




## look into AML1012
library(ArchR)
library(Signac)
AML1012.mito <- loadArchRProject("./data/10026/AML.1012.all.mito//")
AML1012.A.mgatk <- ReadMGATK("./data/10026/Pool148_17.mgatk/final/")

colnames(AML1012.A.mgatk$counts) <- paste0("AML1012_A#", colnames(AML1012.A.mgatk$counts))
rownames(AML1012.A.mgatk$depth) <- paste0("AML1012_A#", rownames(AML1012.A.mgatk$depth))

AML1012.pre.mgatk <- AML1012.A.mgatk
AML1012.pre.mgatk$counts <- AML1012.pre.mgatk$counts[, intersect(colnames(AML1012.pre.mgatk$counts), AML1012.mito$cellNames)]
AML1012.pre.mgatk$depth <- subset(AML1012.pre.mgatk$depth, rownames(AML1012.pre.mgatk$depth) %in% colnames(AML1012.pre.mgatk$counts))
AML1012.pre.mutations <- IdentifyVariants(object = AML1012.pre.mgatk$counts, refallele = AML1012.pre.mgatk$refallele)
AML1012.pre.mutations.high.conf <- subset(AML1012.pre.mutations, subset = n_cells_conf_detected >= 5 & strand_correlation >= 0.65 & vmr > 0.01)
write.table(AML1012.pre.mutations.high.conf, file = "./data/10026/mtDNA/20240625_AML1012_pre_high_confidence_variants.csv")

AML1012.pre.mutations.high.conf <- data.table::fread("./data/10026/mtDNA/20240625_AML1012_pre_high_confidence_variants.csv")

mutational.signature <- get_enrich_mutation_df(AML1012.pre.mutations.high.conf$variant, ref_all)
mutational.signature$fc_called[which(mutational.signature$strand == "L")] <- -mutational.signature$fc_called[which(mutational.signature$strand == "L")]
ggplot(mutational.signature, aes(x = change_plot, fill = strand, color = strand, y = fc_called)) +
  # geom_boxplot(size=0.1, outlier.size = 0.2, outlier.colour = NULL, outlier.stroke = 0.1) +
  geom_col(color = "black") +
  scale_fill_manual(values = c("firebrick", "dodgerblue3")) +
  scale_color_manual(values = c("firebrick", "dodgerblue3")) +
  theme(legend.position = "bottom") +
  scale_y_continuous("Substitution Rate\n(Expected / Observed)", expand = c(0, 0), limits = c(-60, 60)) +
  geom_hline(yintercept = 0, linetype = 2, color = "black", size = 0.25) +
  geom_vline(xintercept = c(16, 32, 48, 64, 80), linetype = 2, size = 0.25) +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(colour = "black", size = 0.25),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title = element_text("Arial", size = 10, colour = "black"),
    axis.text = element_text("Arial", size = 10, colour = "black"),
    axis.ticks = element_line(size = 0.25, colour = "black")
  )
ggsave("./figure_10026/figures/qc/20240625_AML1012_pre_mutational_profile.svg", width = 3.5, height = 1.5)
