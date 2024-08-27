library(dplyr)
library(ggplot2)
library(ggsignif)
library(scales)
library(tidyr)

meta <- readxl::read_excel("./data/CLL_RIC/20231229_CLL_RIC_mtDNA_meta.xlsx")

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

# gather all variants
### read all data
mgatk.data.ATAC <- 0
for (library in list.files("./data/bulk/CLL_RIC/", pattern = "*.mgatk")) {
  boo <- readRDS(paste0("./data/bulk/CLL_RIC/", library, "/final/mgatk.rds"))
  if (class(mgatk.data.ATAC) == "numeric") {
    mgatk.data.ATAC <- boo
  } else {
    mgatk.data.ATAC <- cbind(mgatk.data.ATAC, boo)
  }
}

### read variants
results.df <- data.frame()
for (library in list.files("./data/CLL_RIC/", pattern = "*.mgatk")) {
  message(library)
  results <- variant_calling_bulk(
    sample.bulk = readRDS(paste0("./data/CLL_RIC/", library, "/final/mgatk.rds")),
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

ref_all <- data.table::fread("./data/CRC/CRC-0001-441KN.mgatk/final/MT_refAllele.txt")


mutational.signature <- get_enrich_mutation_df(results.df$variant[which(results.df$heteroplasmy < 0.01)], ref_all)
mutational.signature$fc_called[which(mutational.signature$strand == "L")] <- -mutational.signature$fc_called[which(mutational.signature$strand == "L")]


ggplot(mutational.signature, aes(x = change_plot, fill = strand, color = strand, y = fc_called)) +
  # geom_boxplot(size=0.1, outlier.size = 0.2, outlier.colour = NULL, outlier.stroke = 0.1) +
  geom_col(color = "black") +
  scale_fill_manual(values = c("firebrick", "dodgerblue3")) +
  scale_color_manual(values = c("firebrick", "dodgerblue3")) +
  theme(legend.position = "bottom") +
  scale_y_continuous("Substitution Rate\n(Expected / Observed)", expand = c(0, 0), limits = c(-10, 10)) +
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
ggsave("./figure_CLL_WES_RNA/figures/CLL_RIC/20240104_mutation_signature_not_filtered.svg", width = 4, height = 2)


mutational.signature <- get_enrich_mutation_df(results.df$variant[which(results.df$heteroplasmy > 0.01)], ref_all)
mutational.signature$fc_called[which(mutational.signature$strand == "L")] <- -mutational.signature$fc_called[which(mutational.signature$strand == "L")]

ggplot(mutational.signature, aes(x = change_plot, fill = strand, color = strand, y = fc_called)) +
  # geom_boxplot(size=0.1, outlier.size = 0.2, outlier.colour = NULL, outlier.stroke = 0.1) +
  geom_col(color = "black") +
  scale_fill_manual(values = c("firebrick", "dodgerblue3")) +
  scale_color_manual(values = c("firebrick", "dodgerblue3")) +
  theme(legend.position = "bottom") +
  scale_y_continuous("Substitution Rate\n(Expected / Observed)", expand = c(0, 0), limits = c(-10, 10)) +
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
ggsave("./figure_CLL_WES_RNA/figures/CLL_RIC/20240104_mutation_signature_filtered.svg", width = 4, height = 2)
