library(dplyr)

source("./figure_mixing/bulk/R/00_variant_calling_atac.R")

### read all data
mgatk.data.ATAC <- 0
for (library in list.files("./data/CRC/", pattern = "*.mgatk")) {
  boo <- readRDS(paste0("./data/CRC/", library, "/final/mgatk.rds"))
  if (class(mgatk.data.ATAC) == "numeric") {
    mgatk.data.ATAC <- boo
  } else {
    mgatk.data.ATAC <- cbind(mgatk.data.ATAC, boo)
  }
}

# visualize coverage
boo <- cbind(
  data.frame(pos = seq(1, nrow(mgatk.data.ATAC))),
  as.data.frame(as.matrix(assays(mgatk.data.ATAC)[["coverage"]]))
)
boo <- pivot_longer(boo, cols = colnames(boo)[-1])
p <- ggplot(boo, aes(x = pos, y = value, group = name)) +
  scale_x_continuous("chrM") +
  scale_y_continuous("mean mtDNA coverage") +
  ggrastr::rasterize(geom_line(color = "blue"), dpi = 600) +
  coord_polar() +
  theme_bw() +
  theme(
    axis.title = element_text("Arial", size = 10, color = "black"),
    axis.text = element_text("Arial", size = 10, color = "black")
  )
ggsave("./figure_CLL_WES_RNA/figures/CRC/20231120_bulk_coverage.svg", width = 2, height = 2, plot = p)

# calculate coverage across

# Gene locations
GenePos.tib <- tibble(
  Names = c(
    "MT.ATP6", "MT.ATP8", "MT.COX1", "MT.COX2", "MT.COX3", "MT.CYTB", "MT.ND1", "MT.ND2", "MT.ND3",
    "MT.ND4", "MT.ND4L", "MT.ND5", "MT.ND6", "MT.RNR1", "MT.RNR2"
  ),
  start = c(8527, 8366, 5904, 7586, 9207, 14747, 3307, 4470, 10059, 10760, 10470, 12337, 14149, 648, 1671),
  end = c(9207, 8572, 7445, 8269, 9990, 15887, 4262, 5511, 10404, 12137, 10766, 14148, 14673, 1601, 3229)
)
GenePos.tib <- GenePos.tib %>%
  arrange(start) %>%
  mutate(mid = round((end - start) / 2 + start, 0), ycoord = rep(c(100, 80, 60), length.out = 15))
GenePos.tib$Names <- stringr::str_split_fixed(GenePos.tib$Names, pattern = "\\.", n = 2)[, 2]

cov.df <- boo %>%
  group_by(pos) %>%
  summarize(cov.mean = mean(value)) %>%
  as.data.frame()

gene.pos <-
  unlist(apply(GenePos.tib, MARGIN = 1, FUN = function(x) {
    seq(from = as.numeric(x["start"]), to = as.numeric(x["end"]))
  }))
cov.df$gene <- ifelse(cov.df$pos %in% gene.pos, "yes", "no")

summary(cov.df$cov.mean[which(cov.df$gene == "no")])
summary(cov.df$cov.mean[which(cov.df$gene == "yes")])
