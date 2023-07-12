library(stringr)

results.df.filtered = read.csv(file = './data/mixing/bulk/PRJNA486215_PRJNA563929.csv')

results.df.filtered$ALT = str_split_fixed(results.df.filtered$variant, pattern = '>', n=2)[,2]
results.df.filtered$variant.ALT = paste0(results.df.filtered$pos, results.df.filtered$ALT)

cat ('ID\tRange\tHaplogroup\tPolymorphisms\n', file = './data/mixing/bulk/PRJNA486215_PRJNA563929.hsd')
for (sample in unique(results.df.filtered$sample)) {
  polymorphisms = results.df.filtered$variant.ALT[which(results.df.filtered$sample == sample)]
  polymorphisms = paste(polymorphisms, collapse = '\t')
  cat (paste(c(gsub(sample, pattern = '\\.mgatk', replacement = ''), '1-16569', '?', polymorphisms), collapse = '\t'), sep = '\n',
       file = './data/mixing/bulk/PRJNA486215_PRJNA563929.hsd', append = T)
}


results.df.filtered = read.csv(file = './data/mixing/bulk/PRJNA741686.csv')

results.df.filtered$ALT = str_split_fixed(results.df.filtered$variant, pattern = '>', n=2)[,2]
results.df.filtered$variant.ALT = paste0(results.df.filtered$pos, results.df.filtered$ALT)

cat ('ID\tRange\tHaplogroup\tPolymorphisms\n', file = './data/mixing/bulk/PRJNA741686.hsd')
for (sample in unique(results.df.filtered$sample)) {
  polymorphisms = results.df.filtered$variant.ALT[which(results.df.filtered$sample == sample)]
  polymorphisms = paste(polymorphisms, collapse = '\t')
  cat (paste(c(gsub(sample, pattern = '\\.mgatk', replacement = ''), '1-16569', '?', polymorphisms), collapse = '\t'), sep = '\n',
       file = './data/mixing/bulk/PRJNA741686.hsd', append = T)
}