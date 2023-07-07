# helper functions for reading and merging data

# read vafs from tapestri file
read_vafs = function(file.name, sample.name) {
  vafs = as.data.frame(data.table::fread(file.name))
  vafs$Barcode = paste0(sample.name,'#',vafs$Barcode)
  rownames(vafs) = vafs$Barcode
  vafs = vafs[,-c(1,2)]
  
  vafs
}

# read vafs from mitochondrial vafs and merge with tapestri data
merge_with_mgatk = function(mission_bio_vafs, file.mgatk_vafs, sample.name, truncate_barcodes = F) {
  vafs.mgatk = as.data.frame(data.table::fread(file = file.mgatk_vafs))
  rownames(vafs.mgatk) = paste0(sample.name, '#',vafs.mgatk$V1)
  
  if (truncate_barcodes) {
    rownames(vafs.mgatk) = gsub(rownames(vafs.mgatk), pattern = '-1', replacement = '')
  }
  
  vafs.mgatk = 100*vafs.mgatk[,-1]
  vafs.mgatk$V1 = rownames(vafs.mgatk)
  mission_bio_vafs$V1 = rownames(mission_bio_vafs)
  merged = merge(mission_bio_vafs, vafs.mgatk, by = 'V1')
  rownames(merged) = merged$V1
  merged = merged[,-1]
  
  merged
}