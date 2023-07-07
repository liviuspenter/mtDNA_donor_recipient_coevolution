import_kite_counts <- function(dir){
  mtx <- fread(paste0(dir, "/featurecounts.mtx"), header = FALSE)
  dim <- mtx[1,]
  mtx <- mtx[-1,]
  matx <- sparseMatrix(i = mtx[[1]], j = mtx[[2]], x = mtx[[3]])
  rownames(matx) <- fread(paste0(dir, "/featurecounts.barcodes.txt"), header = FALSE)[[1]]
  colnames(matx) <- paste0(fread(paste0(dir, "/featurecounts.genes.txt"), header = FALSE)[[1]])
  return(t(matx))
}

import_kite_counts2 <- function(dir){
  mtx <- fread(paste0(dir, "/featurecounts.mtx"), header = FALSE)
  dim <- mtx[1,]
  mtx <- mtx[-1,]
  matx <- sparseMatrix(i = mtx[[1]], j = mtx[[2]], x = mtx[[3]])
  rownames(matx) <- fread(paste0(dir, "/featurecounts.barcodes.txt"), header = FALSE)[[1]]
  colnames(matx) <- paste0(fread(paste0(dir, "/featurecounts.genes.txt"), header = FALSE)[[1]])[1:ncol(matx)]
  return(t(matx))
}

# https://stackoverflow.com/questions/43117608/r-binding-sparse-matrices-of-different-sizes-on-rows
merge.sparse <- function(...) {
  
  cnnew <- character()
  rnnew <- character()
  x <- vector()
  i <- numeric()
  j <- numeric()
  
  for (M in list(...)) {
    
    cnold <- colnames(M)
    rnold <- rownames(M)
    
    cnnew <- union(cnnew,cnold)
    rnnew <- union(rnnew,rnold)
    
    cindnew <- match(cnold,cnnew)
    rindnew <- match(rnold,rnnew)
    ind <- unname(which(M != 0,arr.ind=T))
    i <- c(i,rindnew[ind[,1]])
    j <- c(j,cindnew[ind[,2]])
    x <- c(x,M@x)
  }
  
  sparseMatrix(i=i,j=j,x=x,dims=c(length(rnnew),length(cnnew)),dimnames=list(rnnew,cnnew))
}