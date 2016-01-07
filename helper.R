#utils-shrayes.R



###COLOR FUNCTIONS###

#takes a vector of a continuous variable and splits it into ncolor colors in a gradient from start to end

convertcolor <- function(v, ncolors = 20, percentile=FALSE, start="lightgray",end="darkred") {
  if (!require("RColorBrewer")) {
    install.packages("RColorBrewer")
    library(RColorBrewer)
  }
  scaled <- ifelse( rep(percentile,length(v)),
                    rank(v)/length(v),
                    (v - min(v,na.rm=TRUE)) / (max(v,na.rm=TRUE)-min(v,na.rm=TRUE))
  )
  colfunc <- colorRampPalette(c(start, end))
  binned <- cut(scaled, seq(0,1, by = 1/(ncolors)), include.lowest=TRUE, labels = colfunc(ncolors) )
  return (as.character(binned))
}

convertcolor_pal <- function(v, ncolors = 10, percentile=FALSE, palette = "BrBG") {
  if (!require("RColorBrewer")) {
    install.packages("RColorBrewer")
    library(RColorBrewer)
  }
  scaled <- ifelse( rep(percentile,length(v)),
                    rank(v)/length(v),
                    (v - min(v)) / (max(v)-min(v))
  )
  binned <- cut(scaled, seq(0,1, by = 1/(ncolors)), include.lowest=TRUE, labels = brewer.pal(ncolors,palette) )
  return (as.character(binned))
}

###2D Histogram###

#2d manual histogram
hist2d.do <- function (x, y, nbins, xlab="", ylab="", main="", x.bin= NULL, y.bin=NULL) {
  if(is.null(x.bin)) {
    x.bin <- seq(floor(min(x)) , ceiling(max(x)), length = nbins)
  }
  if(is.null(y.bin)) {
    y.bin <- seq(floor(min(y)) , ceiling(max(y)), length = nbins)
  }
  
  xcut <- cut(x,breaks=x.bin,include.lowest=TRUE)
  ycut <- cut(y,breaks=y.bin,include.lowest=TRUE)
  #xcut <- cut(x1,nbins)
  #ycut <- cut(x2,nbins)
  bin2d<- as.matrix(table(xcut,ycut))
  
  #image(x.bin,y.bin,log10(bin2d+1), col=c("black",topo.colors(max(bin2d))), xlab=xlab, ylab=ylab, main=main)
  #contour(x.bin[-length(x.bin)], y.bin[-length(y.bin)],log10(bin2d), add=TRUE, lwd=1, col=rgb(1,1,1,.7))
  
  image(x.bin,y.bin,bin2d, col=c("black",topo.colors(max(bin2d))), xlab=xlab, ylab=ylab, main=main)
  contour(x.bin[-length(x.bin)], y.bin[-length(y.bin)],bin2d, add=TRUE, lwd=1, col=rgb(1,1,1,.7))
  
  list(x.bin = x.bin,
       y.bin = y.bin,
       xcut = xcut,
       ycut = ycut,
       bin2d = bin2d,
       contourlines = contourLines(x.bin[-length(x.bin)], y.bin[-length(y.bin)],log(1+bin2d)))
}




#function to construct cosine-similarity adjacency from dtm:
adj.do <- function(dtm) {
  
  #matrix math is easier with sparseMatrix class in matrix package
  datamatrix<- sparseMatrix(i=dtm$i,j=dtm$j,x=dtm$v)#, dimnames=dimnames(dtm))
  
  #tf weighted
  wordcount <- col_sums(dtm) #note, col_sums requires simple_triplet_matrix version dtm
  normeddtm <- t(t(datamatrix) /wordcount)
  crossp <- tcrossprod(normeddtm)
  
  norm.event <- sqrt(diag(crossp))
  adj.all<- as.matrix(crossp)/(norm.event %o% norm.event)
  
  adj.all
}


####COMPUTE NETWORK LAYOUTS####
layout.do <- function(adj.all, layout.start = NULL){
  forumnet <- graph.adjacency(adj.all, mode="undirected", weighted=TRUE, diag=FALSE)
  print(system.time({
    l4p <- layout.fruchterman.reingold(forumnet, 
                                       area=vcount(forumnet)^2,
                                       repulserad=vcount(forumnet)^3,
                                       weights = E(forumnet)$weight,
                                       start = layout.start)
  }))
  l4p
}


###closest point to selected###
returnpoint<- function(l4p) {
  click<- locator(1)
  dist<- abs(l4p[,1]-click$x) + abs(l4p[,2]-click$y)
  me<- which.min(dist)
  me
}


###WORKING WITH MATRIX,SIMPLE TRIPLET MATRICES, ETC.

as.sparseMatrix <- function(simple_triplet_matrix_sparse) {
  retval <-  sparseMatrix(i=as.numeric(simple_triplet_matrix_sparse$i),
                          j=as.numeric(simple_triplet_matrix_sparse$j),
                          x=as.numeric(as.character(simple_triplet_matrix_sparse$v)),
                          dims=c(simple_triplet_matrix_sparse$nrow, simple_triplet_matrix_sparse$ncol),
                          dimnames = dimnames(simple_triplet_matrix_sparse),
                          giveCsparse = TRUE)
  
}

asSparseMatrix <- function(stm) {
  sparseMatrix(i=stm$i,j=stm$j,x=stm$v, dimnames=dimnames(stm))
}

#long to stm
###given a markov bigram table, compute the markov bigram sparse matrix
long_to_stm <- function(dtm.long, ivar = "i", jvar = "j", var = "x") {
  require(slam)
  stm<- simple_triplet_matrix(i=dtm.long[,ivar],j= dtm.long[,jvar],v= dtm.long[,var], 
                              nrow= nlevels(dtm.long[,ivar]), ncol= nlevels(dtm.long[,jvar]),
                              dimnames = list( ivar= levels(dtm.long[,ivar]), jvar= levels(dtm.long[,jvar]))
  )
}

#factor to SparseMatrix
factor2sparseMatrix <- function(v) {
  sparseMatrix(i=1:length(v),j=as.numeric(v),x=1,
               dims=c(length(v),nlevels(v)),
               dimnames = list(i = names(v), j = levels(v)))
}

#sparse log with Laplacian smoothing
sparselog <- function(sparsem, smoothing =.001) {
  sparsem@x <- log(1+sparsem@x/smoothing)
  sparsem
}