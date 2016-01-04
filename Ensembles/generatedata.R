##new smiley face generator
#generate fake data
nface = 6000
neye1 = 4000
neye2 = 4000
nmouth = 4000
noutliers = 2000

draw_circle <- function(n,center = c(0,0),maxradius = 1) {
  radius <- runif(n,0,maxradius)
  theta <- runif(n,0,2*pi)
  cbind(center[1] + sqrt(radius) * cos(theta),center[2] + sqrt(radius) * sin(theta))
}

draw_arc <- function(n,thetalimits = c(1.25* pi, 1.75*pi ), radius = .35, jitter = .15) {
  radius <- runif(n,radius,radius+jitter)
  theta <- runif(n, thetalimits[1], thetalimits[2])
  cbind(sqrt(radius) * cos(theta),sqrt(radius) * sin(theta))
}

draw_unif <- function(n, radius = 1) {
  cbind( runif(n, -radius, radius), runif(n, -radius, radius) )
}

dataset <- rbind(draw_circle(nface),
                 draw_circle(neye1,c(-.3,.2),.05),
                 draw_circle(neye1,c(.3,.2),.05),
                 draw_arc(nmouth),
                 draw_unif(noutliers))

labels = factor(
  c(rep("face",nface),
    rep("eye1",neye1),
    rep("eye2",neye2),
    rep("mouth",nmouth),
    rep("outliers",noutliers)
    ))

plot(dataset,col=labels,cex=.1)



###TRY SPIRAL DATASET
#dataset <- rbind(draw_arc(2000, c(0,2*pi), 6 , .005),
#                 draw_arc(2000, c(0,2*pi), 15 , .005))
#plot(dataset)






require(slam)
require(Matrix)
#library(methods)
source("C:/Users/sramesh/Documents/GitHub/helper.R") #helper functions


##### MATRIX HELPER FUNCTIONS / TYPE CONVERSIONS #####
factor2sparseMatrix <- function(v) {
  sparseMatrix(i=1:length(v),j=as.numeric(v),x=1,
               dims=c(length(v),nlevels(v)),
               dimnames = list(i = names(v), j = levels(v)))
}


as.sparseMatrix <- function(simple_triplet_matrix_sparse) {
  retval <-  sparseMatrix(i=as.numeric(simple_triplet_matrix_sparse$i),
                          j=as.numeric(simple_triplet_matrix_sparse$j),
                          x=as.numeric(as.character(simple_triplet_matrix_sparse$v)),
                          dims=c(simple_triplet_matrix_sparse$nrow, simple_triplet_matrix_sparse$ncol),
                          dimnames = dimnames(simple_triplet_matrix_sparse),
                          giveCsparse = TRUE)
  
}

#long to stm
long_to_stm <- function(dtm.long, ivar = "i", jvar = "j", var = "x") {
  require(slam)
  stm<- simple_triplet_matrix(i=dtm.long[,ivar],j= dtm.long[,jvar],v= dtm.long[,var], 
                              nrow= nlevels(dtm.long[,ivar]), ncol= nlevels(dtm.long[,jvar]),
                              dimnames = list( ivar= levels(dtm.long[,ivar]), jvar= levels(dtm.long[,jvar]))
  )
}

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


#######START HERE#####

dtm <- data.matrix(dataset)

###PART 1
N <- dim(dtm)[1]
D <- dim(dtm)[2]
datamatrix <- dtm

H <- 80
#png(file=paste("smiley_projections.png"),height=10, width=10, units="in", res=600)
#par(mfrow=c(2,2))
iterationoutput<- sapply(1:H, function(x) {
  print(paste("iteration",x))
  #if sampling, sample here
  s <- 1:N
  #draw random projection
  proj <- matrix(rnorm(D*1,0,1),nrow=D)
  projdata <- datamatrix %*% proj
  K <- 20
  cl <- kmeans(projdata, K)
  #plot(datamatrix[s,],col=cl$cluster, main = paste("iteration",x), xlab= 'x', ylab='y', cex=.1)
  data.frame(partid = x,
             docid = s,
             kid = cl$cluster, row.names=NULL)
}, simplify=FALSE)
#par(mfrow=c(1,1))
#dev.off()
partitions <- do.call(rbind,iterationoutput)


##PART 2
partitions$docid <- factor(partitions$docid,levels=1:N)
partitions$hyperedge <- factor(paste(partitions$partid,partitions$kid))
iterationclusters <- long_to_stm(partitions,ivar="hyperedge",jvar="docid",var="hyperedge")
iterationexists <- as.sparseMatrix(long_to_stm(partitions,ivar="hyperedge",jvar="docid",var="kid")) >0



spectral <- function(m,k) {
  print("constructing laplacian")
  diag(m) <- 0
  Dsqrt = Diagonal(ncol(m),1/sqrt(rowSums(m)))
  m = as.matrix(Dsqrt %*% m %*% Dsqrt)
  print("computing eigenvectors")
  ev <- eigen(m, symmetric=TRUE)
  topk <- ev$vectors[,1:k]
  topk = apply(topk, 2, function(x) (x-mean(x))/sd(x))
  cl <- kmeans(topk,k,nstart=5)$cluster
}

hypergraph <- iterationexists %*% t(iterationexists)
clustercounts <- diag(hypergraph)
norm <- sqrt(clustercounts) %o% sqrt(clustercounts)
lift <- hypergraph/norm

uniongraph <- outer(clustercounts,clustercounts,function(X,Y) X + Y) -hypergraph
jaccard <- hypergraph / uniongraph

#output number of clusters
k<- 4

##schema:
##clusterlabel, x, y

outmat <- do.call(cbind,lapply (list(hyp=hypergraph), function(g){
  #dmat <- as.matrix(1-g)
  dmat <- exp(-g)
  ##mds scale layout for coordinates
  l4p <- cmdscale(dmat)*100
  affinity = g
  ##spectral clustering for metaclusters
  cl <- spectral(affinity,k)
  cbind(cl,l4p)
}))

#fast matrix math to do aggregations and voting
metalabel.sizes <- rowSums(iterationexists)
metalabel.labels <- factor(outmat[,1],levels=1:k)
metalabel.assign <- factor2sparseMatrix(metalabel.labels)
votematrix = t(metalabel.assign) %*% iterationexists
ensemblelabels = apply(votematrix,2,which.max)
entropy = apply(votematrix,2, function(x) -sum(x/sum(x)*log(x/sum(x)),na.rm = TRUE))
colorvar.entropy <- convertcolor(entropy,ncolors=100, start="blue",end="yellow")

density = as.vector(log(metalabel.sizes) %*% iterationexists)
colorvar.density <- convertcolor(-density,ncolors=100, start="blue",end="yellow")

metalabel.coords <- outmat[,2:3]
ensemble.coords = t(t(metalabel.coords) %*% iterationexists)/H

#raw data with labels
#png(file=paste("smiley_labels.png"),height=6, width=10, units="in", res=800)
par(mfrow=c(1,2))
plot(dataset,col=labels,cex=.1, xlab="smiley x", ylab ="smiley y", main = "colored by ground truth")
plot(dataset,col=ensemblelabels,cex=.1, xlab="smiley x", ylab ="smiley y", main = "colored by ensemble labels")
par(mfrow=c(1,1))
#dev.off()

#png(file=paste("smiley_outliers.png"),height=6, width=10, units="in", res=800)
par(mfrow=c(1,2))
plot(dataset,col=colorvar.entropy,cex=.1, xlab="smiley x", ylab ="smiley y", main = "colored by entropy")
plot(dataset,col=colorvar.density,cex=.1, xlab="smiley x", ylab ="smiley y", main = "colored by density")
par(mfrow=c(1,1))
#dev.off()

#png(file=paste("smiley_coords.png"),height=6, width=10, units="in", res=800)
par(mfrow=c(1,2))
plot(ensemble.coords[,1],ensemble.coords[,2],col=labels,cex=.1, xlab = "ensemble coord x", ylab= "ensemble coord y", main = "colored by ground truth")
plot(ensemble.coords[,1],ensemble.coords[,2],col=ensemblelabels,cex=.1, xlab = "ensemble coord x", ylab= "ensemble coord y", main = "colored by ensemble labels")
par(mfrow=c(1,1))
#dev.off()




##metaclustering matrix
#png(file=paste("metaclustering_pre.png"),height=5, width=5, units="in", res=600)
#image(hypergraph)
#dev.off()
#png(file=paste("metaclustering_post.png"),height=5, width=5, units="in", res=600)
#image(hypergraph[order(outmat[,4]),order(outmat[,4])])
#dev.off()


##plotting points
#plist <- c(4300, 10000, 15000)
#png(file=paste("smiley_layoutexplain.png"),height=9, width=6, units="in", res=400)
#par(mfrow=c(3,2))
# sapply(plist,function(p) {
#   plot(datamatrix, xlab="", ylab="",cex=.1, main="point in original data")
#   points(datamatrix[p,1], datamatrix[p,2], pch="X", col="blue", cex=2)
#   hyperedges <- which(iterationexists[,p])
#   
#   plot(metalabel.coords,cex=.3, xlab = "", ylab="", main = "ensemble coordinates")
#   points(metalabel.coords[hyperedges,], pch="X", col="red",cex=1)
#   segments(metalabel.coords[hyperedges,1],metalabel.coords[hyperedges,2],
#            ensemble.coords[p,1],ensemble.coords[p,2],col="yellow",cex=.1)
#   points(ensemble.coords[p,1],ensemble.coords[p,2],col="blue",pch="X", cex=2)
#   TRUE
# })
# par(mfrow=c(1,1))
# dev.off()