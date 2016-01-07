
require(slam)
require(Matrix)

source("C:/Users/sramesh/Documents/GitHub/spkmeans.R")

###HELPER FUNCTIONS

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

###END HELPER FUNCTIONS


#edit this code:
##dtm.long needs a 3 column edgelist:
#two columns (factors) representing the nodes, and one column for the edge weight

dtm.long$node <- factor(dtm.long$node)
dtm.long$node2 <- factor(dtm.long$node2)

#construct sparse adjacency matrix
adj <- as.sparseMatrix(long_to_stm(dtm.long,ivar="node",jvar="node2",var="weight"))



###RUN:

#mindelta: stop/converge if the number of points changing clusters per iteration drops below mindelta
skm = spkmeans(adj,numclusters, mindelta)
#smoothing is laplacian smoothing. set it = 1
unb = unsupervisednb(adj, numclusters, mindelta, smoothing)
