
#stscale computes and estimated triangulation using
#data is a NxK matrix
#dpair is a function that takes two (row) subsets of data and returns a distance in R1
#S is a sampling parameter. When S=N, the model is an exact triangulation of the data.
#computing time is O(S*N)
stscale <- function(data, dpair, S=100) {
  N<- dim(data)[1]
  
  #1. sample estimate of radial distance (mean distance to other points)
  #radial is an estimate of the distance from the origin (0,0) to the current point
  sampleset<- sample(1:N,S,replace=TRUE)
  system.time({
    radialset<- dpair(data,data[sampleset,])
    radial<- apply(radialset,1,mean)
  })
  
  #2. flag the point furthest away
  flag1 <- which.max(radial)
  #flag1 <- which.min(radial)
  
  #triangulate 1st coordinate
  system.time({
    A <- radial
    B <- dpair(data,data[flag1,])
    C <- radial[flag1]
    x1<- as.vector( (A*A + C*C - B*B) / (2 * C) )
  })
  
  #find the point "most orthagonal" to the first flag
  flag2 <- which.min(abs(x1))
  
  #triangulate 2nd coordinate
  system.time({
    A <- radial
    B <- dpair(data,data[flag2,])
    C <- radial[flag2]
    x2<- as.vector( (A*A + C*C - B*B) / (2 * C) )
    rm(list = c("A","B","C"))
  })
  projcoords = cbind(x1,x2)
  projcoords
}


###END HELPER FUNCTIONS

###LOAD DATA
