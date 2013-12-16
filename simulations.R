#generate fake data
# 
# ##METHOD 1:
# #pick three random multivariate normal distributions
# xsample <- list ( g1 = list( xmeans = c(-2,2),
#                              xcovar = rbind( c(1,0),
#                                              c(0,1))
#                              ),
#                   g2 = list( xmeans = c(5,2),
#                              xcovar = rbind( c(1,0),
#                                              c(0,1))
#                              ),
#                   g3 = list( xmeans = c(-3,-4),
#                              xcovar = rbind( c(1,0),
#                                              c(0,1))
#                   )
# )
# 
# require(MASS) #for mvrnorm
# generatedata <- function(N, xpar) {
#   t(replicate(N,  {
#     s<- xpar[[sample(1:length(xpar),1)]]
#     mvrnorm(1,s$xmeans, s$xcovar)
#   }))
# }

##alternate approach to generating data:
##generate 5 clusters with centers aligned to grid points, then add noise
generatedata<- function(groupsize=10, grid=5) {
  Reduce(rbind,sapply(1:grid, function(x1){ 
    Reduce(rbind,sapply( (1:grid), function(x2) {
      mvrnorm(groupsize,c(x1,x2),rbind( c(.01,0),
                                 c(0,.01)))
      },simplify=FALSE))
    },simplify=FALSE))
}

#distance function
d <- function(x1, x2) { sqrt((x1[1]-x2[1])^2 + (x1[2]-x2[2])^2)  }
#pairwise distances (basically, distances from every point in v1 to every point in v2)
#v1 and v2 are each matrices with 2 columns
dpair <- function(v1, v2){
  v1<- matrix(v1,ncol=2)
  v2<- matrix(v2,ncol=2)
  
  as.matrix(apply(v2,1, function(x2) {
    apply(v1,1, function(x1) {
      d(x1,x2)
    })
  }))
}
