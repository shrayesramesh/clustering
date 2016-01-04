#modfied canopy clustering

#### FUNCTIONS

##FOR INITIAL ASSIGNMENTS, CONSTRUCT A CRUDE BALL-TREE
##at each node of points, sample two points, and assign all points in node to the closest point
##split data by assignment, and recursively apply function to each side of the tree
##stop if node is length mingroupsize
assignbins <- function(currentnode, dpair, mingroupsize=300) {
  groupsize<- length(currentnode)
  if (groupsize< mingroupsize) return(currentnode)
  
  #else, generate subnodes (subgroups) and continue
  centroid_id<- sample(currentnode, 2, replace=FALSE)
  
  #find assignments
  distances <- dpair(data[currentnode,],data[centroid_id,])
  subassignments<- (distances[,1] - distances[,2] >0 )+1
  
  #split node by assignments and recurse to those sections
  currentsplit<- split(currentnode,subassignments)
  
  #recurse
  branches<- lapply(currentsplit,assignbins, dpair, mingroupsize)
  
  branches
}
#function to unnest tree
unnest <- function(x) {
  if(is.null(names(x))) {
    list(unname(unlist(x)))
  }
  else {
    do.call(c, lapply(x, unnest))
  }
}


##FLAG FUNCTION TO GENERATE CLUSTERS FROM DATA
##candidate centroids, "flags" are picked in batches of length between 1 and runsize
##points are then picked off and assigned to flags if it's close enough to any
###if a flag is not part of another flag's cluster, emit the flag and cluster
###if two flags are in each other's clusters, merge clusters and keep unique clusters
##repeat runs until no more points are available to be chosen as flags

flagfunction <-function(currentnode, dpair, runsize=1, radius=3)
{
  #all unassigned points; the goal is to assign all
  comparisonpoints <- currentnode
  
  #pick a random flag order
  permutation <- sample(currentnode) #only need to sample once
  
  #keep track of data
  runcount <- 0
  rundata <- list()
  
  #while unassigned points exist
  while(length(comparisonpoints) > 0 ){
    
    #pick k flags (min = 1, max = runsize)
    k = min(runsize,max(1,length(comparisonpoints)/2))
    currentflag <- permutation[1:k]
    runcount<- runcount+1
    
    #compute distances to flags
    currentdistance <- dpair(data[comparisonpoints,],data[currentflag,])
    
    #produce flag and its cluster
    currentclosebin <- sapply(1:length(currentflag), function(f) {
      #calculate close points and close flags
      closepoints <- comparisonpoints[which(as.vector(currentdistance[,f]) < radius)]
      closeflags<- intersect(closepoints,currentflag)

      #if the only close flag is self, return all the points
      #otherwise, we need to do some duplication/reduction
      if(length(closeflags)<=1) {
        flagdata<- list(flag = currentflag[f],
                     closebin = closepoints
                     #currentflag[-f] %in% comparisonpoints[which(as.vector(currentdistance[,f]) < radius)]
                     )
      } else {
        t<-apply(currentdistance[,which(currentflag %in% closeflags)]<radius,2, function(f) comparisonpoints[which(f)])
        flagdata<- list(flag = closeflags,
                        closebin =  Reduce(union,t)
                        )
      }
      flagdata
    }, simplify=FALSE, USE.NAMES=FALSE)

    #toss all points assigned during this run from permutation and comparisonpoints
    #note, this should get rid of current flag too
    closepoints <- Reduce(union,lapply(currentclosebin,function(x) x$closebin))
    comparisonpoints <- setdiff(comparisonpoints,closepoints)
    permutation <- setdiff(permutation,closepoints)

    #remove duplicate flags (two flags that were in each others' clusters)
    rundata[[runcount]] <- unique(currentclosebin)
  }
  
  #strip away list structure and return full list of flags
  unlist(rundata,recursive=FALSE)
}


##MASTER FLAG FUNCTION TO SPLIT,MERGE,AND CALCULATE FINAL CLUSTERS
##takes in:
#   data row ids to cluster
#   a distance function
#   a radius parameter for excluding close points from forming new clusters
#   runsize1 and runsize2 for the size of flagfunction runs (for speed)
masterflagfunction<- function(dpair,radius, runsize1=1, runsize2=1, batchsize=300, plot.batches=FALSE) {
  
  currentnode = 1:dim(data)
  #split all data into batches
  
  #the nieve way works:
  #batch<-split(currentnode, as.numeric(cut(runif(length(currentnode)),breaks=length(currentnode)/batchsize)))
  
  #but the following is better
  #maximize number of matches within batches to minimize time in the merge stage
  print("generating optimal splits")
  print(system.time({
    batch<- unnest(assignbins(currentnode,dpair,mingroupsize=batchsize))
  }))
  print(paste(length(batch),"splits"))

  #if printing
  if(plot.batches==TRUE) { 
    toplot<- do.call(rbind,
                     lapply(batch,function(c) cbind(matrix(data[c,],ncol=2),min(c)))
                     )
    toplot[,3] <- as.numeric(as.factor((toplot[,3])))
    plot(toplot[,1],toplot[,2],col=toplot[,3],xlab="x1",ylab="x2",main="tree splitting")
  }
  

  #conduct flagfunction on each batch of data
  #this is equivalently a MAP job mapping (batch,ids) to (batch, subflags)
  
  #compute subflags using flagfunction
  print("calculating batch flags")
  print(system.time({
    batchflagdata<- lapply(batch, flagfunction, dpair= dpair, runsize=runsize1,radius= radius)
  }))
  #extact subflag data
  #this is equivalently a REDUCE job unioning subflags into list of subflag ids
  subflagdata <- unlist(batchflagdata,recursive=FALSE)
  rm(batchflagdata)
  subflags <- unlist(lapply(subflagdata, function(sf) sf$flag),use.names = FALSE)
  
  #cluster the subflags into masterflags
  #This is a slow step and cannot be parallelized.It can be made faster with:
  #Longer runsize (to a limit)
  #Better splitting of the data to maximize matches within batch
  print("merging batch flags into master flags")
  print(paste(length(subflagdata),"flags to cluster"))
  #calculate masterflags
  print(system.time({
    masterflagdata<- flagfunction(subflags,dpair=dpair, runsize=runsize2, radius = radius)
  }))
  
  
  print("extracting clusters")
  
  #merge data and generate clusters
  print(system.time({
  masterclusters<- lapply(masterflagdata, function(mf) {
    #find the list of masterflag's subflag data
    mysubflagdata<- subflagdata[which(subflags %in% mf$closebin)]
    
    #calculate union of all of my subflag's closebins
    myclosebins<- Reduce(union,lapply(mysubflagdata, function(mysf) mysf$closebin))
    
    #return the masterflag's cluster
    myclosebins
  })
  }))
  
  print(paste(length(masterclusters),"final clusters"))
  
  print("computing final distance matrix")  
  #recompute distances and rbind together
  print(system.time({
    sparsed<- do.call(rbind,lapply(masterclusters, function(myclosebins){
    #return the pairwise distances for the groups in my masterclosebin (note this is inefficient, but works when groups are small)
    sparsify(dpair(data[myclosebins,],data[myclosebins,]),myclosebins,myclosebins)
  }))
  }))
  
  list(clusters=masterclusters,distances = sparsed)
}

###EXECUTE CODE
###CODE
source("sparsify.R")
source("stscale.R")
source("simulations.R")

#generate data

###
###NOTE:: ALL FUNCTIONS in this file REFER TO DATA object as "data"!
###THIS NEEDS TO BE CHANGED IN THE NEXT VERSION
###
data <- generatedata(4,grid=25)

#demo of fast sampled triangulation layout
par(mfrow=c(1,2)) #for plotting batches
layout<- stscale(data,dpair=dpair,S=100)
plot(data, xlab="x1", ylab="x2", main="real coordinates")
plot(layout,xlab="x1", ylab="x2", main="stscale coordinates")

#comparison of computation of full NxN matrix
# system.time({
# dmat<- dpair(data,data)
# })

system.time({
  all.output<-masterflagfunction(dpair=dpair,radius=.5,runsize1=1,runsize2=1,batchsize=1000, plot.batches=TRUE)
  output<- all.output$distances
  clusters<-all.output$clusters
})

#assign points to clusters for printing
#toplot[,3] contains cluster ids
toplot<- do.call(rbind,
                 lapply(clusters,function(c) cbind(matrix(data[c,],ncol=2),min(c)))
)
toplot[,3] <- as.numeric(as.factor((toplot[,3])))
plot(toplot[,1],toplot[,2],col=toplot[,3],xlab="x1",ylab="x2",main="final clusters")
par(mfrow=c(1,1))
