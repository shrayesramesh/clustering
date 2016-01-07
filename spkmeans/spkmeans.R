#spherical kmeans
require(Matrix)
require(slam)

rowNormalize <- function(datamatrix, l2=TRUE) {
  #datamatrix is a sparseMatrix from the Matrix package
  if(l2) {
    row.norm<- as.vector(sqrt( (datamatrix * datamatrix) %*% rep(1,dim(datamatrix)[2])))
  } else {
    row.norm<- as.vector(datamatrix %*% rep(1,dim(datamatrix)[2]))
  }
  D.norm <- Diagonal(dim(datamatrix)[1], 1/row.norm)
  D.norm %*% datamatrix
}

spkmeans<- function(dtm, K, mindelta = 2, clusterid=NULL) {

  
  #initialize
  N <- dim(dtm)[1]
  W <- dim(dtm)[2]
  
  #check class of input data
  #normalize dtm for cosine similarity
  if( length(which(c("DocumentTermMatrix","simple_triplet_matrix") %in% class(dtm))) > 0 ) {
    data<- rowNormalize(sparseMatrix(i=dtm$i,j=dtm$j,x=dtm$v, dimnames=dimnames(dtm)))
  } else {
    data <- rowNormalize(dtm)
  }
  #random start initialization
  if (is.null(clusterid)) {
    clusterid <- sample(1:K,N, replace=TRUE)
  }
  
  delta <- N

  while(delta >mindelta) {
    
    #split data into clusters
    clusters <- split(1:N, clusterid)
    
    #update cluster centroids
    centroids <- do.call(rbind,lapply(clusters, function(docs){
      if(length(docs) == 0) {
        print("empty cluster")
        s<- sample(1:N,K,prob = -(rowSums(similarities)/W))
        colSums(data[s,])
      } else if (length(docs) == 1) {
        data[docs,]
      } else {
        colSums(data[docs,])
      }
    }))
    #row normalize centroids
    centroids <- rowNormalize(centroids)
    
    #compute new clusterids
    similarities <- data %*% t(centroids)
    #for each row, find the cluster with the maximum similarity
    newclusterid <- apply(similarities,1, which.max)
    
    #convergence check
    delta <- length(which(newclusterid != clusterid))
    
    #update clusterids
    clusterid <- newclusterid
    
    print(paste("changes this iteration:", delta))
  }
  list(clusterid = clusterid,
       centroids = centroids
       #totalsimilarity = try(sum(similarities[,clusterid]))
  )
}


unsupervisednb <- function(datamatrix, K, mindelta=2, smoothing = .001) {
  
  N<-dim(datamatrix)[1]
  W<- rowSums(datamatrix)
  #initial color matrix
  clusterid<- sample(1:K,N,replace=TRUE)
  
  delta <- N
  while(delta >mindelta) {
    clusters <- split(1:N, factor(clusterid,levels=1:K))
    
    #update clusters
    centroids <- do.call(rbind,lapply(clusters, function(edges){
      if(length(edges) == 0) {
        print("empty cluster")
        s<- sample(1:N,K,prob = -(rowSums(clusterll)/W))
        colSums(datamatrix[s,])
      } else if (length(edges) == 1) {
        datamatrix[edges,]
      } else {
        colSums(datamatrix[edges,])
      }
    }))
    centroids <- rowNormalize(centroids,l2=FALSE)
    
    clusterll <- datamatrix %*% t(sparselog(centroids,smoothing))
    
    newclusterid <-  apply(clusterll,1, which.max)
    delta <- length(which(newclusterid != clusterid))
    clusterid <- newclusterid
    
    print(paste(delta," changes this iteration."))
    gc()
  }
  list(clusterid = clusterid,
       centroids = centroids
  )
}


#explain clusters
explain.clusters <- function(pr.w.k,nwords=10, idf = rep(1,dim(pr.w.k)[2]), min=.02) {
  
  data<- apply(pr.w.k, 1, function(topic) {
    o<- order(topic*idf, decreasing = TRUE)[1:nwords]
    #o<- order(topic, decreasing=TRUE)[1:15]
    weight <- unname((topic*idf)[o])
    term  <- dimnames(pr.w.k)[[2]][o]                 
    list(term = ifelse(weight>min, term, ""),
         weight = weight
         )
  })
  terms<-  do.call(cbind,lapply(data,function(x) x$term))
  weight<- do.call(cbind,lapply(data,function(x) x$weight))
  
  list(terms = terms,weights = weight)
}

#helper function
asSparseMatrix <- function(stm) {
  sparseMatrix(i=stm$i,j=stm$j,x=stm$v, dimnames=dimnames(stm))
}

#for output

terms.output <- function(topic.terms) {
  K <- dim(topic.terms[[1]])[2]
  base <- data.frame(topic.terms)
  temp <- rep(1:K, each = 2)
  temp[(1:K)*2] <- 1:K + K
  base <- base[,temp]
  base
}

