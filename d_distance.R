
#
# R code to compute 'Dissimilarity' metric
# source: https://github.com/tischieber/Quantifying-Network-Structural-Dissimilarities
# paper: http://www.nature.com/articles/ncomms13928
#

library(igraph)
library(Matrix)
#library('MatrixStats')


# Shannon entropy

entropia <- function(a) {
  a <- a[which(a > 0)]
  -sum(a * log(a))
}


#returns the node distance matrix

node_distance <- function(g) {
  n <- vcount(g)
  retorno <- NA
  if( n == 1 ) {
    retorno = 1
  }
  if( n > 1 ) {
    a <- Matrix(0, nrow=n, ncol=n, sparse=TRUE)
    m <- shortest.paths(g, algorithm=c("unweighted"))
    m[which(m == "Inf")] <- n
    quem <- setdiff(intersect(m,m),0)
    for(j in (1:length(quem))) {
      l <- which(m==quem[j]) / n
      linhas <- floor(l) + 1
      posicoesm1 <- which(l == floor(l))
      if(length(posicoesm1) > 0) {
        linhas[posicoesm1] <- linhas[posicoesm1] - 1
      }
      a[1:n, quem[j]] <- hist(linhas, plot=FALSE, breaks=(0:n))$counts
    }
    #m<-c()
    retorno = a / (n-1)
  }
  return(retorno)
}


# nnd

nnd <- function(g) {
  N <- vcount(g)
  nd <- node_distance(g)
  pdfm <- colMeans(as.matrix(nd))
  norm <- log(max(c(2, length(which(pdfm[1:(N-1)] > 0)) + 1)))
  return(c(pdfm, max(c(0, entropia(pdfm) - entropia(nd) / N)) / norm))
}

#function

alpha <- function(g) {
  N <- vcount(g)
  if( N == 0 )
    return(c(0,1))
  r <- sort(alpha.centrality(g, exo=degree(g) / (N - 1), alpha=1 / N)) / ((N^2))
  return(c(r, max(c(0, 1 - sum(r)))))
}

#function

D <- function(g, h, w1, w2, w3) {
  first <- 0
  second <- 0
  third <- 0
  g <- read.graph(g, format=c("edgelist"), directed=FALSE)
  h <- read.graph(h, format=c("edgelist"), directed=FALSE)
  N <- vcount(g)
  M <- vcount(h)
  PM <- matrix(0, ncol=max(c(M, N)))
  if(w1 + w2 > 0){
    pg=nnd(g)
    PM[1:(N-1)] = pg[1:(N-1)]
    PM[length(PM)] <- pg[N]
    ph = nnd(h)
    PM[1:(M-1)] = PM[1:(M-1)] + ph[1:(M-1)]
    PM[length(PM)] <- PM[length(PM)] + ph[M]
    PM <- PM / 2
    first <- sqrt(max(c((entropia(PM) - (entropia(pg[1:N]) + entropia(ph[1:M])) / 2) / log(2), 0)))
    second <- abs(sqrt(pg[N+1]) - sqrt(ph[M+1]))
  }
  if(w3 > 0){
    pg <- alpha(g)
    ph <- alpha(h)
    m <- max(c(length(pg), length(ph)))
    Pg <- matrix(0, ncol=m)
    Ph <- matrix(0,ncol=m)
    Pg[(m - length(pg) + 1):m] <- pg
    Ph[(m - length(ph) + 1):m] <- ph
    third <- third + sqrt((entropia((Pg + Ph) / 2) - (entropia(pg) + entropia(ph)) / 2) / log(2)) / 2
    g <- graph.complementer(g)
    h <- graph.complementer(h)
    pg <- alpha(g)
    ph <- alpha(h)
    m <- max(c(length(pg), length(ph)))
    Pg <- matrix(0, ncol=m)
    Ph <- matrix(0, ncol=m)
    Pg[(m - length(pg) + 1):m] <- pg
    Ph[(m - length(ph) + 1):m] <- ph
    third <- third + sqrt((entropia((Pg + Ph) / 2) - (entropia(pg) + entropia(ph)) / 2) / log(2)) / 2
  }
  
  return(w1 * first + w2 * second + w3 * third)
}

#function

d <- function(g, h, w1, w2, w3) {
  first <- 0
  second <- 0
  third <- 0
  g <- read.graph(g, format=c("edgelist"), directed=FALSE)
  h <- read.graph(h, format=c("edgelist"), directed=FALSE)
  N <- vcount(g)
  M <- vcount(h)
  PM <- matrix(0, ncol=max(c(M, N)))
  if( w1 + w2 > 0){
    pg = nnd(g)
    PM[1:(N-1)] = pg[1:(N-1)]
    PM[length(PM)] <- pg[N]
    ph = nnd(h)
    PM[1:(M-1)] = PM[1:(M-1)] + ph[1:(M-1)]
    PM[length(PM)] <- PM[length(PM)] + ph[M]
    PM <- PM / 2
    first <- sqrt(max(c((entropia(PM) - (entropia(pg[1:N]) + entropia(ph[1:M])) / 2) / log(2), 0)))
    second <- abs(sqrt(pg[N+1]) - sqrt(ph[M+1]))
  }
  
  if(w3 > 0){
    pg <- alpha(g)
    ph <- alpha(h)
    m <- max(c(length(pg), length(ph)))
    Pg <- matrix(0, ncol=m)
    Ph <- matrix(0, ncol=m)
    Pg[(m - length(pg) + 1):m] <- pg
    Ph[(m - length(ph) + 1):m] <- ph
    third <- third + sqrt((entropia((Pg + Ph) / 2) - (entropia(pg) + entropia(ph)) / 2) / log(2))
  }
  return(w1 * first + w2 * second + w3 * third)
}


# Function: centrality term

central.term <- function(g, h, third) {
  pg <- alpha(g)
  ph <- alpha(h)
  m <- max(c(length(pg), length(ph)))
  Pg <- matrix(0, ncol=m)
  Ph <- matrix(0, ncol=m)
  Pg[(m - length(pg) + 1):m] <- pg
  Ph[(m - length(ph) + 1):m] <- ph
  third <- third + sqrt((entropia((Pg + Ph) / 2) - (entropia(pg) + entropia(ph)) / 2) / log(2)) / 2
}

# Function: dissimilarity joint function

dissimilarity <- function(g, h, w1=.45, w2=.45, w3=.1, use.complement=FALSE) {
  g <- simplify(g)
  h <- simplify(h)
  first <- 0
  second <- 0
  third <- 0
  N <- vcount(g)
  M <- vcount(h)
#print(paste0("N=", N, " M=", M))
  PM <- matrix(0, ncol=max(c(M, N, 1)))
  if(w1 + w2 > 0){
    if( N == 0 || M == 0 ){
      zero.result <- 1
      if( N + M == 0 ){
        zero.result <- 0
      }
      first <- second <- zero.result
    } else {
      pg <- nnd(g)
      PM[1:(N-1)] <- pg[1:(N-1)]
      PM[length(PM)] <- pg[N]
      ph <- nnd(h)
      PM[1:(M-1)] <- PM[1:(M-1)] + ph[1:(M-1)]
      PM[length(PM)] <- PM[length(PM)] + ph[M]
      PM <- PM / 2
      first <- sqrt(max(c((entropia(PM) - (entropia(pg[1:N]) + entropia(ph[1:M])) / 2) / log(2), 0)))
      second <- abs(sqrt(pg[N+1]) - sqrt(ph[M+1]))
    }
  }
  if(w3 > 0){
    third <- central.term(g, h, third)
    if( use.complement ){
      g <- graph.complementer(g)
      h <- graph.complementer(h)
      third <- central.term(g, h, third)
    }
  }
  return(w1 * first + w2 * second + w3 * third)
}


#rm(list=setdiff(ls(),lsf.str()))

# g3 <- graph.graphdb("si6_r005_s100.B99") database de graph isomorphism 
#write(
