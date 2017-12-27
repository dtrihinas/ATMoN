
#
# Author: Luis F. Chiroque
# e-mail: lf.chiroque@imdea.org
# IMDEA Networks Institute
#

library("igraph") # graph metrics

source("d_distance.R") # dissmilarity metric


if (!exists(gstreamObjectName))
  stop(paste0("'", gstreamObjectName, "' object not loaded"))

# lobby rank or mu-PCI (mu=1)
#   muPCI(v) = k, such that mu*k nodes in the mu-hop neighbourhood
# has at least degree k
muPCI.rank <- function(graph, vids=V(graph), mu=1) {
  graph <- simplify(graph)
  sapply(vids, function(v) {
    nbhd <- neighborhood(graph, mu, v, mode="out")[[1]][-1]
    nbhd.deg <- sort(degree(graph, nbhd, mode="out"), decreasing=T)
    v.deg <- length(nbhd)
    int.div <- v.deg %/% mu
    k.set <- 1:int.div
    table(factor(nbhd.deg[k.set * mu] >= k.set, levels=c(T,F)))["TRUE"]
  })
}
muPCI.rank.alt <- function(graph, vids=V(graph), mu=1) {
  graph <- simplify(graph)
  nbhds <- neighborhood(graph, mu, mode="out")
  sapply(nbhds, function(nbhd) {
    #nbhd <- neighborhood(graph, mu, v, mode="out")[[1]][-1]
    v.deg <- nbhd[1]
    nbhd <- nbhd[-1]
    nbhd.deg <- sort(degree(graph, nbhd, mode="out"), decreasing=T)
    v.deg <- length(nbhd)
    int.div <- v.deg %/% mu
    k.set <- 1:int.div
    table(factor(nbhd.deg[k.set * mu] >= k.set, levels=c(T,F)))["TRUE"]
  })
}

# Function: given a datastream object (ds.obj) compute its list of metrics

compute_gstream_metrics <- function(ds.obj) {
  if (is.null(ds.obj$metrics))
    ds.obj$metrics <- list()
  # Regular metrics
  if (! "diameter" %in% names(ds.obj$metrics))
    ds.obj$metrics$diameter <- sapply(ds.obj$gstream, diameter)
  if (! "pageRank" %in% names(ds.obj$metrics))
    ds.obj$metrics$pageRank <- sapply(ds.obj$gstream, function(gg) page.rank(gg)$vector)
  if (! "outDegree" %in% names(ds.obj$metrics))
    ds.obj$metrics$outDegree <- sapply(ds.obj$gstream, function(gg) degree(gg, mode='out'))
  if (! "betweenness" %in% names(ds.obj$metrics))
    ds.obj$metrics$betweenness <- sapply(ds.obj$gstream, betweenness)
  if (! "l.rank" %in% names(ds.obj$metrics))
    ds.obj$metrics$l.rank <- sapply(ds.obj$gstream, muPCI.rank)
  if (! "cluster.mem" %in% names(ds.obj$metrics))
    ds.obj$metrics$cluster.mem <- sapply(ds.obj$gstream, function(gg) clusters(gg)$membership)
  if (! "nofComponents" %in% names(ds.obj$metrics))
    ds.obj$metrics$nofComponents <- sapply(ds.obj$gstream, function(gg) components(gg)$no)
  if( ! "max.clique" %in% names(ds.obj$metrics))
    ds.obj$metrics$max.clique <- sapply(ds.obj$gstream, function(gg) {
      max(unlist(sapply(max_cliques(gg), length)), 0, na.rm = TRUE)
    })
  if( ! "deg.dist.pow" %in% names(ds.obj$metrics))
    ds.obj$metrics$deg.dist.pow <- sapply(ds.obj$gstream, function(gg) {
      if( vcount(gg) == 0 ) {
        res <- 0
      }else{
        res <- power.law.fit(degree_distribution(gg, cumulative = T, mode="out"))$alpha
      }
      if (res == Inf )
        res <- 0
      res
    })
  if (! "vsize" %in% names(ds.obj$metrics))
    ds.obj$metrics$vsize <- sapply(ds.obj$gstream, vcount)
  if (! "esize" %in% names(ds.obj$metrics))
    ds.obj$metrics$esize <- sapply(ds.obj$gstream, ecount)
  if (! "gcRatio" %in% names(ds.obj$metrics))
    ds.obj$metrics$gcRatio <- sapply(ds.obj$gstream, function(gg) {
      cc <- components(gg)
      ifelse(cc$no > 0
             , table(cc$membership)[as.character(which.max(cc$csize))] / vcount(gg)
             , 0)
    })
  # average metrics
  if (! "avgPageRank" %in% names(ds.obj$metrics))
    ds.obj$metrics$avgPageRank <- sapply(ds.obj$metrics$pageRank, mean)
  if (! "avgOutDegree" %in% names(ds.obj$metrics))
    ds.obj$metrics$avgOutDegree <- sapply(ds.obj$metrics$outDegree, mean)
  # Topology Metrics
  if (! "assortativity" %in% names(ds.obj$metrics))
    ds.obj$metrics$assortativity <- sapply(ds.obj$gstream, assortativity_degree, directed=T)
  if (! "transitivity" %in% names(ds.obj$metrics))
    ds.obj$metrics$transitivity <- sapply(ds.obj$gstream, transitivity, type="global")
  if( !"dissimilarity" %in% names(ds.obj$metrics) )
    ds.obj$metrics$dissimilarity <- c(0, mapply(dissimilarity
                                          , ds.obj$gstream[-length(ds.obj$gstream)]
                                          , ds.obj$gstream[-1]
                                          , use.complement=FALSE))
  
  # RETURN
  ds.obj
}


# Code
## compute the metrics for each data stream
gstream.collection <- lapply(gstream.collection, compute_gstream_metrics)
## save the data stream with the computed metrics
save(gstream.collection, file=Robj.fileName)
