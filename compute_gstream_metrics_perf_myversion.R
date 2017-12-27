
compute_gstream_metrics <- function(g.graph, metric.name, h.graph=g.graph) {
  tstart <- proc.time()

  # Regular metrics
  if ( metric.name == "diameter" )
    diameter(g.graph)
  else if ( metric.name == "pageRank" )
    page.rank(g.graph)$vector
  else if ( metric.name == "outDegree" )
    degree(g.graph, mode='out')
  else if ( metric.name == "betweenness" )
    betweenness(g.graph)
  else if ( metric.name == "l.rank" )
    muPCI.rank(g.graph)
  else if ( metric.name == "cluster.mem" )
    clusters(g.graph)$membership
  else if ( metric.name == "nofComponents" )
    components(g.graph)$no
  else if(  metric.name == "max.clique" )
    max(unlist(sapply(max_cliques(g.graph), length)), 0, na.rm = TRUE)
  else if( metric.name == "deg.dist.pow" ) {
    #if( vcount(gg) != 0 ) {
    #  res <- power.law.fit(degree_distribution(gg, cumulative = T, mode="out"))$alpha
    #}
    if( vcount(g.graph) != 0 ) {
        res <- power.law.fit(degree_distribution(g.graph, cumulative = T, mode="out"))$alpha
    }
  } else if ( metric.name == "vsize" )
    vcount(g.graph)
  else if ( metric.name == "esize" )
    ecount(g.graph)
  else if ( metric.name == "gcRatio" ) {
      cc <- components(g.graph)
      ifelse(cc$no > 0
             , table(cc$membership)[as.character(which.max(cc$csize))] / vcount(g.graph)
             , 0)
  # average metrics
  } else if ( metric.name == "avgPageRank" )
    mean(page.rank(g.graph)$vector)
  else if ( metric.name == "avgOutDegree" )
    mean(degree(g.graph, mode='out'))
  # Topology Metrics
  else if ( metric.name == "assortativity" )
    assortativity_degree(g.graph, directed=T)
  else if ( metric.name == "transitivity" )
    transitivity(g.graph, type="global")
  else if( metric.name == "dissimilarity" )
    dissimilarity(g.graph, h.graph, use.complement=FALSE)
  else
    warning(paste0("unrecognized metric: ", metric.name))
  
  #some disk io for emulator
  write.graph(g.graph, "mygraph")
  fz <- file.size("mygraph")
  
  tend <- proc.time() - tstart

  # RETURN
  return(list(t=tend['elapsed'], v=vcount(g.graph), e=ecount(g.graph), fz=fz))
}

run_performance_test_loading_graph <- function(graph.name, metric.name) {
  graphInfo <- f2f.dataset[[graph.name]]
  
  dd <- read.csv(paste0(datasets.folderName, graphInfo$fileName, ".csv"), sep=graphInfo$sep, header=F)
  names(dd) <- c("time", "a", "b", "a_class", "b_class")[1:length(dd)]
  md <- read.table(paste0(datasets.folderName, "metadata_", graphInfo$fileName, ".txt"))
  names(md) <- c("nodes", "class", "gender")[1:length(md)]
  
  #nodes <- data.frame(a=md$nodes, b=md$nodes) # id and labels match
  edges <- dd[,c("a", "b")]
  #edges <- data.frame(t(apply(edges,1,function(pair) if(pair[1] < pair[2]) pair else rev(pair))))
  edges.dyn <- data.frame(start=(dd$time - graphInfo$timeOffset), end=dd$time)
  
  #dup.filter <- !duplicated(cbind(edges, edges.dyn))
  #dup.filter <- !duplicated(edges.dyn)
  
  snapsSize <- floor(diff(range(dd$time)) / graphInfo$timeOffset)
  thrs <- min(edges.dyn$start) + 0:snapsSize * graphInfo$timeOffset
  if ( !is.null(graphInfo$filters) ) {
    zz <- as.POSIXct(thrs, tz="GMT", origin="1970-01-01")
    dow <- weekdays(zz)
    hour <- hour(zz)
    # opening hours: Tue-Sun 09-19
    sample <- which(! dow %in% graphInfo$filters$dow
                    & hour >= graphInfo$filters$open.hour
                    & hour < graphInfo$filters$close.hour)
    thrs <- thrs[sample]
  }
  lapply(thrs, function(thrs) {
    #thrs <- min(edges.dyn$start) + i*timeUnit
    filter <- which(edges.dyn$start > thrs & edges.dyn$start <= (thrs + graphInfo$timeOffset))
    gg <- graph.data.frame(edges[filter,], directed=F
                     #, vertices=nodes
    )
    compute_gstream_metrics(gg, metric.name)
    invisible()
  })
  
  return(0)
}

# gstream.collection must exist!
run_performance_test <- function(graph.name, metric.name)
{
  gstream <- gstream.collection[[graph.name]]$gstream
  
  print(paste0("graph instances: ", length(gstream)))

  res <- mapply(compute_gstream_metrics
          , gstream
          , metric.name
          , c(gstream[1], gstream[-length(gstream)]))

  t_list <- unlist(res)[c(TRUE, FALSE, FALSE, FALSE)]
  #print(t_list)
  avg_t <- mean(t_list)
  
  v_list <- unlist(res)[c(FALSE, TRUE, FALSE, FALSE)]
  #print(v_list)
  
  e_list <- unlist(res)[c(FALSE, FALSE,TRUE, FALSE)]
  #print(e_list)
  
  static_payload = 192 #assuming json monitoring struct
  edge_payload = 74 #assuming we have an edgelist per node
  ginst_payload <- v_list*(static_payload + e_list * edge_payload)
  #print(ginst_payload)
  total_payload <- sum(ginst_payload)
  
  fsto <- sum(unlist(res)[c(FALSE, FALSE, FALSE, TRUE)])
  
  print(paste0("avg metric cc time: ", avg_t, ", total ingress load (MB): ", total_payload/1000000, ", cloud storage low est (MB): ", fsto/1000000))#, ", avg nodes: ", mean(res$v), ", avg edges", mean(res$e)))
  
  return(0)
}


# example
## load gstream.collection
source("load_configuration.R")
source("load_graphStreams.R")

#run_performance_test("thiersHS.12", "pageRank")
#run_performance_test_loading_graph("thiersHS.12", "pageRank")

run_performance_test_loading_graph("sg.09", "diameter")
#run_performance_test("sg.09", "diameter")
