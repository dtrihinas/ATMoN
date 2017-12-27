
#
# Author: Luis F. Chiroque
# e-mail: lf.chiroque@imdea.org
# IMDEA Networks Institute
#

library("parallel") # mclapply()
library("igraph")
library("lubridate") # hour()

# Variables
## Autonomous System datasets
as.datasets <- list("as.733"=c(folder="as-733/")
                    , "as.15y"=c(folder="as_topology_daily/"))
## Face 2 Face networks
f2f.dataset <- list("sg.09"=list(fileName="listcontacts_2009", sep="\t", timeOffset=20
                                 , filters=list(dow="Monday", open.hour=10, close.hour=18))
                    , "thiersHS.12"=list(fileName="thiers_2012", sep="\t", timeOffset=20
                                         , filters=list(dow=c("Saturday", "Sunday"),  open.hour=6, close.hour=14))
                    , "thiersHS.11"=list(fileName="thiers_2011", sep="\t", timeOffset=10)
                    , "primaryschool"=list(fileName="primaryschool", sep="\t", timeOffset=20)
                    , "highschool.13"=list(fileName="HS2013", sep=" ", timeOffset=20)
                    , "hospital"=list(fileName="hospital", sep="\t", timeOffset=20)
                    )
## Options list
graphOpt <- list()
graphOpt$aggregateNodeSet <- FALSE


# Functions
## graph stream from set of files
get_gstream_from_files <- function(graphFolder, aggregateNodeSet = F) {
  # set of file names
  fileNames <- list.files(paste0(datasets.folderName, graphFolder))
  if (aggregateNodeSet) {
    # compute aggregated node set
    #nodes <- unique(do.call(c, mclapply(fileNames, function(fileName){
    nodes <- unique(do.call(c, lapply(fileNames, function(fileName){
      edges <- read.table(paste0(dirName, graphFolder, fileName))
      union(edges[,1], edges[,2])
    })))
    # compute graph stream
    #gstream <- mclapply(fileNames, function(fileName){
    gstream <- lapply(fileNames, function(fileName){
      edges <- read.table(paste0(dirName, graphFolder, fileName))
      graph.data.frame(edges, directed = T, vertices = data.frame(name=nodes))
    })
  } else {
    # compute graph stream
    #gstream <- mclapply(fileNames, function(fileName){
    gstream <- lapply(fileNames, function(fileName){
      edges <- read.table(paste0(datasets.folderName, graphFolder, fileName))
      graph.data.frame(edges, directed = T)
    })
  }
  # RETURN
  list(gstream=gstream)
}

get_gstream_from_eventList <- function(graphInfo) {
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
  gstream <- mclapply(thrs, function(thrs) {
    #thrs <- min(edges.dyn$start) + i*timeUnit
    filter <- which(edges.dyn$start > thrs & edges.dyn$start <= (thrs + graphInfo$timeOffset))
    graph.data.frame(edges[filter,], directed=F
                     #, vertices=nodes
    )
  })
  list(gstream=gstream)
}


# Code

if ( ! file.exists(Robj.fileName) ) {
  # init the collection list
  gstream.collection <- list()
  
  # add AS gstreams
  #gstream.collection <- c(gstream.collection
  #                           , lapply(as.datasets, get_gstream_from_files
  #                                    , aggregateNodeSet=graphOpt$aggregateNodeSet))
  # add face2face networks
  gstream.collection <- c(gstream.collection
                         , lapply(f2f.dataset, get_gstream_from_eventList))
  
  # write the gstream.collection obj
  save(gstream.collection, file=Robj.fileName)
} else {
  load(Robj.fileName)
}

