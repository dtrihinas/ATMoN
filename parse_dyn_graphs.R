
# R-code to parse temporal graph files

library("rgexf")
library("igraph")
library("lubridate")
library("parallel")

source("load_configuration.R")


# as.uniqueDay
ts2uniqueDay <- function(z) paste0(year(z), month(z, label=T), day(z))

##############################
### SocioPatterns datasets ###
##############################

datasetNames <- c("primaryschool", "HS2013", "thiers_2012", "thiers_2011", "hospital", "listcontacts_2009")
ds.sep <- c("\t", " ", "\t", "\t", "\t", "\t")
timeOffset <- c(20, 20, 20, 10, 20, 20)

ds.id <- 3 # choose one id within datasetNames

dd <- read.csv(paste0(dataFolder, datasetNames[ds.id], ".csv"), sep=ds.sep[ds.id], header=F)
names(dd) <- c("time", "a", "b", "a_class", "b_class")[1:length(dd)]
md <- read.table(paste0(dataFolder, "metadata_", datasetNames[ds.id], ".txt"))
names(md) <- c("nodes", "class", "gender")[1:length(md)]

nodes <- data.frame(a=md$nodes, b=md$nodes) # id and labels match
edges <- dd[,c("a", "b")]
#edges <- data.frame(t(apply(edges,1,function(pair) if(pair[1] < pair[2]) pair else rev(pair))))
edges.dyn <- data.frame(start=dd$time-timeOffset[ds.id], end=dd$time)

# visualize
zz <- as.POSIXct(edges.dyn$start, tz="GMT", origin="1970-01-01")
hist(zz, "hours", format="%H")

timeUnit <- timeOffset[ds.id]

max.slots <- floor(diff(range(dd$time))/timeUnit)
thrs <- min(edges.dyn$start) + 0:max.slots * timeUnit
zz <- as.POSIXct(thrs, tz = "GMT", origin="1970-01-01")
dow <- weekdays(zz)
hour <- hour(zz)
# DS 6; opening hours: Tue-Sun 09-19
sample <- which(dow != "Monday" & hour >= 10 & hour < 18)
# DS 3; opening hours: Mon-Fry 05-18
sample <- which(!dow %in% c("Saturday", "Sunday") & hour >= 5 & hour < 16)
gstream <- lapply(thrs[sample], function(thrs){
  #thrs <- min(edges.dyn$start) + i*timeUnit
  filter <- which(edges.dyn$start >= thrs & edges.dyn$start < (thrs + timeOffset[ds.id]))
  graph.data.frame(edges[filter,], directed=F
                   #, vertices=nodes
  )
})
metric.vsize <- sapply(gstream,vcount)
metric.diameter <- sapply(gstream, diameter)
metric.diameter <- do.call(c, mclapply(gstream, diameter))

#### daily gstream
day.granularity <- factor(ts2uniqueDay(zz))
filters <- split(1:nrow(edges), day.granularity)
fil.ord <- order(sapply(filters, head, n=1))
## DS 6; opening hours: Tue-Sun 09-19
#dow <- weekdays(zz[sapply(filters, head, n=1)])
#sample <- which(dow != "Monday")
gstream <- lapply(filters[fil.ord], function(filter){
  graph.data.frame(edges[filter,], directed=F
                   #, vertices=nodes
  )
})


###########################
### MIT mining datasets ###
###########################

konect.mit.dsInfo <- list(
  list(fileName="reality-mit.txt", timeUnit=60*10, max.volume=1)
  , list(fileName="contact-mit.txt", timeUnit=60*2, max.volume=.8)
)

dataset.info <- konect.mit.dsInfo[[2]]

edges.dyn <- read.table(paste0(datasets.folderName, dataset.info$fileName)
                        , comment.char = "%")
names(edges.dyn) <- c("a", "b", "w", "end")

# visualize
zz <- as.POSIXct(edges.dyn$end, tz="GMT", origin="1970-01-01")
hist(zz, "hours")

edges.dyn$start <- edges.dyn$end - dataset.info$timeUnit
edges <- edges.dyn[,c(1,2)]
snapsSize <- with(edges.dyn
                  , floor(diff(range(c(start, end))) / dataset.info$timeUnit * dataset.info$max.volume))
thrs <- min(edges.dyn$start) + 0:snapsSize * dataset.info$timeUnit
gstream <- lapply(thrs, function(thrs){
  #thrs <- min(edges.dyn$start) + i*timeUnit
  filter <- which(edges.dyn$start > thrs & edges.dyn$start <= (thrs + dataset.info$timeUnit))
  graph.data.frame(edges[filter,], directed=F
                   #, vertices=nodes
  )
})
vcount <- sapply(gstream, vcount)
ds.obj <- compute_gstream_metrics(list(gstream=gstream))


########################
### MIT sms networks ###
########################

library("XML")

datasets.folderName <- "datasets/"
gephi.fileName <- "SMSdataset_graph.gexf"

xdoc <- xmlTreeParse(paste0(datasets.folderName, gephi.fileName))
top <- xmlRoot(xdoc)
edges <- top[["graph"]][["edges"]]
edges.dyn <- data.frame(t(sapply(xmlChildren(edges), function(edge) xmlAttrs(edge)[-1])), stringsAsFactors = F)
edges.dyn <- data.frame(sapply(edges.dyn, as.numeric))

date.filter <- as.numeric(strptime("2008-10-01", format = "%Y-%m-%d", tz = "GMT"))
subdata <- which(edges.dyn$start >= date.filter)
edges.dyn <- edges.dyn[subdata,]
edges <- edges.dyn[,c(1,2)]

zz <- as.POSIXct(edges.dyn$start, tz="GMT", origin="1970-01-01")

timeUnit <- 360 # *10*24
edges.supp <- with(edges.dyn
                   , min(start) + 0:floor(diff(range(c(start, end)))/timeUnit) * timeUnit)
#zz <- as.POSIXct(thrs, tz = "GMT", origin="1970-01-01")
gstream <- lapply(edges.supp, function(thrs){
  #thrs <- min(edges.dyn$start) + i*timeUnit
  filter <- which(edges.dyn$start >= thrs & edges.dyn$start < (thrs + timeUnit))
  graph.data.frame(edges[filter,], directed=F
                   #, vertices=nodes
  )
})
metric.vsize <- sapply(gstream, vcount)
ds.obj <- compute_gstream_metrics(list(gstream=gstream))

#### daily gstream
day.granularity <- factor(ts2uniqueDay(zz))
filters <- split(1:nrow(edges.dyn), day.granularity)
fil.ord <- order(sapply(filters, head, n=1))
## DS 6; opening hours: Tue-Sun 09-19
#dow <- weekdays(zz[sapply(filters, head, n=1)])
#sample <- which(dow != "Monday")
gstream <- lapply(filters[fil.ord], function(filter){
  graph.data.frame(edges[filter,], directed=F
                   #, vertices=nodes
  )
})


##########################
### Vehicular Networks ###
##########################

dirName <- paste0(dataFolder, "vehicular.net/")
fileNames <- list.files(dirName)


# gstream without node set
gstream <- lapply(fileNames, function(fileName){
  gg <- read.graph(paste0(dirName, fileName), format="pajek")
  V(gg)$name <- V(gg)
  gg
})


#################################################
### Spatio-Temporal Public Transport Networks ###
#################################################

dataFileName <- "underground_stnet.csv"

dataflow <- read.csv(paste0(dataFolder, dataFileName))
data.lst <- split(dataflow, dataflow$layer1)
gstream <- lapply(data.lst, function(edges)
  graph.edgelist(as.matrix(edges[,c(1,2)])))


##########################
### Autonomous Systems ###
##########################

# AS gstream
library("igraph")
library("scatterplot3d")
#library("rgl")
library("plot3Drgl")
library("parallel")

options("mc.cores"=2) #=10)
options(scipen=999) # to deal with large node IDs

dirName <- "~/Documents/as-733/"
dirName <- "~/Downloads/as-733/"
dirName <- "~/Research/as-733/"
dirName <- "~/Downloads/as_topology_daily/"
dirName <- "~/Research/as_topology_daily/"
fileNames <- list.files(dirName)


# gstream without node set
gstream <- mclapply(fileNames, function(fileName){
  #gstream <- lapply(fileNames, function(fileName){
  edges <- read.table(paste0(dirName, fileName))
  graph.data.frame(edges, directed = T)
})
# gstream with aggregated node set
nodes <- unique(do.call(c, mclapply(fileNames, function(fileName){
  #nodes <- unique(do.call(c, lapply(fileNames, function(fileName){
  edges <- read.table(paste0(dirName, fileName))
  union(edges[,1], edges[,2])
})))
gstream <- mclapply(fileNames, function(fileName){
  #gstream <- lapply(fileNames, function(fileName){
  edges <- read.table(paste0(dirName, fileName))
  graph.data.frame(edges, directed = T, vertices = data.frame(name=nodes))
})


####


dataFolder <- "datasets/"
datasetNames <- c("primaryschool", "HS2013", "thiers_2012", "thiers_2011", "hospital", "listcontacts_2009")
ds.sep <- c("\t", " ", "\t", "\t", "\t", "\t")
timeOffset <- c(20, 20, 20, 10, 20, 20)

ds.id <- 6


dd <- read.csv(paste0(dataFolder, datasetNames[ds.id], ".csv"), sep=ds.sep[ds.id], header=F)
names(dd) <- c("time", "a", "b", "a_class", "b_class")[1:length(dd)]
md <- read.table(paste0(dataFolder, "metadata_", datasetNames[ds.id], ".txt"))
names(md) <- c("nodes", "class", "gender")[1:length(md)]

nodes <- data.frame(a=md$nodes, b=md$nodes) # id and labels match
edges <- dd[,c("a", "b")]
#edges <- data.frame(t(apply(edges,1,function(pair) if(pair[1] < pair[2]) pair else rev(pair))))
edges.dyn <- data.frame(start=dd$time-timeOffset[ds.id], end=dd$time) * 1000

dup.filter <- !duplicated(cbind(edges, edges.dyn))
#dup.filter <- !duplicated(edges.dyn)

gg <- write.gexf(nodes, edges[dup.filter,]
                 , edgeDynamic = edges.dyn[dup.filter,]
                 , tFormat = "integer"
                 , defaultedgetype = "undirected"
                 , meta = list(creator="UCY"
                               , description="A graph file writing in R using \"rgexf\""
                               ,keywords="gexf graph, R, rgexf")
                 , output = paste0(dataFolder, datasetNames[ds.id], "_output.gexf")
                 )

#gg <- read.gexf("~/Downloads/SMSdataset_graph.gexf")

