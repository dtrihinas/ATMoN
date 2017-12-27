
# R-code to compute and plot metrics

#library("rgexf")
library("igraph")
library("lubridate")
library("parallel")


# Functions
topology.fun <- function(metric, addTerm=0, multTerm=1) {
  abs(metric - addTerm) * multTerm
}

compute_sensor <- function(metric, metric.size, sensor.lim.thrs=1) {
  metric[which(is.na(metric))] <- 0
  metric.timeMean <- sapply(1:length(metric)
                            , function(i) mean(metric[1:i], na.rm=T))
  metric.timeSD <- sapply(1:length(metric)
                          , function(i) sd(metric[1:i], na.rm=T))
  metric.sensor <- topology.fun(metric, metric.timeMean, metric.size)
  metric.lim <- topology.fun(metric,  metric.timeMean + sqrt(metric.timeSD) / sensor.lim.thrs, metric.size)
  #metric.sensor <- topology.fun(metric)
  #metric.lim <- topology.fun(0, -metric.timeMean - sqrt(metric.timeSD) / sensor.lim.thrs)
  metric.alarm <- which(metric.lim < metric.sensor)
  list(metric=metric.sensor, metric.lim=metric.lim, metric.alarm=metric.alarm)
}

# dist functions
## f(x, y) => R
## one-dim; x,y numeric
abs.diff <- function(x, y) abs(x - y)
## multi-dim; x,y labels
rank.diff <- function(x, y, k.top, eq.weight=0) {
  k <- min(length(x), length(y), k.top)
  if (k==0) return(0)
  x <- names(sort(x, decreasing=T))
  y <- names(sort(y, decreasing=T))
  bool.lvl <- c("TRUE", "FALSE")
  set.diff <- table(factor(x[1:k] %in% y[1:k], levels=bool.lvl))["FALSE"] / k.top
  eq.diff <- table(factor(x[1:k]==y[1:k], levels=bool.lvl))["FALSE"] / k.top
  set.diff * (1 - eq.weight) + eq.diff * eq.weight
}
## clusters diff
clus.diff <- function(x, y) {
  # if (is.null(names(x))) {
  #   names(x) <- x
  # }
  # if (is.null(names(y))) {
  #   names(y) <- y
  # }
  names <- names(x)[names(x) %in% names(y)]
  if( length(names) > 0 ) {
    #1 - comparing.Partitions(x[names], y[names])
    1 - compare(x[names], y[names], method="nmi")
  } else { # no common vertex
    1
  }
}

multiobj.na.fill <- function(z, fill) {
  if (class(z) == "list"){
    z[which(is.na(z))] <- fill
  } else { # assume vector
    z <- na.fill(z, fill)
  }
  z
}

adam_sampling.multidim <- function(metric, sensor.obj=NULL, topology.weight=-1
                                   , k.top=5, rank.eq.weight=.5) {
  if (class(metric) == "list") {
    dist.fun <- function(x, y) rank.diff(x, y, k.top = k.top, eq.weight = rank.eq.weight)
  } else {
    dist.fun <- abs.diff 
  }
  metric <- multiobj.na.fill(metric, 0)
  metric.df <- data.frame(t=seq(0,length(metric)-1), v=I(metric))
  if (!is.null(sensor.obj)) {
    sensor.df <- data.frame(t=seq(0,length(sensor.obj$metric)-1), v=sensor.obj$metric)
    alarms <- rep(0,length(metric))
    alarms[sensor.obj$metric.alarm] <- 1
  }
  cc <- 1
  ii <- 1
  adam.sample.collection <- list()
  adam.sample.collection[[cc]] <- adam_function.generic(c(list(t=metric.df[1,1]), v=metric.df[1,2])
                                                        , dist.fun = dist.fun)
  if (!is.null(sensor.obj)) {
    adam.sample.collection.sensor <- list()
    adam.sample.collection.sensor[[cc]] <- adam_function.generic(as.list(sensor.df[1,]))
    counter.decision <- list(metric=0, topology=0, both=0)
  }
  while( ii + adam.sample.collection[[cc]]$TNext - 1 < length(metric) ){
    prev.ii <- ii
    cc <- cc + 1
    ii <- ii + adam.sample.collection[[cc-1]]$TNext
    sample <- c(list(t=metric.df[ii,1]), v=metric.df[ii,2])
    adam.obj <- adam.sample.collection[[cc-1]]
    if (!is.null(sensor.obj)) {
      sample.sensor <- as.list(sensor.df[ii,])
      adam.obj.sensor <- adam.sample.collection.sensor[[cc-1]]
      sumAlarms <- sum(alarms[(prev.ii+1):ii])
    }
    ##adam.obj$v.prev <- metric.lst[[ii-1]]$v
    adam.sample.collection[[cc]] <- adam_function.generic(sample, adam.obj, dist.fun = dist.fun) #, sumAlarms)
    if (!is.null(sensor.obj)) {
      adam.sample.collection.sensor[[cc]] <- adam_function.generic(sample.sensor, adam.obj.sensor)
      if (adam.sample.collection[[cc]]$TNext < adam.sample.collection.sensor[[cc]]$TNext) {
        counter.decision$metric <- counter.decision$metric + 1
      } else if (adam.sample.collection[[cc]]$TNext > adam.sample.collection.sensor[[cc]]$TNext) {
        counter.decision$topology <- counter.decision$topology + 1
      } else {
        counter.decision$both <- counter.decision$both + 1
      }
      if( topology.weight >= 0 ) {
        TNext <- (1-topology.weight) * adam.sample.collection[[cc]]$TNext + topology.weight * adam.sample.collection.sensor[[cc]]$TNext
      } else {
        TNext <- min(adam.sample.collection[[cc]]$TNext
                     , adam.sample.collection.sensor[[cc]]$TNext)
      }
      adam.sample.collection[[cc]]$TNext <- TNext
      adam.sample.collection.sensor[[cc]]$TNext <- TNext
    }
  }
  adam.res <- unlist(adam.sample.collection, recursive = F)
  dim(adam.res) <- c(4, length(adam.sample.collection))
  adam.sample.df <- data.frame(t(adam.res))
  names(adam.sample.df) <- names(adam.sample.collection[[1]])
  adam.sample.df
}

adam_sampling <- function(metric, sensor.obj=NULL, topology.weight=-1) {
  counter.decision <- NULL
  dist.fun <- abs.diff
  #metric <- multiobj.na.fill(metric, 0)
  metric.df <- data.frame(t=seq(0,length(metric)-1), v=I(metric))
  if (!is.null(sensor.obj)) {
    sensor.df <- data.frame(t=seq(0,length(sensor.obj$metric)-1), v=sensor.obj$metric)
    alarms <- rep(0,length(metric))
    alarms[sensor.obj$metric.alarm] <- 1
  }
  cc <- 1
  ii <- 1
  adam.sample.collection <- list()
  adam.sample.collection[[cc]] <- adam_function.generic(c(list(t=metric.df[1,1]), v=metric.df[1,2])
                                                        , dist.fun = dist.fun)
  if (!is.null(sensor.obj)) {
    adam.sample.collection.sensor <- list()
    adam.sample.collection.sensor[[cc]] <- adam_function.generic(as.list(sensor.df[1,]))
    counter.decision <- list(metric=0, topology=0, both=0)
  }
  while( ii + adam.sample.collection[[cc]]$TNext - 1 < length(metric) ){
    prev.ii <- ii
    cc <- cc + 1
    ii <- ii + adam.sample.collection[[cc-1]]$TNext
    sample <- c(list(t=metric.df[ii,1]), v=metric.df[ii,2])
    adam.obj <- adam.sample.collection[[cc-1]]
    if (!is.null(sensor.obj)) {
      sample.sensor <- as.list(sensor.df[ii,])
      adam.obj.sensor <- adam.sample.collection.sensor[[cc-1]]
      sumAlarms <- sum(alarms[(prev.ii+1):ii])
    }
    ##adam.obj$v.prev <- metric.lst[[ii-1]]$v
    adam.sample.collection[[cc]] <- adam_function.generic(sample, adam.obj, dist.fun = dist.fun) #, sumAlarms)
    if (!is.null(sensor.obj)) {
      adam.sample.collection.sensor[[cc]] <- adam_function.generic(sample.sensor, adam.obj.sensor)
      if (adam.sample.collection[[cc]]$TNext < adam.sample.collection.sensor[[cc]]$TNext) {
        counter.decision$metric <- counter.decision$metric + 1
      } else if (adam.sample.collection[[cc]]$TNext > adam.sample.collection.sensor[[cc]]$TNext) {
        counter.decision$topology <- counter.decision$topology + 1
      } else {
        counter.decision$both <- counter.decision$both + 1
      }
      if( topology.weight >= 0 ) {
        TNext <- (1-topology.weight) * adam.sample.collection[[cc]]$TNext + topology.weight * adam.sample.collection.sensor[[cc]]$TNext
      } else {
        TNext <- min(adam.sample.collection[[cc]]$TNext
                     , adam.sample.collection.sensor[[cc]]$TNext)
      }
      adam.sample.collection[[cc]]$TNext <- TNext
      adam.sample.collection.sensor[[cc]]$TNext <- TNext
    }
  }
  list(adam.sample.df=do.call(rbind.data.frame, adam.sample.collection)
       , counter.decision=counter.decision)
}

adam_sampling_result <- function(metric, sensor.obj=NULL, topology.weight=-1
                                 , multid.metrics=NULL, k.top=5, rank.eq.weight=0
                                 , cluster.mem = NULL, plot.quality=F, metricName = "") {
  metric <- multiobj.na.fill(metric, 0)
  list[adam.sample.df, counter.decision] <- adam_sampling(metric, sensor.obj, topology.weight)
  #adam.sample.obj <- adam_sampling(metric, sensor.obj, topology.weight)
  #adam.sample.df <- adam.sample.obj$adam.sample.df
  #counter.decision <- adam.sample.obj$counter.decision
  adam.xseq <- cumsum(adam.sample.df$TNext)
  adam.xseq <- adam.xseq[which(adam.xseq <= length(metric))]
  
  # Error measurement
  adam.metric <- rep(NA, length(metric))
  adam.metric[adam.xseq] <- metric[adam.xseq]
  adam.metric <- na.locf(as.data.frame(I(adam.metric)))
  
  #zoom.area <- c(0, length(metric)-1)
  #plot(metric, type="l", xlim=zoom.area)
  #lines(adam.metric, lty=2, col="red")
  
  # sMAPE (error)
  diff <- sum(abs(metric - adam.metric)) / sum(metric + adam.metric)
  compression <- 1 - length(adam.xseq) / length(metric)
  ## harmonic mean
  # 1/mean(1/c(1-error/100, 1-(length(adam.xseq) / nrow(metric) - 1/time.thrs$TMax)))
  efficiency <- 1 / mean(1 / c(1-diff, compression))
  result <- list(error=diff, compression=compression, efficiency=efficiency)
  if (!is.null(sensor.obj) && topology.weight < 0) {
    result <- c(result, list(decision=counter.decision))
  }
  
  # extra multi-dim metrics
  if (!is.null(multid.metrics)) {
    qualities <- sapply(names(multid.metrics), function(multid.metric.name) {
      multid.metric <- multid.metrics[[multid.metric.name]]
      quality <- sapply(2:length(adam.xseq)-1, function(i) {
        1 - rank.diff(multid.metric[[adam.xseq[i]]], multid.metric[[adam.xseq[i+1]-1]]
                      , k.top = k.top, eq.weight = rank.eq.weight)
      })
      if (plot.quality) {
        suffix <- ifelse(is.null(sensor.obj), "", " + tplgy")
        par(mfrow=c(1,2))
        plot(quality, type="l", ylim=c(0,1)
             , main=multid.metric.name
             , xlab=paste0(metricName, " sampling", suffix))
        plot(sapply(1:length(quality), function(k) mean(quality[1:k])), type="l", ylim=c(0,1)
             , main=multid.metric.name
             , ylab="quality avg"
             , xlab=paste0(metricName, " sampling", suffix))
        par(mfrow=c(1,1))
      }
      mean(quality)
    })
    result <- c(result, qualities)
  }
  if (!is.null(cluster.mem)) {
    cluster.perf <- sapply(2:length(adam.xseq)-1, function(i) {
      1 - clus.diff(cluster.mem[[adam.xseq[i]]], cluster.mem[[adam.xseq[i+1]-1]])
    })
    result <- c(result, cluster.perf=mean(cluster.perf))
  }
  # RETURN
  result
}

benchmark_sampling_result <- function(metric, compression.ratio=.5) {
  # average
  avgM <- mean(metric, na.rm=T)
  avg.err <- sum(abs(metric - avgM), na.rm=T) / sum(metric + avgM, na.rm = T)
  
  # equidist
  avgT <- 1 / (1-compression.ratio)
  equid.xseq <- cumsum(sample(floor(avgT):ceiling(avgT)
                              , round(length(metric)*(1-compression.ratio))
                              , replace=T, prob=c(ceiling(avgT) - avgT, avgT - floor(avgT))))
  equid.xseq <- equid.xseq[which(equid.xseq <= length(metric))]
  equid.metric <- rep(NA, length(metric))
  equid.metric[c(1,equid.xseq)] <- metric[c(1,equid.xseq)]
  equid.metric <- na.locf(equid.metric)
  equid.err <- sum(abs(metric - equid.metric), na.rm=T) / sum(metric + equid.metric, na.rm=T)
  
  # random
  M <- 50 # TODO: add $M as a funciton parameter
  err.M <- sapply(1:M, function(XXX) {
    rnd.xseq <- cumsum(sample(seq(time.thrs$TMin, time.thrs$TMax)
                              , length(metric), replace=T))
    rnd.xseq <- rnd.xseq[which(rnd.xseq <= length(metric))]
    rnd.metric <- rep(NA, length(metric))
    rnd.metric[c(1,rnd.xseq)] <- metric[c(1,rnd.xseq)]
    rnd.metric <- na.locf(rnd.metric)
    sum(abs(metric - rnd.metric), na.rm=T) / sum(metric + rnd.metric, na.rm=T)
  })
  rnd.err <- mean(err.M) # average error
  rnd.se <- sd(err.M) / sqrt(M) # standard error
  # RETURN
  list(avg=list(error=avg.err)
       , equid=list(error=equid.err)
       , rnd=list(error=rnd.err, sd=rnd.se))
}

compute_adam_params <- function(metricName, dataset, topology.obj
                                , topology.weight=-1, k.top=5, rank.eq.weight=.5
                                , plots = FALSE) {
  multid.metrics <- list(pageRank.quality=dataset$metrics$pageRank
                         , outDegree.quality=dataset$metrics$outDegree
                         , betweenness.quality=dataset$metrics$betweenness
                         , l.rank.quality=dataset$metrics$l.rank)
  # AdaM
  adam.result <- adam_sampling_result(dataset$metrics[[metricName]]
                                      , multid.metrics = multid.metrics, k.top = k.top, rank.eq.weight = rank.eq.weight
                                      , cluster.mem = dataset$metrics$cluster.mem
                                      , plot.quality = plots
                                      , metricName = metricName)
  # AdaM + topology
  adam.topology.result <- adam_sampling_result(dataset$metrics[[metricName]]
                                               , topology.obj, topology.weight
                                               , multid.metrics = multid.metrics, k.top = k.top, rank.eq.weight = rank.eq.weight
                                               , cluster.mem = dataset$metrics$cluster.mem
                                               , plot.quality = plots
                                               , metricName = metricName)
  # benchmark
  bench.results <- benchmark_sampling_result(dataset$metrics[[metricName]]
                                             , adam.result$compression)
  # build results vector
  result.vector <- unlist(adam.topology.result)
  names(result.vector) <- paste0(names(result.vector), ".topology")
  c(unlist(adam.result), result.vector, unlist(bench.results))
}

compute_performance <- function(target.topology, dataset
                                , target.metrics = c("diameter", "nofComponents", "deg.dist.pow"
                                                     , "gcRatio", "esize", "max.clique"
                                                     , "avgPageRank", "avgOutDegree")
                                , topology.thrs = 1, topology.weight=-1
                                , k.top=5, rank.eq.weight=0
                                , plots = FALSE) {
  topology.obj <- compute_sensor(dataset$metrics[[target.topology]]
                                 , dataset$metrics$vsize
                                 , topology.thrs)
  results <- t(sapply(target.metrics, compute_adam_params
                      , dataset, topology.obj, topology.weight, k.top, rank.eq.weight, plots))
}

# one-dim
compute_dataset_performance <- function(dataset, target.topologies = c("assortativity", "transitivity")
                                        , target.metrics = c("diameter", "nofComponents", "deg.dist.pow"
                                                             , "gcRatio", "esize", "max.clique"
                                                             , "avgPageRank", "avgOutDegree")
                                        , topology.thrs = 1, topology.weight=-1
                                        , k.top=10, rank.eq.weight=0
                                        , plots = FALSE) {
  results <- lapply(target.topologies, compute_performance
                    , dataset, target.metrics = target.metrics
                    , topology.thrs = topology.thrs, topology.weight = topology.weight
                    , k.top = k.top, rank.eq.weight = rank.eq.weight
                    , plots = plots)
  names(results) <- target.topologies
  results
}

# multi-dim (quality)
compute_dataset_qPerformance <- function(dataset, target.topologies = c("assortativity", "transitivity")
                                         , topology.thrs = 1, topology.weight=-1
                                         , k.top=10, rank.eq.weight=.5) {
  results <- lapply(target.topologies, compute_performance
                    , dataset, target.metrics = c("pageRank", "outDegree"
                                                  , "betweenness", "l.rank")
                    , topology.thrs = topology.thrs, topology.weight = topology.weight
                    , k.top=k.top, rank.eq.weight=rank.eq.weight)
  names(results) <- target.topologies
  results
}



tops <- c(5,10,15,20)
performance <- lapply(tops, function(top.k) {
  #lapply(gstream.collection, compute_dataset_performance, k.top=top.k, plots=T)
  lapply(perf.config, function(config) {
    print(config$ds.name)
    compute_dataset_performance(gstream.collection[[config$ds.name]]
                                #, target.topologies = c("vsize", "transitivity", "assortativity", "dissimilarity")
                                , target.topologies = c("transitivity", "assortativity", "dissimilarity")
                                , target.metrics = config$target.metrics
                                , k.top = top.k, plots=T)
  })
})
names(performance) <- paste0("top", tops)

perf.objFileName <- "performance.Robj"
if (!file.exists(perf.objFileName)) {
  save(performance, file=perf.objFileName)
}

target.multimetrics = c("pageRank", "outDegree"
                        , "betweenness", "l.rank"
                        , "cluster"
)
sapply(names(gstream.collection), function(ds.name) {
  sapply(target.multimetrics, function(multid.metricName) {
    sapply(perf.config[[ds.name]]$target.metrics, function(metricName){
      values <- sapply(performance, function(perf.top) {
        affix <- ifelse(multid.metricName == "cluster", ".perf", ".quality")
        colNames <- paste0(multid.metricName, affix, c("", ".topology"))
        #perf.top[[ds.name]]$assortativity[metricName,colNames]
        perf.top[[ds.name]]$vsize[metricName,colNames]
      })
      plot(tops, values[1,], type="l", lty=2, ylim=c(.25,1)
           , xlab="top", ylab="quality"
           , main=paste0(ds.name, " ", multid.metricName, " @", metricName))
      lines(tops, values[2,])
    })
  })
  invisible()
})

performance.top5 <- lapply(gstream.collection, compute_dataset_performance, k.top=5, plots=T)
performance.top10 <- lapply(gstream.collection, compute_dataset_performance, k.top=10)
performance.top15 <- lapply(gstream.collection, compute_dataset_performance, k.top=15)
performance.top20 <- lapply(gstream.collection, compute_dataset_performance, k.top=20)

performance.topology.w1 <- lapply(gstream.collection, compute_dataset_performance
                                  , topology.weight=1)
performance.topology.w08 <- lapply(gstream.collection, compute_dataset_performance
                                   , topology.weight=.8)
performance.topology.w06 <- lapply(gstream.collection, compute_dataset_performance
                                   , topology.weight=.6)
performance.topology.w04 <- lapply(gstream.collection, compute_dataset_performance
                                   , topology.weight=.4)
performance.topology.w02 <- lapply(gstream.collection, compute_dataset_performance
                                   , topology.weight=.2)
performance.topology.w0 <- lapply(gstream.collection, compute_dataset_performance
                                  , topology.weight=0)


## quality
qPerformance.top5.rankw0 <- lapply(list("sg.24hrs.09", "sms.mit.24hrs.09")
                                   , function(ds.name) compute_dataset_qPerformance(gstream.collection[[ds.name]]
                                                                                    , k.top = 5, rank.eq.weight = 0))
qPerformance.top5.rankw1 <- lapply(list("sg.24hrs.09", "sms.mit.24hrs.09")
                                   , function(ds.name) compute_dataset_qPerformance(gstream.collection[[ds.name]]
                                                                                    , k.top = 5, rank.eq.weight = 1))
qPerformance.top5 <- lapply(list("sg.24hrs.09", "sms.mit.24hrs.09")
                            , function(ds.name) compute_dataset_qPerformance(gstream.collection[[ds.name]]
                                                                             , k.top = 5, rank.eq.weight = .5))
qPerformance.top5 <- lapply(list("sg.24hrs.09", "sms.mit.24hrs.09")
                            , function(ds.name) compute_dataset_qPerformance(gstream.collection[[ds.name]]
                                                                             , k.top = 5, rank.eq.weight = .5
                                                                             , topology.weight = 1))

library("xlsx")

resultFileName <- "adam-performance.xlsx"

df2xlsx <- function(df.name, performance, resultFileName="result.xlsx") {
  target.topologies <- names(performance[[df.name]])
  for( topologyName in target.topologies ) {
    write.xlsx(performance[[df.name]][[topologyName]]
               , resultFileName, sheetName=paste0(df.name, ".", topologyName)
               , append=T)
  }
}

trash <- sapply(names(performance), df2xlsx, performance, resultFileName)



get_diameters <- function(gg){
  dists <- distances(gg, mode="out")
  dists[dists==Inf] <- vcount(gg)+1
  apply(dists, 1, max) # regardless 1 or 2 if unidirected
}
percentile_diameter <- function(gg, percentile=.9) {
  diams <- get_diameters(gg)
  distr <- cumsum(table(diams))/sum(table(diams))
  ffun <- splinefun(distr, names(distr), method="hyman")
  ffun(quantile(distr, percentile))
}
metric <- sapply(gstream, percentile_diameter)
#
metric <- sapply(gstream, closeness)
#
metric.diameter <- sapply(gstream, diameter)
metric.pr <- sapply(gstream, function(gg) page.rank(gg)$vector)
metric.outDeg <- sapply(gstream, function(gg) degree(gg, mode='out'))
#metric <- sapply(gstream, diameter)
metric.noComp <- sapply(gstream, function(gg) components(gg)$no)
#
metric.avgPR <- sapply(metric.pr, mean)
metric.avgOutDeg <- sapply(metric.outDeg, mean)
metric.size <- sapply(gstream, vcount)
metric.gcRatio <- sapply(gstream, function(gg) {
  cc <- components(gg)
  ifelse(cc$no > 0
         , table(cc$membership)[as.character(which.max(cc$csize))] / vcount(gg)
         , 0)
})

# Topology Metrics
metric.assortativity <- sapply(gstream, assortativity_degree, directed=T)
metric.transitivity <- sapply(gstream, transitivity, type="global")


lim <- 90 #ncol(metric)
ss <- sort(sample(1:nrow(metric), 250)) #1:nrow(metric)
#persp3Drgl(
persp(1:length(ss)
      , 1:lim
      #, metric[order(metric[,1], decreasing=T),1:lim]
      #, metric[ss,1:lim]
      , metric[ss,1:lim + ncol(metric) - lim]
      , box=F, theta=15, expand=.5)

#scatterplot3d(kk[1:(nrow(metric)*2),], box=F, cex.symbols=.25)

scatterplot3d(rep(1:lim, each=nrow(metric))
              , rep(1:nrow(metric), lim)
              , as.vector(metric)[1:(nrow(metric)*lim)]
              , box=F, cex.symbols=.25)
segments(metric[1,1], 1, metric[7716,1], 1)
##

x.seq <- seq(0, 999999, by=5)
metric <- sapply(x.seq, function(i, metric_fun, timeWindow){
  thrs <- (min(edges.dyn$start) + i*timeUnit)
  filter <- which(edges.dyn$start >= thrs & edges.dyn$start < (thrs + timeWindow))
  gg <- simplify(graph.data.frame(unique(edges[filter,]), directed=F, vertices=nodes))
  metric_fun(gg)
}, diameter, timeOffset[ds.id])

plot(x.seq, metric, type="l", ylim=c(0,max(metric)))

# Export metric
metricName <- "degree"
graphName <- "as-733"
## multi-dim
df.metric <- data.frame(ts=0:(ncol(metric)-1)
                        , metric=colMeans(metric))
# uni-dim
df.metric <- data.frame(ts=0:(length(metric)-1)
                        , metric=metric)
write.table(df.metric, paste0(graphName, "_", metricName, ".csv"), row.names=F, col.names=F, sep=",")


## plot time correlations

ws <- ncol(metric) #20 # windows size
plot(NULL, type='n', ylim=c(-1,1), xlim=c(1,ws), ylab="cor")
sapply(3:ncol(metric), function(i, ws=20) {
  corr <- cor(metric[,max(1,(i-1-ws)):(i-1)], metric[,i])
  lines(rev(corr[,1]), lty=2)
  invisible()
}, ws)



### plots

plotsFolderName <- "plots/"
gstreamName <- "" #"sg09", "as-733"

plot_sensor <- function(metric.sensor, sensorName
                        , plotsFolderName="", gstreamName=""
                        , zoom.area=c(1,length(metric.sensor))
                        , toFile=F) {
  if( toFile )
    png(paste0(plotsFolderName, sensorName, "_onTime-", gstreamName, ".png")
        , width = 800, height = 1200)
  par(mfrow=c(2,1))
    plot(metric.sensor, type="l"
         , main=paste0(sensorName, " weighted mean distance"))
    plot(metric.sensor, type="l", main="zoomed in", xlim=zoom.area)
  par(mfrow=c(1,1))
  if( toFile )
    dev.off()
}

plot_metric_sensor <- function(metric.sensor, metric.test, sensorName, testMetricName
                               , plotsFolderName="", gstreamName=""
                               , zoom.area=c(1, length(metric.sensor))
                               , toFile=F){
  metric.cor <- cor(metric.sensor, metric.test, use="complete.obs")
  if( toFile ) 
    png(paste0(plotsFolderName, testMetricName, "-", sensorName, "_cor-", gstreamName, ".png")
        , width = 800, height = 1200)
  par(mfrow=c(2,1))
    plot(metric.test, type="l", xlim=zoom.area
         , main=paste0(testMetricName, " cor=", signif(metric.cor,4)))
    plot(metric.sensor, type="l", xlim=zoom.area, main=sensorName)
  par(mfrow=c(1,1))
  if( toFile )
    dev.off()
}


if( ! file.exists(plotsFolderName) )
  dir.create(plotsFolderName)

metric.assortativity[which(is.na(metric.assortativity))] <- 0
metric.transitivity[which(is.na(metric.transitivity))] <- 0
metric.assort.timeMean <- sapply(1:length(metric.assortativity)
                                 , function(i) mean(metric.assortativity[1:i], na.rm=T))
metric.assort.timeSD <- sapply(1:length(metric.assortativity)
                                 , function(i) sd(metric.assortativity[1:i], na.rm=T))
metric.transi.timeMean <- sapply(1:length(metric.transitivity)
                                 , function(i) mean(metric.transitivity[1:i], na.rm=T))
metric.transi.timeSD <- sapply(1:length(metric.transitivity)
                                 , function(i) sd(metric.transitivity[1:i], na.rm=T))
#topology mean evolution
png(paste0(plotsFolderName, "topology_evolution-", gstreamName, ".png")
    , width = 800, height = 1200)
par(mfrow=c(2,1))
  plot(metric.assort.timeMean, type="l")
  plot(metric.transi.timeMean, type="l")
par(mfrow=c(1,1))
dev.off()

zoom.area <- c(300,600) # 1:length(gstream)
sensor.lim.thrs <- 1
# Assortativity
sensorName <- "assortativity"
metric.sensor <- abs(metric.assortativity - metric.assort.timeMean) * metric.size
metric.sensor.lim <- abs(metric.assortativity - metric.assort.timeMean - sqrt(metric.assort.timeSD) / sensor.lim.thrs) * metric.size
metric.sensor.alarm <- which(metric.sensor.lim < metric.sensor)
## update metric obj
metric$assort.alarm <- rep(0,length(gstream))
metric$assort.alarm[metric.sensor.alarm] <- 1
### generic
## assortativity weighted mean distance
plot_sensor(metric.sensor, sensorName
            , plotsFolderName = plotsFolderName, gstreamName = gstreamName
            , zoom.area = zoom.area, toFile = T)

# Transitivty
sensorName <- "transitivity"
metric.sensor <- abs(metric.transitivity - metric.transi.timeMean) * metric.size
metric.sensor.lim <- abs(metric.transitivity - metric.transi.timeMean - sqrt(metric.transi.timeSD) / sensor.lim.thrs) * metric.size
metric.sensor.alarm <- which(metric.sensor.lim < metric.sensor)
## update metric obj
metric$transi.alarm <- rep(0,length(gstream))
metric$transi.alarm[metric.sensor.alarm] <- 1
### generic
## assortativity weighted mean distance
plot_sensor(metric.sensor, sensorName
            , plotsFolderName = plotsFolderName, gstreamName = gstreamName
            , zoom.area = zoom.area, toFile = T)


## checking no. components
testMetricName <- "no. components"
metric.test <- metric.noComp
plot_metric_sensor(metric.sensor, metric.test, sensorName, testMetricName
                   , plotsFolderName = plotsFolderName, gstreamName = gstreamName
                   , zoom.area = zoom.area, toFile=T)

## checking diameter
testMetricName <- "diameter"
metric.test <- metric.diameter
plot_metric_sensor(metric.sensor, metric.test, sensorName, testMetricName
                   , plotsFolderName = plotsFolderName, gstreamName = gstreamName
                   , zoom.area = zoom.area, toFile=T)

## checking avg out degree
testMetricName <- "avgOutDegree"
metric.test <- metric.avgOutDeg
plot_metric_sensor(metric.sensor, metric.test, sensorName, testMetricName
                   , plotsFolderName = plotsFolderName, gstreamName = gstreamName
                   , zoom.area = zoom.area, toFile=T)

## checking avg pageRank
testMetricName <- "avgPR"
metric.test <- metric.avgPR
plot_metric_sensor(metric.sensor, metric.test, sensorName, testMetricName
                   , plotsFolderName = plotsFolderName, gstreamName = gstreamName
                   , zoom.area = zoom.area, toFile=T)


### optimal ###

obs.ratios <- seq(.1, .9, by=.2)
err <- sapply(obs.ratios, function(obs.ratio, metric) {
  metric.obs <- getOptimalObsMetric(obs.ratio, metric)
  errFun(metric, metric.obs)
}, metric)

plot(obs.ratios, err, type="l", lty=2)

sapply(names(perf.top), function(ds.name) {
  performance.ds <- perf.top[[ds.name]]
  sapply(rownames(performance.ds$transitivity), function(metric.name) {
    #performance.ds
    #metric.name
    print(paste0(ds.name, " [", metric.name, "]"))
    metric <- gstream.collection[[ds.name]]$metrics[[metric.name]]
    metric[which(is.na(metric))] <- 0
    obs.ratios <- c(.05, seq(.1, .9, by=.2))
    err <- sapply(obs.ratios, function(obs.ratio, metric) {
      metric.obs <- getOptimalObsMetric(obs.ratio, metric)
      errFun(metric, metric.obs)
    }, metric)
    err.equid <- sapply(obs.ratios, function(obs.ratio, metric) {
      metric.obs <- getEquidObsMetric(obs.ratio, metric)
      errFun(metric, metric.obs)
    }, metric)
    err.rand <- sapply(obs.ratios, function(obs.ratio, metric) {
      mean(sapply(1:100, function(s) {
        metric.obs <- getRandObsMetric(obs.ratio, metric)
        errFun(metric, metric.obs)
      }))
    }, metric)
    err.avg <- errFun(metric, mean(metric))
    
    perf.points <- rbind(temporal=performance.ds[[1]][metric.name, c("error", "compression")]
                         , t(sapply(names(performance.ds), function(topology.name) {
                           perf.obj <- performance.ds[[topology.name]]
                           perf.obj[metric.name, c("error.topology", "compression.topology")]
                         })))
    
    plot(obs.ratios, err, type="l", lty=2
         , ylim=range(c(err, err.equid, err.rand, err.avg, perf.points[,1]))
         , main=paste0(ds.name, " [", metric.name, "]"))
    lines(obs.ratios, err.equid, lty=3)
    lines(obs.ratios, err.rand, lty=4)
    abline(h=err.avg, lty=5)
    sapply(seq_len(nrow(perf.points)), function(i) {
      points(1-perf.points[i,2], perf.points[i,1], pch=i)
    })
    legend("topright", pch=seq_len(nrow(perf.points))
           , legend=rownames(perf.points)
           , text.width = .15)
    invisible()
  })
})
