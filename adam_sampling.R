
#
# Author: Luis F. Chiroque
# e-mail: lf.chiroque@imdea.org
# IMDEA Networks Institute
#

library("zoo") # na.locf()
library("devtools") # 'list[...]<-' assign operator
source_url("https://raw.githubusercontent.com/ggrothendieck/gsubfn/master/R/list.R")
library("clusterSim") # comparing.Partitions
source("adamImpl.R") # time.thrs := [Tmin, TMax]
source("optimal_adaption.R")


# performance configuration
perf.config <- list(sg.09 = list(ds.name="sg.09"
                                 , target.metrics=c("diameter", "nofComponents", "gcRatio", "deg.dist.pow"))
                    , sg.24hrs.09 = list(ds.name="sg.24hrs.09"
                                         , target.metrics=c("diameter", "nofComponents", "gcRatio", "deg.dist.pow"))
                    , thiersHS.12 = list(ds.name="thiersHS.12"
                                         , target.metrics=c("esize", "max.clique", "avgOutDegree"))
                    , reality.mit = list(ds.name="reality.mit"
                                         , target.metrics=c("esize", "deg.dist.pow"))
                    , contact.mit.05 = list(ds.name="contact.mit.05"
                                            , target.metrics=c("nofComponents", "gcRatio"))
                    , sms.mit.09 = list(ds.name="sms.mit.09"
                                        , target.metrics=c("esize"))
                    , sms.mit.24hrs.09 = list(ds.name="sms.mit.24hrs.09"
                                              , target.metrics=c("esize"))
                    #, underground = list(ds.name="underground", target.metrics=c(""))
                    , vehicular.net = list(ds.name="vehicular.net"
                                           , target.metrics=c("diameter", "nofComponents", "gcRatio"))
                    )

tplgy.sensor.lst <- list(assortativity="assortativity"
                         , transitivity="transitivity"
                         , dissimilarity="dissimilarity")

#for each dataset
##for each topology sensor
### for each metric

## get dynamic sampling results; observed metrics
# balancerFun = function(TNextA, TNextB, wei=.5) {
#   TNextA * wei + TNextB * (1 - wei)
# }
balancerFun = min

results <- lapply(perf.config, function(config) {
  lapply(tplgy.sensor.lst, function(tplgy.name) {
    sapply(config$target.metrics, function(metric.name) {
      print(paste0("dataset=", config$ds.name, " metric=", metric.name, " tplgy=", tplgy.name))
      # set up metric and topology metric
      metric <- gstream.collection[[config$ds.name]]$metrics[[metric.name]]
      metric[is.na(metric)] <- 0
      tplgy.metric <- gstream.collection[[config$ds.name]]$metrics[[tplgy.name]]
      tplgy.metric[is.na(tplgy.metric)] <- 0
      # get dyn monitoring result
      df.sampling <- dynSampling(metric, tplgy.metric, balancerFun)
      # prepare observed metric
      ss <- rep(NA, length(metric))
      x.seq <- cumsum(df.sampling$TNext)
      x.seq <- x.seq[which(x.seq <= length(metric))]
      ss[x.seq] <- metric[x.seq]
      # return: list(obs.metric, obs.ratio)
      list(obs.metric=na.locf(ss)
           , obs.ratio=length(df.sampling$TNext) / length(metric))
    }, simplify=FALSE, USE.NAMES=TRUE)
  })
})


obs.ratios <- seq(.1, .9, by=.2)

## get optimal err. vs obs. ratio
## DO NOT RUN for all config's. takes O^3 space on length(metric)
## run on metrics which length(metric) ~= 300
res.optim <- lapply(perf.config, function(config) {
  sapply(config$target.metrics, function(metric.name) {
    metric <- gstream.collection[[config$ds.name]]$metrics[[metric.name]]
    
    if( length(metric) < 500 ) {
      ll <- lapply(obs.ratios, getOptimalObsMetric, metric)
      #plot(obs.ratios, sapply(ll, errFun, metric), type="l")
      #sapply(ll, errFun, metric)
    } else {
      NULL
    }
  }, simplify=FALSE)
})


## get random err. vs obs. ratio
res.rand <- lapply(perf.config, function(config) {
  sapply(config$target.metrics, function(metric.name) {
    metric <- gstream.collection[[config$ds.name]]$metrics[[metric.name]]
    
    ll <- lapply(obs.ratios, getRandObsMetric, metric)
    #plot(obs.ratios, sapply(ll, errFun, metric), type="l")
    #sapply(ll, errFun, metric)
  }, simplify=FALSE)
})


## get equid err. vs obs. ratio
res.equid <- lapply(perf.config, function(config) {
  sapply(config$target.metrics, function(metric.name) {
    metric <- gstream.collection[[config$ds.name]]$metrics[[metric.name]]
    
    ll <- lapply(obs.ratios, getEquidObsMetric, metric)
    #plot(obs.ratios, sapply(ll, errFun, metric), type="l")
    #sapply(ll, errFun, metric)
  }, simplify=FALSE)
})


## TODO: add code below to paper_plots.R
sapply(perf.config, function(config) {
  sapply(config$target.metrics, function(metric.name) {
    # prepare metric
    metric <- gstream.collection[[config$ds.name]]$metrics[[metric.name]]
    # prepare error for optimal, random and equidistant
    err.optim <- sapply(res.optim[[config$ds.name]][[metric.name]]
                  , errFun, metric)
    err.rand <- sapply(res.rand[[config$ds.name]][[metric.name]]
                  , errFun, metric)
    err.equid <- sapply(res.equid[[config$ds.name]][[metric.name]]
                       , errFun, metric)
    # plot error lines vs. obs.ratio (1 - compression)
    plot(NA, type='n', xlim=c(0,1), ylim=c(0,.35)
         , main=paste0("dyn monitoring [", config$ds.name, "::", metric.name, "]")
         , xlab="obs. ratio", ylab="err")
    if( length(err.optim) > 0)
      lines(obs.ratios, err.optim)
    lines(obs.ratios, err.rand, lty=2)
    lines(obs.ratios, err.equid, lty=3)
    # plot our result
    ## compute error for our results
    err <- sapply(tplgy.sensor.lst, function(tplgy.name) {
      errFun(results[[config$ds.name]][[tplgy.name]][[metric.name]]$obs.metric
             , metric)
    }, simplify=FALSE)
    ## plot our error
    points(results[[config$ds.name]]$assortativity[[metric.name]]$obs.ratio
           , err$assortativity
           , pch=2)
    points(results[[config$ds.name]]$transitivity[[metric.name]]$obs.ratio
           , err$transitivity
           , pch=3)
    points(results[[config$ds.name]]$dissimilarity[[metric.name]]$obs.ratio
           , err$dissimilarity
           , pch=4)
    invisible()
  })
})
