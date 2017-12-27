
#
# Author: Luis F. Chiroque
# e-mail: lf.chiroque@imdea.org
# IMDEA Networks Institute
#

library(zoo) # na.locf()
library(Rcpp) # sourceCpp()

sourceCpp("memoization_optimal.cpp")

min.err <- function(env, pos, i, m) {
  #print(paste0("pos=", pos, " i=", i, " m=", m))
  if ( pos > env$N )
    return(0)
  if ( is.na(env$e[pos,i,m]) ) {
    env$e[pos,i,m] <- env$distFun(env$metric[pos], env$metric[i]) + min.err(env, pos+1, i, m)
    if( m < env$M ) {
      alt <- min.err(env, pos+1, pos, m+1)
      if ( env$e[pos,i,m] > alt ) {
        env$e[pos,i,m] <- alt
        env$obs[pos,i,m] <- TRUE
      }
    }
  }
  return(env$e[pos,i,m])
}

opt_obs <- function(env, pos, i, m) {
  if ( pos > env$N )
    return()
  if ( env$obs[pos,i,m] ) {
    env$opt.points[pos] <- TRUE
    opt_obs(env, pos+1, pos, m+1)
  } else {
    opt_obs(env, pos+1, i, m)
  }
}

getOptimalObsPoints.env <- function(obs.ratio, metric) {
  env <- new.env(hash = FALSE)
  env$metric <- metric
  env$distFun <- function(a, b) abs(a-b)
  env$N <- length(metric)
  env$M <- floor(N * obs.ratio)
  edim <- c(env$N, env$N - 1, env$M)
  env$e <- rep(NA, prod(edim))
  dim(env$e) <- edim
  env$obs <- rep(FALSE, prod(edim))
  dim(env$obs) <- edim
  
  opt.err <- min.err(env, 2, 1, 1)
  
  env$opt.points <- rep(FALSE, env$N)
  env$opt.points[1] <- TRUE
  
  opt_obs(env, 2, 1, 1)
  which(env$opt.points) # RETURN
}

getOptimalObsPoints <- function(obs.ratio, metric) {
  min.err <- function(pos, i, m, distFun = function(a, b) abs(a-b)) {
    #print(paste0("pos=", pos, " i=", i, " m=", m))
    if ( pos > N )
      return(0)
    cur <- e[pos,i,m]
    if ( is.na(cur) ) {
      cur <- distFun(metric[pos], metric[i]) + min.err(pos+1, i, m, distFun)
      if( m < M ) {
        alt <- min.err(pos+1, pos, m+1, distFun)
        if ( cur > alt ) {
          cur <- alt
          obs[pos,i,m] <<- TRUE
        }
      }
      e[pos,i,m] <<- cur
    }
    return(cur)
  }
  opt_obs <- function(pos, i, m) {
    if ( pos > N )
      return()
    if ( obs[pos,i,m] ) {
      opt.points[pos] <<- TRUE
      opt_obs(pos+1, pos, m+1)
    } else {
      opt_obs(pos+1, i, m)
    }
  }
  
  N <- length(metric)
  M <- floor(N * obs.ratio)
  edim <- c(N, N - 1, M)
  if (prod(edim) > 1e8)
    return(1)
  e <- rep(NA, prod(edim))
  dim(e) <- edim
  obs <- rep(FALSE, prod(edim))
  dim(obs) <- edim
  
  opt.err <- min.err(2, 1, 1) #OPTIONAL: distFun()
  
  opt.points <- rep(FALSE, N)
  opt.points[1] <- TRUE
  
  opt_obs(2, 1, 1)
  which(opt.points) # RETURN
}

getOptimalObsPoints.rcpp <- function(obs.ratio, metric) {
  N <- length(metric)
  M <- floor(N * obs.ratio)
  opt.points <- computeOptimalObservations(metric, N, M, 2, 1, 1)
  which(opt.points) # RETURN
}

getOptimalObsMetric <- function(obs.ratio, metric, plot=FALSE) {
  N <- length(metric)
  ## obtain optimal points
  #x.seq <- getOptimalObsPoints(obs.ratio, metric)
  x.seq <- getOptimalObsPoints.rcpp(obs.ratio, metric)
  ## obtain optimal observed metric
  metric.obs <- rep(NA, N)
  metric.obs[x.seq] <- metric[x.seq]
  metric.obs <- na.locf(metric.obs)
  if ( plot ) {
    plot(metric, type="l")
    lines(metric.obs, type="s", lty=2)
    points(x.seq, metric[x.seq], pch=2)
  }
  metric.obs # RETURN
}

getEquidObsMetric <- function(obs.ratio, metric, plot=FALSE) {
  N <- length(metric)
  M <- round(N * obs.ratio)
  ## obtain equid points
  x.seq <- c(1, unique(round(1:(M-1) * 1/obs.ratio)) + 1)
  ## obtain optimal observed metric
  metric.obs <- rep(NA, N)
  metric.obs[x.seq] <- metric[x.seq]
  metric.obs <- na.locf(metric.obs)
  if ( plot ) {
    plot(metric, type="l")
    lines(metric.obs, type="s", lty=2)
    points(x.seq, metric[x.seq], pch=2)
  }
  metric.obs # RETURN
}

getRandObsMetric <- function(obs.ratio, metric, plot=FALSE) {
  N <- length(metric)
  M <- round(N * obs.ratio)
  ## obtain equid points
  x.seq <- c(1, sort(sample(2:N, M-1)))
  ## obtain optimal observed metric
  metric.obs <- rep(NA, N)
  metric.obs[x.seq] <- metric[x.seq]
  metric.obs <- na.locf(metric.obs)
  if ( plot ) {
    plot(metric, type="l")
    lines(metric.obs, type="s", lty=2)
    points(x.seq, metric[x.seq], pch=2)
  }
  metric.obs # RETURN
}

## SMAPE: error formula
#sum(abs(metric - metric.obs)) / sum(metric + metric.obs)
errFun <- function(metric, metric.obs) {
  sum(abs(metric - metric.obs)) / sum(metric + metric.obs)
}

