
#
# Author: Luis F. Chiroque
# e-mail: lf.chiroque@imdea.org
# IMDEA Networks Institute
#
# AdaM - R's implementation
# source: http://dx.doi.org/10.1109/BigData.2015.7363816
#


# Functions #
abs.diff <- function(x, y) abs(x - y)

probDistro <- function(dist, dist.hat, sd.hat) {
  z <- (dist - dist.hat) / sd.hat
  if( is.na(z) ) z <- 0 # line added for solving 0/0
  1 / sqrt(2 * pi) * exp(-z^2/2)
}

PEWMA <- function(prob, dist, adam.obj, alpha, beta) {
  a <- alpha * (1 - beta * prob)
  s1 <- a * adam.obj$mean + (1 - a) * dist
  s2 <- a * adam.obj$sd + (1 - a) * dist^2
  list(mean=s1, sd=sqrt(abs(s2 - s1^2))) # TODO: review
}

calcConfidence <- function(sd, sd.hat) {
  abs.diff <- abs(sd.hat - sd)
  1 - abs.diff / ifelse(abs.diff==0, 1, sd)
}

# sample := [t, v]
# adam.obj := [mean, sd, TNext, prev.v]
adam_function.original <- function(sample, adam.obj=list()
                                  , alpha=.4, beta=1, gamma=.9, lambda=1
                                  , dist.fun=abs.diff){
  if( sample$t > 0 ){
    dist <- dist.fun(sample$v, adam.obj$v.prev)
    prob <- probDistro(dist, adam.obj$mean, adam.obj$sd)
    new.adam.obj <- PEWMA(prob, dist, adam.obj, alpha, beta)
    conf <- calcConfidence(new.adam.obj$sd, adam.obj$sd)
    
    if( conf >= (1 - gamma) ){
      TNext <- adam.obj$TNext + ceiling(lambda * (1 + (conf - gamma) / conf))
      new.adam.obj$TNext <- min(max(TNext, time.thrs$TMin), time.thrs$TMax)
      #added conf to adam object
      new.adam.obj$conf <- conf
    }else{
      new.adam.obj$TNext <- time.thrs$TMin
      #added conf to adam object
      new.adam.obj$conf <- -1
    }
  }else{
    new.adam.obj <- list(mean=ifelse(length(sample$v)>1, 1, sample$v)
                         , sd=0
                         , TNext=time.thrs$TMin, conf = 0)
  }
  new.adam.obj$v.prev <- sample$v
  # RETURN
  new.adam.obj
}

adam_function.simple <- function(adam.obj=list(mean=0, sd=0, TNext=time.thrs$TMin, v.prev=0), sample
                                   , alpha=.4, beta=1, gamma=.9, lambda=1
                                   , dist.fun=abs.diff){
  dist <- dist.fun(sample, adam.obj$v.prev)
  prob <- probDistro(dist, adam.obj$mean, adam.obj$sd)
  new.adam.obj <- PEWMA(prob, dist, adam.obj, alpha, beta)
  conf <- calcConfidence(new.adam.obj$sd, adam.obj$sd)
    
  if( conf >= (1 - gamma) ){
    TNext <- adam.obj$TNext + ceiling(lambda * (1 + (conf - gamma) / conf))
    new.adam.obj$TNext <- min(max(TNext, time.thrs$TMin), time.thrs$TMax)
    #added conf to adam object
    new.adam.obj$conf <- conf
  }else{
    new.adam.obj$TNext <- time.thrs$TMin
    #added conf to adam object
    new.adam.obj$conf <- -1
  }
  new.adam.obj$v.prev <- sample
  new.adam.obj # RETURN
}

# extension which supports (topology) alarms and configuration
adam_function.generic <- function(sample, adam.obj=list(), sumAlarms=0
                                  , alpha=.4, beta=1, gamma=.9, lambda=1
                                  , adam.conf=list()
                                  , dist.fun=abs.diff){
  if( sample$t > 0 ){
    dist <- dist.fun(sample$v, adam.obj$v.prev)
    prob <- probDistro(dist, adam.obj$mean, adam.obj$sd)
    new.adam.obj <- PEWMA(prob, dist, adam.obj, alpha, beta)
    conf <- calcConfidence(new.adam.obj$sd, adam.obj$sd)
    
    if( !is.null(adam.conf$confIncr) && adam.conf$confIncr && sumAlarms == 0 )
      conf <- conf*2
    
    if( conf >= (1 - gamma) ){
      new.adam.obj$TNext <- min(max(adam.obj$TNext + ceiling(lambda * (1 + (conf - gamma) / conf)), time.thrs$TMin)
                                , time.thrs$TMax)
      if( !is.null(adam.conf$Tdecr) && sumAlarms > 0 )
        if (adam.conf$Tdecr == "proportional") {
          new.adam.obj$TNext <- max(floor(new.adam.obj$TNext * (1 - sumAlarms / adam.obj$TNext))
                                    , time.thrs$TMin)
        } else if (adam.conf$Tdecr == "multiplicative") {
          new.adam.obj$TNext <- floor(new.adam.obj$TNext / 2)
        } else if (adam.conf$Tdecr == "hard") {
          new.adam.obj$TNext <- time.thrs$TMin
        }
      #added conf to adam object
      new.adam.obj$conf <- conf
    }else{
      new.adam.obj$TNext <- time.thrs$TMin
      #added conf to adam object
      new.adam.obj$conf <- -1
    }
  }else{
    new.adam.obj <- list(mean=ifelse(length(sample$v)>1, 1, sample$v)
                         , sd=0
                         , TNext=time.thrs$TMin, conf = 0)
  }
  new.adam.obj$v.prev <- sample$v
  new.adam.obj
}

# metric:= <vector>   the full computed metric
adam_sampling.standalone <- function(metric) {
  dist.fun <- abs.diff
  #metric <- multiobj.na.fill(metric, 0)
  t.seq <- seq(0, length(metric)-1)
  metric.lst <- split(data.frame(t=t.seq, v=metric), t.seq)
  #metric.df <- data.frame(t=seq(0,length(metric)-1), v=metric)
  cc <- 1
  ii <- 1
  adam.sample.collection <- list()
  sample <- metric.lst[[1]]
  adam.sample.collection[[cc]] <- adam_function.original(sample
                                                        , dist.fun = dist.fun)
  while( ii + adam.sample.collection[[cc]]$TNext - 1 < length(metric) ){
    prev.ii <- ii
    cc <- cc + 1
    ii <- ii + adam.sample.collection[[cc-1]]$TNext
    sample <- metric.lst[[ii]]
    adam.obj <- adam.sample.collection[[cc-1]]
    adam.sample.collection[[cc]] <- adam_function.generic(sample
                                                          , adam.obj
                                                          , dist.fun = dist.fun)
  }
  # RETURN
  do.call(rbind.data.frame, adam.sample.collection)
}

# USAGE EXAMPLE #
# UNCOMMENT FOR RUNNING EXAMPLE #
# N <- 100
# metric <- sin(seq(0,10*pi, length=N)) + rnorm(N, 0,.25)
# metric.adam.df <- adam_sampling.standalone(metric)
# # plot
# plot(metric, type="l")
# x.seq <- cumsum(metric.adam.df$TNext)
# lines(x.seq, metric[x.seq], lty=3) # interpolation
# lines(x.seq, metric[x.seq], lty=2, type="s") # actual values for measuring error


dynSampling <- function(metricA, metricB, balancerFun=min) {
  # assert(length(metricA) == length(metricB))
  dist.fun <- abs.diff
  t.seq <- seq(0, length(metricA)-1)
  metricA.lst <- split(data.frame(t=t.seq, v=metricA), t.seq)
  metricB.lst <- split(data.frame(t=t.seq, v=metricB), t.seq)
  cc <- 1
  ii <- 1
  adam.sample.collectionA <- list()
  adam.sample.collectionB <- list()
  sampleA <- metricA.lst[[1]]
  sampleB <- metricB.lst[[1]]
  adam.sample.collectionA[[cc]] <- adam_function.original(sampleA
                                                         , dist.fun = dist.fun)
  adam.sample.collectionB[[cc]] <- adam_function.original(sampleB
                                                          , dist.fun = dist.fun)
  TNext <- balancerFun(adam.sample.collectionA[[cc]]$TNext
                       , adam.sample.collectionB[[cc]]$TNext)
  adam.sample.collectionA[[cc]]$TNext <- TNext
  adam.sample.collectionB[[cc]]$TNext <- TNext
  while( ii + adam.sample.collectionA[[cc]]$TNext - 1 < length(metricA) ){
    prev.ii <- ii
    cc <- cc + 1
    ii <- ii + adam.sample.collectionA[[cc-1]]$TNext
    sampleA <- metricA.lst[[ii]]
    sampleB <- metricB.lst[[ii]]
    adam.objA <- adam.sample.collectionA[[cc-1]]
    adam.objB <- adam.sample.collectionB[[cc-1]]
    adam.sample.collectionA[[cc]] <- adam_function.generic(sampleA
                                                          , adam.objA
                                                          , dist.fun = dist.fun)
    adam.sample.collectionB[[cc]] <- adam_function.generic(sampleB
                                                           , adam.objB
                                                           , dist.fun = dist.fun)
    TNext <- balancerFun(adam.sample.collectionA[[cc]]$TNext
                         , adam.sample.collectionB[[cc]]$TNext)
    adam.sample.collectionA[[cc]]$TNext <- TNext
    adam.sample.collectionB[[cc]]$TNext <- TNext
  }
  # RETURN
  do.call(rbind.data.frame, adam.sample.collectionA)
}

#
#bandit balancing
#
AdvancedBandit1 <- function (n) {
  if (n < 0) stop("narms must be > 0")
  value <- list(narms=n, cnts = matrix(0, ncol = n), vals = matrix(0, ncol = n), alpha = 0.4)
  
  attr(value, "class") <- "AdvancedBandit1"
  
  return(value)
}

selectArm <- function(bandit) {UseMethod("selectArm")}
updateBandit <- function(bandit, chosen_arm, reward) {UseMethod("updateBandit")}

selectArm.AdvancedBandit1 <- function(bandit) {
#  for(arm in seq_along(bandit$cnts)) 
#    if (bandit$cnts[arm] == 0)
#      return(arm)
  
  ucb_values <- matrix(0, ncol = bandit$narms)
  total_cnts = sum(bandit$cnts)
  
  for(arm in seq_along(bandit$cnts)){
    bonus = sqrt((2 * log(total_cnts)) / bandit$cnts[arm])
    if (bonus == Inf)
      bonus = 0
    ucb_values[arm] = bandit$vals[arm] + bonus
  }
  
  print(bandit$vals)
  print(ucb_values)
  print(which.max(ucb_values))
  print(" ")
  return (which.max(ucb_values))
}

updateBandit.AdvancedBandit1 <- function(bandit, chosen_arm, reward) {
  bandit$cnts[chosen_arm] <- bandit$cnts[chosen_arm] + 1
  n <- bandit$cnts[chosen_arm]
  
  old_v <- bandit$vals[chosen_arm]
  #new_v <- ((n - 1) / n) * old_v + (1 / n) * reward
  #EWMA
  new_v <- (1 - bandit$alpha) * old_v + bandit$alpha* reward
  bandit$vals[chosen_arm] <- new_v
  
  return(bandit)
}

AdvancedBandit2 <- function (n, tau = 0.1) {
  if (n < 0 || tau < 0.0) stop("narms must be > 0 and tau must be in range [0,1]")
  value <- list(narms=n, tau=tau, cnts = matrix(0, ncol = n), vals = matrix(0, ncol = n), alpha = 0.3)
  
  attr(value, "class") <- "AdvancedBandit2"
  
  return(value)
}

selectArm.AdvancedBandit2 <- function(bandit) {
  t = sum(bandit$cnts) + 1
  bandit$tau = 1 / log(t + 0.0000001)
  
  z <- sum(exp(bandit$vals/bandit$tau))
  probs <- exp(bandit$vals/bandit$tau)/z
  
  arm <- categorical_draw(probs)
  
  print(bandit$vals)
  print(probs)
  print(arm)
  print(" ")
  
  return (arm)
}

updateBandit.AdvancedBandit2 <- function(bandit, chosen_arm, reward) {
  bandit$cnts[chosen_arm] <- bandit$cnts[chosen_arm] + 1
  n <- bandit$cnts[chosen_arm]
  
  old_v <- bandit$vals[chosen_arm]
  new_v <- (1 - bandit$alpha) * old_v + bandit$alpha* reward
  
  bandit$vals[chosen_arm] <- new_v
  
  return(bandit)
}

categorical_draw <- function(probs) {
  z <- runif(1)
  cum_prob <- 0.0
  for(i in seq_along(probs)) {
    cum_prob <- cum_prob + probs[i]
    if (cum_prob > z)
      return(i)
  }
  
  return (length(probs))
  
}

dynSampling2 <- function(metricA, metricB, balancerFun=min) {
  dist.fun <- abs.diff
  t.seq <- seq(0, length(metricA)-1)
  metricA.lst <- split(data.frame(t=t.seq, v=metricA), t.seq)
  metricB.lst <- split(data.frame(t=t.seq, v=metricB), t.seq)
  cc <- 1
  ii <- 1
  adam.sample.collectionA <- list()
  adam.sample.collectionB <- list()
  sampleA <- metricA.lst[[1]]
  sampleB <- metricB.lst[[1]]
  adam.sample.collectionA[[cc]] <- adam_function.original(sampleA
                                                          , dist.fun = dist.fun)
  adam.sample.collectionB[[cc]] <- adam_function.original(sampleB
                                                          , dist.fun = dist.fun)
  
  TNext <- balancerFun(adam.sample.collectionA[[cc]]$TNext
                       , adam.sample.collectionB[[cc]]$TNext)
  adam.sample.collectionA[[cc]]$TNext <- TNext
  adam.sample.collectionB[[cc]]$TNext <- TNext
  
  training <- 10
  change_thres <- 0.5
  tt <- 1
  
  while( ii + adam.sample.collectionA[[cc]]$TNext - 1 < length(metricA) ){
    prev.ii <- ii
    cc <- cc + 1
    ii <- ii + adam.sample.collectionA[[cc-1]]$TNext
    sampleA <- metricA.lst[[ii]]
    sampleB <- metricB.lst[[ii]]
    adam.objA <- adam.sample.collectionA[[cc-1]]
    adam.objB <- adam.sample.collectionB[[cc-1]]
    adam.sample.collectionA[[cc]] <- adam_function.generic(sampleA
                                                           , adam.objA
                                                           , dist.fun = dist.fun)
    adam.sample.collectionB[[cc]] <- adam_function.generic(sampleB
                                                           , adam.objB
                                                           , dist.fun = dist.fun)
    
    #suggested period based on temporal evolution
    mu_period <- adam.sample.collectionA[[cc]]$TNext
    mu_conf <- adam.sample.collectionA[[cc]]$conf
    
    tau_evolution <- adam.sample.collectionB[[cc]]$mean
    tau_period <- adam.sample.collectionB[[cc]]$TNext
    tau_conf <- adam.sample.collectionB[[cc]]$conf
    
    #TODO instead of change_tres use ADMin online change detection
    #compare ADMin to proportional approach with round(mu_period * (1 - tau_evolution))
    #ADMin will have larger compression rate but slightly increased error opposed to proportional approach
    if (tt < training || tau_evolution > change_thres)
      TNext <- 1
    else
      TNext <- round(mu_period * (1 - tau_evolution))
      
    
    tt <- tt + 1
    
    #print(paste("t: ", ii, " T: ", TNext,  " mu_period: ", mu_period, " mu_conf: ", mu_conf, " tau_evolution: ", tau_evolution, " tau_conf: ", tau_conf))
    
   # TNext <- balancerFun(mu_period, tau_period)
    
    adam.sample.collectionA[[cc]]$TNext <- TNext
    adam.sample.collectionB[[cc]]$TNext <- TNext
  }
  
  # RETURN
  do.call(rbind.data.frame, adam.sample.collectionA)
}