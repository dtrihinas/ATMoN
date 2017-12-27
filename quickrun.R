#new script execution

library("zoo") # na.locf()


#a_data <- list(dataset = "sg.09", metric = c("diameter", "nofComponents", "gcRatio", "deg.dist.pow"), topo="dissimilarity")
#a_data <- list(dataset = "reality.mit", metric = c("deg.dist.pow", "esize", "avgPageRank", "max.clique", "gcRatio"), topo="dissimilarity")
a_data <- list(dataset = "vehicular.net", metric = c("deg.dist.pow", "diameter", "avgPageRank", "max.clique", "gcRatio"), topo="dissimilarity")


#plot options
par(mfrow=c(5,1))

get_comp_ratio <- function (metricFrameSample, metricFrameBase) {
  (1 - length(metricFrameSample) / length(metricFrameBase))*100
}

get_mape_error <- function(metricFrameSample, metricFrameBase) {
  f <- rep(NA, length(metricFrameBase))
  s <- cumsum(metricFrameSample)
  s <- s[which(s <= length(metricFrameBase))]
  f[s] <- metricFrameBase[s]
  a <- na.locf(f)
  
  sum(abs(metricFrameBase - a)) / sum(metricFrameBase + a) * 100
}

for (metric in a_data$metric) {

  #original stream of dataset
  a_original <- gstream.collection[[a_data$dataset]]$metrics[[metric]]
  
  #for easy visualization
  timespan = list(tstart = 1, tend = length(a_original))
  
  plot(a_original[timespan$tstart:timespan$tend], xlab="time intervals", ylab=metric, type="l",
       main=paste("dataset: ", a_data$dataset, "\tmetric: ", metric)
       )
  
  #apply AdaM to get metric stream
  a_mstream <- adam_sampling.standalone(a_original[timespan$tstart:timespan$tend])
  x_mseq <- cumsum(a_mstream$TNext)
  
  mstream_optratio <- get_comp_ratio(a_mstream$TNext, a_original[timespan$tstart:timespan$tend])
  mstream_err <- get_mape_error(a_mstream$TNext, a_original[timespan$tstart:timespan$tend])
  
  plot(x_mseq, a_original[x_mseq], xlab="time intervals", ylab=metric, type="l",
       main=paste("AdaM applied\n", 
                  "compression ratio: ", round(mstream_optratio), "%\t",
                  "error: ", round(mstream_err), "%"
                  )
        )
  
  #original topology stream of dataset and remove NaN
  a_original_topo <- gstream.collection[[a_data$dataset]]$metrics[[a_data$topo]]
  a_original_topo[is.na(a_original_topo)] <- 0
  
  #apply AdaM to get topology stream
  a_tstream <- adam_sampling.standalone(a_original_topo[timespan$tstart:timespan$tend])
  x_tseq <- cumsum(a_tstream$TNext)
  
  tstream_optratio <- get_comp_ratio(a_tstream$TNext, a_original[timespan$tstart:timespan$tend])
  tstream_err <- get_mape_error(a_tstream$TNext, a_original[timespan$tstart:timespan$tend])
  
  plot(x_tseq, a_original[x_tseq], xlab="time intervals", ylab=metric, type="l",
       main=paste("AdaM applied on topology\n",
                  "compression ratio: ", round(tstream_optratio), "%\t",
                  "error: ", round(tstream_err), "%"
                  )
       )
  
  #apply AdaM with balancing function on metric and topology stream
  a_mtfinal <-dynSampling(a_original[timespan$tstart:timespan$tend], a_original_topo[timespan$tstart:timespan$tend], balancerFun = min)
  x_mtfinal_seq <- cumsum(a_mtfinal$TNext)
  
  mtstream_optratio <- get_comp_ratio(a_mtfinal$TNext, a_original[timespan$tstart:timespan$tend])
  mtstream_err <- get_mape_error(a_mtfinal$TNext, a_original[timespan$tstart:timespan$tend])
  
  plot(x_mtfinal_seq, a_original[x_mtfinal_seq], xlab="time intervals", ylab=metric, type="l",
       main=paste("AdaM applied on both metric and topology stream with balancing (min)\n", 
                  "compression ratio: ", round(mtstream_optratio), "%\t",
                  "error: ", round(mtstream_err), "%"
                  )
       )
  
  
  #apply AdaM with advanced balancing function on metric and topology stream
  a_mtfinal2 <-dynSampling2(a_original[timespan$tstart:timespan$tend], a_original_topo[timespan$tstart:timespan$tend], balancerFun = min)
  x_mtfinal_seq2 <- cumsum(a_mtfinal2$TNext)
  
  mtstream_optratio2 <- get_comp_ratio(a_mtfinal2$TNext, a_original[timespan$tstart:timespan$tend])
  mtstream_err2 <- get_mape_error(a_mtfinal2$TNext, a_original[timespan$tstart:timespan$tend])
  
  plot(x_mtfinal_seq2, a_original[x_mtfinal_seq2], xlab="time intervals", ylab=metric, type="l",
       main=paste("Advanced balancing with AdaM incorportating topology structure change predictor (D-Measure)\n", 
                  "compression ratio: ", round(mtstream_optratio2), "%\t",
                  "error: ", round(mtstream_err2), "%"
       )
  )
  
  dev.copy(png, paste(a_data$dataset,"-", metric, ".png"))
  dev.off()
  
  
  a_method <- c("AdaM_orignal", "AdaM_only_topo", "AdaM_balancing_min", "AdaM_topo_integrated")
  a_comp <- round(c(mstream_optratio, tstream_optratio, mtstream_optratio, mtstream_optratio2), digits=2)
  a_error <- round(c(mstream_err, tstream_err, mtstream_err, mtstream_err2), digits=2)
  a_title <-paste("dataset: ", a_data$dataset, " metric: ", metric, " Tmin = 1, Tmax = 15")
  
  df <- data.frame(a_method, a_comp, a_error)
  colnames(df) <- c("method", "compression (%)", "error (%)")
#  print(a_title)
#  print(df)
  
  sink(paste(a_data$dataset,".txt"), append = TRUE)
  cat("\n\n")
  cat(a_title)
  cat("\n\n")
  print(df)
  sink()
  
}
