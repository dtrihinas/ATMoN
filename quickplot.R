library("zoo") # na.locf()

#quick plot original datasets and metrics

#plot options
par(mfrow=c(2,1))

#original stream of dataset
qdataset <- "reality.mit"
qdataset.pretty <- "mit"
qmetric <- "dissimilarity"
qmetric.pretty <- "graph topology structure volatility"

a_original <- gstream.collection[[qdataset]]$metrics[[qmetric]]
a_original[is.na(a_original)] <- 0

#for easy visualization
timespan = list(tstart = 1, tend = length(a_original))
#timespan = list(tstart = 1, tend = 50000)


plot(a_original[timespan$tstart:timespan$tend], xlab="time intervals", ylab=qmetric, type="l",
     main=paste("dataset: ", qdataset.pretty, "\tmetric: ", qmetric.pretty) 
)

a_mstream <- adam_sampling.standalone(a_original[timespan$tstart:timespan$tend])
x_mseq <- cumsum(a_mstream$TNext)

mstream_optratio <- get_comp_ratio(a_mstream$TNext, a_original[timespan$tstart:timespan$tend])
mstream_err <- get_mape_error(a_mstream$TNext, a_original[timespan$tstart:timespan$tend])

plot(x_mseq, a_original[x_mseq], xlab="time intervals", ylab=qmetric, type="l",
     main=paste("dataset: ", qdataset.pretty, "\tmetric: ", qmetric.pretty) 
)


a_topo <- gstream.collection[[qdataset]]$metrics[["dissimilarity"]]

plot(a_topo[timespan$tstart:timespan$tend], xlab="time intervals", ylab="structure volatility", type="l",
     main=paste("dataset: ", qdataset.pretty, "\tmetric: ", "Graph Structure Volatility (D-measure)") 
)
#dev.copy(png, paste(qdataset.pretty,"-", qmetric.pretty, ".png"))
#dev.off()