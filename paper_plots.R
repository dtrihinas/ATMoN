
#
# Author: Luis F. Chiroque
# e-mail: lf.chiroque@imdea.org
# IMDEA Networks Institute
#

# TODO: code not working 100%

library(zoo) # na.fill()
# library("tikzDevice") # tikz()
# options(tikzMetricPackages = c("\\usepackage[utf8]{inputenc}",
#                                "\\usepackage[T1]{fontenc}", "\\usetikzlibrary{calc}",
#                                "\\usepackage{amsmath}", "\\usepackage{amsfonts}", "\\usepackage{amssymb}"))

outputDir <- "figures/"
Robj.name <- "gstream.collection"
Robj.fileName <- "gstream.collection-lst.Robj"

if (!file.exists(outputDir)) {
  dir.create(outputDir)
}

if (!exists(Robj.name)) {
  load(Robj.fileName)
}


# Main figures

## Topology correlations

get_correlations <- function(topologyName, ds.obj) {
  tplgy <- ds.obj$metrics[[topologyName]]
  tplgy.mean <- sapply(1:length(tplgy), function(k) mean(tplgy[1:k], na.rm=T))
  tplgy.n <- ds.obj$metrics$vsize
  #kk.sd <- sapply(1:length(tplgy), function(k) sd(tplgy[1:k], na.rm=T))
  tplgy.absmeandiff <- abs(tplgy - tplgy.mean)
  tplgy.prop <- tplgy.absmeandiff * tplgy.n
  tplgy.norm <- tplgy.absmeandiff / tplgy.n
  tplgy.df <- data.frame(tplgy, tplgy.mean, tplgy.absmeandiff, tplgy.prop, tplgy.norm, tplgy.n)
  tplgy.df <- data.frame(na.fill(tplgy.df, 0))
  sapply(tplgy.df, function(tplgy)
    sapply(ds.obj$metrics[which(!sapply(ds.obj$metrics, is.list))], function(mm)
      cor(tplgy, mm, use="complete.obs")
    )
  )
}

get_topology_correlations <- function(ds.obj, target.topologies = c("assortativity", "transitivity", "dissimilarity")) {
  lapply(lapply(target.topologies, get_correlations, ds.obj)
         , function(corr.tbl) {
           metrics.filter <- which(!rownames(corr.tbl) %in% c(target.topologies, "vsize", "pageRank", "outDegree", "cluster.mem"))
           corr.tbl[metrics.filter,]
         })
}

target.topologies <- c("assortativity", "transitivity", "dissimilarity")
metric.corrs <- lapply(gstream.collection, get_topology_correlations
                       , target.topologies)
metric.corrs.aggreg <- lapply(seq_along(target.topologies), function(ii) {
  do.call(rbind, lapply(metric.corrs, function(mcorr) abs(mcorr[[ii]])))
})
names(metric.corrs.aggreg) <- target.topologies

par.bak <- par()
# tikz(paste0(outputDir, "corr-tplgy_boxplot.tex")
#      , standAlone = TRUE, width=50, height=25
#      , packages = c("\\usepackage{tikz}",
#                     "\\usepackage[active,tightpage,psfixbb]{preview}",
#                     "\\PreviewEnvironment{pgfpicture}",
#                     "\\setlength\\PreviewBorder{0pt}",
#                     "\\usepackage{amssymb}"))
pdf(paste0(outputDir, "corr-tplgy_boxplot.pdf"), width = 12, height = 6)
par(las=1, mfrow=c(1, length(target.topologies)), mar=c(7,4,4,2)+.1)
sapply(target.topologies, function(topologyName) {
  m <- ncol(metric.corrs.aggreg[[topologyName]])
  boxplot(metric.corrs.aggreg[[topologyName]], ylim=c(0,1)
          , main=topologyName, ylab="correlation"
          , xaxt='n'
          #, names=rep("$\\mathcal{T}$", 5))#, "$\\overline{\\mathcal{T}}$"
          #, "$|\\mathcal{T} - \\overline{\\mathcal{T}}|$"
          #, "$|\\mathcal{T} - \\overline{\\mathcal{T}}| / |V|$"
          #, "$|\\mathcal{T} - \\overline{\\mathcal{T}}| · |V|$"))
  )
  abline(h=.5, lty=2)
  axis(1,1:m, labels = NA)
  text(1:m, par("usr")[3] - 0.1
       , colnames(metric.corrs.aggreg[[1]])
       , srt=45, xpd=T, adj=c(1,0))
  invisible()
})
dev.off()
# tools::texi2dvi(paste0(outputDir, "corr-tplgy_boxplot.tex"), pdf=T)
par(par.bak)



## compression; data volume reduction
compression.lst <- do.call(c, sapply(perf.top
                                     , sapply
                                     , function(tplgy) tplgy[,"compression.topology"])
                           )
hist(compression.lst)
summary(compression.lst)
## accuracy; (1 - error)
error.lst <- do.call(c, sapply(perf.top
                               , sapply
                               , function(tplgy) 1 - tplgy[,"error.topology"])
                     )
hist(error.lst)
summary(error.lst)
##
eff.lst <- do.call(c, sapply(perf.top
                             , sapply
                             , function(tplgy) tplgy[,"efficiency.topology"])
                   )


# Illustrative adaptive sampling
set.seed(5)
N <- 50

metric <- rnorm(N,2) + 5
#metric <- sin(seq(0,3*pi, length=N) + .2) + rnorm(N,0,.2) + 5
metric.adam.df <- adam_sampling.standalone(metric)

metric.adam.xseq <- cumsum(metric.adam.df$TNext)
metric.adam.xseq <- metric.adam.xseq[which(metric.adam.xseq <= length(metric))]
avgT <- mean(metric.adam.df$TNext)
equid.xseq <- seq(1,100, by=round(avgT))
              # cumsum(sample(floor(avgT):ceiling(avgT)
              #               , nrow(metric.adam.df)
              #               , replace=T, prob=c(ceiling(avgT) - avgT, avgT - floor(avgT))))
equid.xseq <- equid.xseq[which(equid.xseq <= length(metric))]
## plot comparing the sampling points
plot(metric, type="l", ylab="", xaxt="n", yaxt="n", xlab="sample metric")
abline(h=mean(metric), lty=2)
lines(metric.adam.xseq, metric[metric.adam.xseq], lty=3, col="blue")
lines(equid.xseq, metric[equid.xseq], lty=3, col="red")

adam.metric <- rep(NA, length(metric))
adam.metric[metric.adam.xseq] <- metric[metric.adam.xseq]
adam.metric <- na.locf(adam.metric)
equid.metric <- rep(NA, length(metric))
equid.metric[c(1,equid.xseq)] <- metric[c(1,equid.xseq)]
equid.metric <- na.locf(equid.metric)
## plot actual comparison
pdf(paste0(outputDir, "benchmark_comparison.pdf"), width=16, height=9)
par.bak <- par()
par(mar=c(5,4,4,15) + .1)
plot(metric, type="o", lwd=2, cex=.75
     , ylab="", xaxt="n", yaxt="n"
     , xlab="time")
abline(h=mean(metric), lty=6, lwd=1.5)
lines(adam.metric, lty=2, col="blue", lwd=1.5, type="s")
lines(equid.metric, lty=2, col="red", lwd=1.5, type="s")
par(xpd=TRUE)
legend("right", lty=c(1,6,2,2), inset=c(-.2, 0), y.intersp=1
       , legend=c("stream", "average", "AdaM", "equidist")
       , col=c(1,1,"blue", "red"), lwd=3, cex=1.5, text.width = 5)
dev.off()
par(par.bak)

sum(abs(metric - adam.metric)) / sum(metric + adam.metric)
sum(abs(metric - equid.metric)) / sum(metric + equid.metric)
sum(abs(metric - mean(metric))) / sum(metric + mean(metric))




# Error barplots
error_barplot <- function(ds.name, perf.top, legend=F, args.leg=list(cex=.5, x="top", ncol=2)) {
  err <- lapply(names(perf.top[[ds.name]]), function(tplgy.name) {
    sapply(rownames(perf.top[[ds.name]][[tplgy.name]]), function(metric.name) {
      perf.top[[ds.name]][[tplgy.name]][metric.name, c("error", "error.topology"
                                                       , "avg.error", "equid.error", "rnd.error", "rnd.sd")]
    })
  })
  labels <- colnames(err[[1]])
  err.se <- err[[1]][6,]
  if (length(err) == 1 || ncol(err[[1]]) == 1) {
    #TODO: special case length(err) == 1
    err <- t(t(c(as.matrix(err[[1]][c(3,5,4),])
               , as.matrix(err[[1]][1,])
               , sapply(err, function(errMtx) {errMtx[2,]}))))
  } else {
    err <- rbind(as.matrix(err[[1]][c(3,5,4),])
                 , t(as.matrix(err[[1]][1,]))
                 , t(sapply(err, function(errMtx) {errMtx[2,]})))
  }
  rownames(err) <- c("Avg. heuristic", "random", "equidistant"
                      , paste0("temporal", c("", paste0("+", names(perf.top[[ds.name]])))))
  M <- nrow(err)
  barplot(err*100, beside=T, ylim=c(0,55), ylab="error (%)"
          , legend.text=legend, args.legend=args.leg
          , main=ds.name, xaxt="n")
  text(1:ncol(err)*M-.5, par("usr")[3] - 3
       , labels, cex=.75
       , srt=45, xpd=T, adj=c(1,0))
  centers <- (1:ncol(err) - 1) * (M + 1) + 2.5
  segments(centers, err["random",]*100 + err.se*200
           , centers, err["random",]*100 - err.se*200
           , lwd=1.5)
  arrows(centers, err["random",]*100 + err.se*200
         , centers, err["random",]*100 - err.se*200
         , angle=90, lwd=1.5, code = 3, length=.09 * .9 / ncol(err))
  invisible()
}

#sapply(performance, function(perf.top) {
  perf.top <- performance$top5
  sapply(names(perf.top), error_barplot, perf.top
         , legend=T, args.leg=list(cex=.5, x="topleft", ncol=2))
#})

par.bak <- par()
pdf(paste0(outputDir, "error_barplots.pdf"))
par(mfrow=c(3,3), las=2)
sapply(names(perf.top)[1:4], error_barplot, perf.top)
par(xpd=T)
frame()
lgd.text <- c("Avg. heuristic", "equidistant", "random"
              , paste0("temporal", c("", paste0("+", names(perf.top[[ds.name]])))))
legend("center", lgd.text, ncol = 2, cex=.75
       , fill=grey.colors(8))
par(xpd=par.bak$xpd)
sapply(names(perf.top)[5:8], error_barplot, perf.top)
dev.off()
par(par.bak)


adam.topology.result.df <- lapply(names(perf.top[[1]]), function(topology.name) {
  do.call(rbind
         , sapply(perf.top, function(ds.obj) {
    as.matrix(as.data.frame(ds.obj[[topology.name]])[,c("error.topology", "compression.topology")])
  }))
})

par(mfrow=c(1,1))
plot(NULL, type='n', xlim=c(0,1), ylim=c(0,.22)
     , ylab="error", xlab="redundancy")
sapply(1:length(adam.topology.result.df), function(i) {
  xy <- adam.topology.result.df[[i]]
  print(names(xy))
  points(1-xy[,2], xy[,1], pch=i)
})
