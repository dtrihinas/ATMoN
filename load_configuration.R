
#
# Author: Luis F. Chiroque
# e-mail: lf.chiroque@imdea.org
# IMDEA Networks Institute
#


# Global variables
## whether the metrics should be recomputed on graphs list
redo <- FALSE # not implemented yet!
## number of cores to be used for computing metrics
options("mc.cores"=10)
## R-object file name to store graphs list and metrics
Robj.fileName <- "gstream.collection-lst.Robj"
## relative path to datasets
datasets.folderName <- "datasets/"
## R-object name for graphs list collection
gstreamObjectName <- "gstream.collection"
## T_min and T_max default
time.thrs <- list(TMin=1, TMax=15)
