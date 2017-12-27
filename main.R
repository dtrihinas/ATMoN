
source("load_configuration.R")

# Load datasets/graphs list
source("load_graphStreams.R")
# Compute metrics on graphs list
source("compute_gstream_metrics.R")
# temporal graphs (dataset) are loaded as a list of graphs
# and 'gstream.collection' has a list of temporal graphs
# plus a list of computed metrics for each dataset

# Dyn monitoring
source("adam_sampling.R") # it yields: 'results', 'res.optim', 'res.rand' and 'res.equid'

