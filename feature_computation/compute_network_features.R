#!/usr/bin/env Rscript
library(tidyverse)
library(bipartite)
source("/groups/itay_mayrose/halabikeren/tmp/plant_pollinator_inter/feature_computation/utils.R")

args = commandArgs(trailingOnly=TRUE)

if (length(args)<2) {
  stop("Two arguments (input path and output path) must be supplied", call.=FALSE)
}

input_path = args[1]
output_path = args[2]

start_time <- Sys.time()
network_features = data.frame(matrix(ncol = 6, nrow = 0))
colnames(network_features) = c("network", "connectance", "avg_plant_inter", "avg_pollinator_inter", "NODF", "relative_NODF")
features = get_network_features(input_path)
network_features[nrow(network_features) + 1,] = c(c(basename(input_path)), features)
write.csv(network_features, output_path, row.names = TRUE)
end_time <- Sys.time()
duration = end_time - start_time
print(str_glue("duration = {duration}"))