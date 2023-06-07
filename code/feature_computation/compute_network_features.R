#!/usr/bin/env Rscript
library(tidyverse)
library(bipartite)
source("utils.R")

args = commandArgs(trailingOnly=TRUE)

if (length(args)<2) {
  stop("Two arguments (input path and output path) must be supplied", call.=FALSE)
}

input_path = args[1]
null_dir = args[2]
output_path = args[3]
only_input = args[4]

start_time = Sys.time()
if (only_input) {
  network = process_network(input_path)
  is_weighted = any(network > 1)
  network_features = networklevel(web=network, weighted=is_weighted)
  network_features[["modularity"]] = get_modularity(network)
  extinction_features = simulate_network_extinction(network, nsim=nsim)
  features = c(network_features, extinction_features)
} else {
  null_sim_features_path = str_replace(output_path, ".csv", "_across_null_networks.csv")
  features = get_network_features(input_path, null_sim_features_path, null_dir)
}
features["network"] = basename(input_path)
network_features = t(data.frame(Reduce(rbind, features)))
colnames(network_features) = names(features)
write.csv(network_features, output_path, row.names = TRUE)
end_time = Sys.time()
duration = end_time - start_time
print(str_glue("duration = {duration}"))