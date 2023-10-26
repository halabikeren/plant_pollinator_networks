#!/usr/bin/env Rscript
require(tidyverse)
require(bipartite)
require(glue)
source("utils.R")

args = commandArgs(trailingOnly=TRUE)

if (length(args)<2) {
  stop("Two arguments (input path and output path) must be supplied", call.=FALSE)
}

input_path = args[1]
null_dir = args[2]
null_method = args[3]
output_path = args[4]
only_input = args[5]
print(glue("input_path=", input_path))
print(glue("output_path=", output_path))

start_time = Sys.time()
if (only_input) {
  network = process_network(input_path)
  is_weighted = any(network > 1)
  network_features = networklevel(web=network, weighted=is_weighted)
  network_features[["modularity"]] = get_modularity(network)
} else {
  null_str = glue("_across_null_networks_", null_method, ".csv")
  null_sim_features_path = str_replace(output_path, ".csv", null_str)
  features = get_network_features(input_path, null_sim_features_path, null_dir)
}
features["network"] = basename(input_path)
network_features = t(data.frame(Reduce(rbind, features)))
colnames(network_features) = names(features)
write.csv(network_features, output_path, row.names = TRUE)
end_time = Sys.time()
duration = end_time - start_time
print(str_glue("duration = {duration}"))