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
null_sim_features_path = str_replace(output_path, ".csv", "_across_null_networks.csv")
species_features = get_species_features(input_path, null_sim_features_path, level="lower", nsim=1000)
species_features["network"] = basename(input_path)
write.csv(species_features, output_path, row.names = TRUE)
end_time <- Sys.time()
duration = end_time - start_time
print(str_glue("duration = {duration}"))
