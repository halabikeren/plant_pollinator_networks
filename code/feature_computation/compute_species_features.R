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

start_time <- Sys.time()
null_sim_features_path = str_replace(output_path, ".csv", "_across_null_networks.csv")
species_features = get_species_features(input_path, null_sim_features_path, null_dir, level="lower")
species_features["network"] = basename(input_path)
write.csv(species_features, output_path, row.names = TRUE)
end_time <- Sys.time()
duration = end_time - start_time
print(str_glue("duration = {duration}"))
