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
  network[is.na(network)] = 0
  is_weighted = any(network > 1)

  if (level == "lower")
  {
    row_names = rownames(network)
  } else {
    row_names = colnames(network)
  }
  species_features = data.frame(specieslevel(web=network, level=level))
  row.names(species_features) = row_names
  if (level == "lower")
  {
    species_features["nestedness_contribution"] = get_plant_nestedness_contribution(web=network, nsimul=nsim)
  } else {
    species_features["nestedness_contribution"] = get_pollinator_nestedness_contribution(web=network, nsimul=nsim)
  }
} else {
  null_sim_features_path = str_replace(output_path, ".csv", "_across_null_networks.csv")
  species_features = get_species_features(input_path, null_sim_features_path, null_dir, level="lower")
}
species_features["network"] = basename(input_path)
write.csv(species_features, output_path, row.names = TRUE)
end_time = Sys.time()
duration = end_time - start_time
print(str_glue("duration = {duration}"))
