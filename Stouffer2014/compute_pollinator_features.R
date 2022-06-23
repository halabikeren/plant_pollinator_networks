#!/usr/bin/env Rscript
library(tidyverse)
library(bipartite)
source("/groups/itay_mayrose/halabikeren/tmp/plant_pollinator_inter/Stouffer2014/utils.R")

args = commandArgs(trailingOnly=TRUE)

if (length(args)<2) {
  stop("Two arguments (input path and output path) must be supplied", call.=FALSE)
}

input_path = args[1]
output_path = args[2]
species_classification_path = args[3]

start_time <- Sys.time()
species_features = data.frame(matrix(ncol = 5, nrow = 0))
colnames(species_features) = c("network","pollinator", "ranked_degree", "ranked_nestedness_contribution", "exotic_tendency")
network_species_features = get_pollinator_features(input_path, species_classification_path)
species_features = rbind(species_features, network_species_features)
write.csv(species_features, output_path, row.names = TRUE)
end_time <- Sys.time()
duration = end_time - start_time
print(str_glue("duration = {duration}"))