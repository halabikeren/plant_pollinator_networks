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
community_features = data.frame(matrix(ncol = 5, nrow = 0))
colnames(community_features) = c("network", "plant_richness", "pollinator_richness", "total_species_richness", "richness_ratio")
features = get_community_features(input_path)
community_features[nrow(community_features) + 1,] = c(c(basename(input_path)), features)
write.csv(community_features, output_path, row.names = TRUE)
end_time <- Sys.time()
duration = end_time - start_time
print(str_glue("duration = {duration}"))