#!/usr/bin/env Rscript
library(tidyverse)
library(bipartite)
source("utils.R")

args = commandArgs(trailingOnly=TRUE)

if (length(args)<2) {
  stop("Two arguments (input path and output path) must be supplied", call.=FALSE)
}

input_path = args[1]
output_path = args[2]

start_time <- Sys.time()
interaction_features = data.frame(matrix(ncol = 7, nrow = 0))
colnames(interaction_features) = c("network", "plant", "pollinator", "plant_dependence_on_pollinator", "plant_interactions_num", "pollinator_dependence_on_plant", "pollinator_interactions_num")
features = get_interaction_features(input_path)
interaction_features[nrow(interaction_features) + 1,] = c(c(basename(input_path)), features)
write.csv(interaction_features, output_path, row.names = TRUE)
end_time <- Sys.time()
duration = end_time - start_time
print(str_glue("duration = {duration}"))