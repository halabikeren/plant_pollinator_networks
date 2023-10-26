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

start_time = Sys.time()

network = process_network(input_path)
network[is.na(network)] = 0
row_names = rownames(network)
species_features = data.frame(specieslevel(web=network, level="lower"))
row.names(species_features) = row_names
species_features["network"] = basename(input_path)
write.csv(species_features, output_path, row.names = TRUE)
end_time = Sys.time()
duration = end_time - start_time
print(str_glue("duration = {duration}"))
