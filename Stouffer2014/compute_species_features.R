#!/usr/bin/env Rscript
library(tidyverse)
library(bipartite)
source("utils.R")

args = commandArgs(trailingOnly=TRUE)

if (length(args)<2) {
  stop("Two arguments (input path and output path) must be supplied", call.=FALSE)

input_path = args[1]
output_path = args[2]

species_features = data.frame(matrix(ncol = 4, nrow = 0))
colnames(species_features) = c("network", "species", "rank", "nestedness_contribution")
network_species_features = get_species_features(input_path)
species_features = rbind(species_features, network_species_features)
write.csv(species_features, output_path, row.names = TRUE)