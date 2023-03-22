#!/usr/bin/env Rscript
library(tidyverse)
library(bipartite)
source("utils.R")

args = commandArgs(trailingOnly=TRUE)

if (length(args)<1) {
  stop("One argument must be supplied", call.=FALSE)
}



get_module_sharing <-function(member1, member2)
{
  m1_modules = modules_matrix[, which(member1 == rownames(network))+2] > 0
  m2_modules = modules_matrix[, which(member2 == rownames(network))+2] > 0
  frac_shared_modules = sum(m1_modules&m2_modules) / ncol(modules_matrix)
  return(frac_shared_modules)
}

pairs_path = args[1]

start_time <- Sys.time()

pairs_data = read.csv(pairs_path)
network_path = str_replace(pairs_data$network_path[1], "../../../", "../../")
network = process_network(network_path)

modules = computeModules(network)
modules_matrix = as.matrix(slot(modules, "modules"))
pairs_data["n_modules"] = nrow(modules_matrix)
pairs_data["frac_shared_modules"] = mapply(get_module_sharing,
                                           pairs_data$member1,
                                           pairs_data$member2)
write.csv(pairs_data, pairs_path, row.names = TRUE)

end_time <- Sys.time()
duration = end_time - start_time
print(str_glue("duration = {duration}"))
