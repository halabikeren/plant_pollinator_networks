library(tidyverse)
library(bipartite)

process_network <- function(network_path)
{
  network <- read_csv(file = network_path, show_col_types = FALSE)
  network$is_exotic = NULL
  network = network %>% select(., -"...1") %>% column_to_rownames(., var = "plant_name")
  network = 1*(network >= 1)
  return(network)
}

get_community_features <- function(network_path)
{
  network <- process_network(network_path)
  plant_richness = nrow(network)
  pollinator_richness = ncol(network)-1 # last column is exotic classification
  total_species_richness = plant_richness + pollinator_richness
  richness_ratio = pollinator_richness/plant_richness
  result <- list("plant_richness" = plant_richness,
                 "pollinator_richness" = pollinator_richness,
                 "total_species_richness" = total_species_richness,
                 "richness_ratio" = richness_ratio)
  return(result)
}

get_interaction_based_features <- function(network)
{
  L = sum(network > 0)
  P = nrow(network)
  A = ncol(network)
  res = list("connectance" = L/(P*A), "avg_plant_inter" = L/P, "avg_pollinator_inter" = L/A)
  return(res)
}

sample_networks <- function(network, nsim)
{
  row_probs = matrix(rowSums(network)/ncol(network), nrow=nrow(network), ncol=ncol(network))
  col_probs = matrix(colSums(network)/nrow(network), nrow=nrow(network), ncol=ncol(network), byrow=TRUE)
  cell_probs = (row_probs + col_probs)/2

  num_ones = matrix(0, nrow=nrow(network), ncol=ncol(network))

  for (i in 1:nsim)
  {
    null_network = matrix(NA, ncol=ncol(network), nrow=nrow(network))
    null_network[] <- rbinom(n=ncol(network)*nrow(network), size=1, prob = cell_probs)
    num_ones = num_ones + null_network

  }
  sampled_probs = num_ones/nsim
  return(list(cell_probs, sampled_probs))
}

get_nestedness <- function (network, nsim=1000)
{
  row_probs = matrix(rowSums(network)/ncol(network), nrow=nrow(network), ncol=ncol(network))
  col_probs = matrix(colSums(network)/nrow(network), nrow=nrow(network), ncol=ncol(network), byrow=TRUE)
  cell_probs = (row_probs + col_probs)/2
  null_nodf_vals = c()
  for (i in 1:nsim)
  {
    null_network = matrix(NA, ncol=ncol(network), nrow=nrow(network))
    null_network[] <- rbinom(n=ncol(network)*nrow(network), size=1, prob = cell_probs)
    null_nodf = vegan::nestednodf(null_network)$statistic["NODF"]
    null_nodf_vals = append(null_nodf_vals, null_nodf)
  }
  nodf = vegan::nestednodf(network)$statistic["NODF"]
  res = list("NODF" = nodf, "relative_NODF" = (nodf - mean(null_nodf_vals)) / sd(null_nodf_vals))
  return(res)
}

get_network_features <- function(network_path, nsim=1000)
{
  network <- process_network(network_path)
  interaction_based_features = get_interaction_based_features(network)
  nestedness_features = get_nestedness(network, nsim)
  newtwork_features = append(interaction_based_features, nestedness_features)
  return(newtwork_features)
}

get_species_features <- function(network_path, nsim=1000)
{
  network <- process_network(network_path)
  species_features = data.frame(matrix(ncol = 4, nrow = nrow(network)))
  colnames(species_features) = c("network", "species", "rank", "nestedness_contribution")
  species_features$network = network_path
  species_features$species = rownames(network)
  species_features$rank = rank(rowSums(network), ties.method="average")-1
  species_features$rank = species_features$rank / max(species_features$rank)
  species_features$nestedness_contribution = unlist(nestedcontribution(network, nsimul=nsim)['lower level'], use.names=FALSE)
  return(species_features)
}

get_interaction_features <- function(network_path, nsim=1000)
{
  network <- process_network(network_path)
  interaction_features = data.frame(matrix(ncol = 6, nrow = nrow(network)*ncol(network)))
  colnames(interaction_features) = c("network", "plant", "pollinator", "plant_dependance_on_pollinator", "pollinator_dependance_on_plant", "preference")
  plant_dependence_values = network/rowSums(network)
  pollinator_dependence_values = t(t(network)/colSums(network))
  interaction_features$network = network_path
  interaction_features$plant = rep(rownames(network), ncol(network))
  interaction_features$pollinator = rep(colnames(network), each=nrow(network))
  interaction_features$plant_dependance_on_pollinator = as.vector(plant_dependence_values)
  interaction_features$pollinator_dependance_on_plant = as.vector(pollinator_dependence_values)
  interaction_features$preference = 
  return(interaction_features)
}