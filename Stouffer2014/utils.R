library(tidyverse)
library(bipartite)
library(matrixStats)

min_max_scale <- function(v1)
{
  v2 = (v1-min(v1))/(max(v1)-min(v1))
  return(v2)
}

process_network <- function(network_path)
{
  network <- read_csv(file = network_path, show_col_types = FALSE)
  network$is_exotic = NULL
  network = select(network, -1)
  network = network %>% column_to_rownames(., var = "plant_name")
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

get_plant_nestedness_contribution <- function (web, nsimul = 1000) 
{
  web <- ifelse(web > 0, 1, 0)
  if (is.null(rownames(web))) 
    rownames(web) <- paste0("L", seq.int(nrow(web)))
  lower.out <- data.frame(row.names = rownames(web))
  lower.out$nestedcontribution <- NA
  if (any(dim(web) < 2)) {
    warning("Your web is too small for a meaningful computation of nestedcontrrank (and probably other indices)!")
  }
  else {
    nested.orig <- vegan::nestednodf(web)$statistic["NODF"]
    for (i in rownames(web)) {
      message(i)
      probs <- (rowSums(web)[i]/ncol(web) + colSums(web)/nrow(web))/2
      nested.null <- sapply(1:nsimul, function(x) {
        web.null <- web
        web.null[i, ] <- rbinom(ncol(web), 1, probs)
        vegan::nestednodf(web.null)$statistic["NODF"]
      })
      lower.out[i, "nestedcontribution"] <- (nested.orig - 
                                               mean(nested.null))/sd(nested.null)
    }
  }
  return(lower.out)
}


get_pollinator_nestedness_contribution <- function (web, nsimul = 1000) 
{
  web <- ifelse(web > 0, 1, 0)
  if (is.null(colnames(web))) 
    colnames(web) <- paste0("H", seq.int(ncol(web)))
  higher.out <- data.frame(row.names = colnames(web))
  higher.out$nestedcontribution <- NA
  if (any(dim(web) < 2)) {
    warning("Your web is too small for a meaningful computation of nestedcontrrank (and probably other indices)!")
  }
  else {
    nested.orig <- vegan::nestednodf(web)$statistic["NODF"]
    for (i in colnames(web)) {
      message(i)
      probs <- (rowSums(web)/ncol(web) + colSums(web)[i]/nrow(web))/2
      nested.null <- sapply(1:nsimul, function(x) {
        web.null <- web
        web.null[, i] <- rbinom(nrow(web), 1, probs)
        vegan::nestednodf(web.null)$statistic["NODF"]
      })
      higher.out[i, "nestedcontribution"] <- (nested.orig - 
                                                mean(nested.null))/sd(nested.null)
    }
  }
  return(higher.out)
}

get_species_features <- function(network_path, nsim=1000)
{
  network <- process_network(network_path)
  species_features = data.frame(matrix(ncol = 4, nrow = nrow(network)))
  colnames(species_features) = c("network", "species", "ranked_degree", "ranked_nestedness_contribution")
  species_features$network = basename(network_path)
  species_features$species = rownames(network)
  species_features$ranked_degree = min_max_scale(rank(rowSums(network), ties.method="average"))
  species_features$ranked_nestedness_contribution = min_max_scale(rank(unlist(get_plant_nestedness_contribution(network, nsimul=nsim), use.names=FALSE)))
  return(species_features)
}

get_interaction_features <- function(network_path, nsim=1000)
{
  network <- process_network(network_path)
  interaction_features = data.frame(matrix(ncol = 7, nrow = nrow(network)*ncol(network)))
  colnames(interaction_features) = c("network", "plant", "pollinator", "plant_dependance_on_pollinator", "plant_interactions_num", "pollinator_dependance_on_plant", "pollinator_interactions_num")
  plant_dependence_values = network/rowSums(network)
  pollinator_dependence_values = t(t(network)/colSums(network))
  interaction_features$network = basename(network_path)
  interaction_features$plant = rep(rownames(network), ncol(network))
  interaction_features$pollinator = rep(colnames(network), each=nrow(network))
  interaction_features$plant_dependance_on_pollinator = as.vector(plant_dependence_values)
  interaction_features$plant_interactions_num = rep(rowSums(network), ncol(network))
  interaction_features$pollinator_dependance_on_plant = as.vector(pollinator_dependence_values)
  interaction_features$pollinator_interactions_num = rep(colSums(network), each=nrow(network))
  return(interaction_features)
}


get_exotic_tendency <- function(network, classification)
{
  network = as_tibble(network, rownames=NA)
  species = colnames(network)
  exotic_plant_species = c(classification %>% filter(is_exotic == TRUE) %>% pull(species))
  num_interactions = colSums(network) # vector of number per pollinator
  only_exotic_network = network %>% filter(row.names(network) %in% exotic_plant_species)
  num_interactions_with_exotic = colSums(only_exotic_network)
  empirical_exotic_tendencies = num_interactions_with_exotic/num_interactions
  return(empirical_exotic_tendencies)
}


get_pollinator_features <- function(network, plant_species_classification_path, nsim = 1000)
{
  plant_classification = read_csv(plant_species_classification_path, show_col_types = FALSE)
  pollinator_features = data.frame(matrix(ncol = 5, nrow = ncol(network)))
  colnames(pollinator_features) = c("network", "pollinator", "ranked_degree", "ranked_nestedness_contribution", "exotic_tendency")
  pollinator_features$network = basename(network_path)
  pollinator_features$pollinator = colnames(network)
  pollinator_features$ranked_degree = min_max_scale(rank(colSums(network), ties.method="average"))
  nestedness_contribution = get_pollinator_nestedness_contribution(network, nsimul=1000)
  pollinator_features$ranked_nestedness_contribution = min_max_scale(rank(unlist(nestedness_contribution, use.names=FALSE)))
  pollinator_features$exotic_tendency = get_exotic_tendency(network, plant_classification)
  return (pollinator_features)
}
