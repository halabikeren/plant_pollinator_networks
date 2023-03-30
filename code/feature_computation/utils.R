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
  if (! file.exists(network_path))
  {
    print(paste("network ", network_path, "does not exist"))
    return (NA)
  }
  network <- read_csv(file = network_path, show_col_types = FALSE)
  if ("Plant" %in% colnames(network))
  {
    network = network %>% column_to_rownames(., var = "Plant")
  }
  else if ("...1" %in% colnames(network))
  {
    network = network %>% column_to_rownames(., var="...1")
  }
  else
  {
    print("no relevant columns")
  }
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


get_feature_values <- function(networks_features, feature_name)
{
  feature_vals = list()
  for (i in 1:length((networks_features)))
    {
      feature_vals[[i]] = networks_features[[i]][feature_name]
  }
  feature_vals = unlist(feature_vals, use.names=FALSE)
  return(feature_vals)
}


get_modularity <- function(network)
{
	modularity = NA
    out <- tryCatch(
            {
            modularity = computeModules(network)@likelihood
            },
            error=function(cond)
            {
            message("failed to compute modularity due to error")
            message(cond)
            },
            finally=
            {
            message("done with modularity")
            }
    )
	return(modularity)
}

get_network_features <- function(network_path, null_sim_features_path, null_dir, nsim = 100)
{
  network <- process_network(network_path)
  is_weighted = any(network > 1)
  null_network_paths = list.files(null_dir)
  null_networks_features = list()
  null_networks_features_dfs = list()

  for (i in 1:length(null_network_paths))
  {
    path = null_network_paths[i]
    null_network = process_network(paste(null_dir, path, sep=""))
    null_network_features = networklevel(web=null_network, weighted=is_weighted)
    null_network_features[["modularity"]] = get_modularity(null_network)
    null_networks_features[[i]] = null_network_features
    null_network_features[["network_index"]] = i
    null_network_features[["observed_network"]] = basename(network_path)
    df = t(data.frame(Reduce(rbind, null_network_features)))
    colnames(df) = names(null_network_features)
    null_networks_features_dfs[[i]] = df

  }
  null_networks_features_df = do.call(rbind, null_networks_features_dfs)
  write.csv(null_networks_features_df, null_sim_features_path, row.names = TRUE)

  network_features = networklevel(web=network, weighted=is_weighted)

  network_features[["modularity"]] = get_modularity(network)
  extinction_features = simulate_network_extinction(network, nsim=nsim)
  network_features = c(network_features, extinction_features)
  non_transformable_feature_names = list("H2",
                                    "connectance",
                                    "weighted connectance",
                                    "links per species",
                                    "number of species",
                                    "number.of.species.HL",
                                    "number.of.species.LL",
                                    "mean number of shared partners",
                                    "partner diversity",
                                    "robustness_mean",
                                    "robustness_median",
                                    "robustness_min",
                                    "robustness_median",
                                    "robustness_var",
                                    "robustness_std")
  features_names = names(network_features)
  for (i in 1:length(features_names))
  {
    if ((!features_names[i] %in% non_transformable_feature_names) & (!is.na(network_features[features_names[i]])))
    {
          null_values = get_feature_values(networks_features=null_network_features)
          val = network_features[features_names[i]] - mean(null_values)
          network_features[paste("standardized_", features_names[i])] = val
    }
  }

  return(network_features)
}

get_plant_nestedness_contribution <- function (web, nsimul = 1000)
{
  web <- ifelse(web > 0, 1, 0)
  if (is.null(rownames(web)))
    rownames(web) <- paste0("L", seq.int(nrow(web)))
  lower.out <- data.frame(row.names = rownames(web)) # lower - rows - plants
  lower.out$nestedcontribution <- NA
  if (any(dim(web) < 2)) {
    warning("Your web is too small for a meaningful computation of nested contr rank (and probably other indices)!")
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

# switch with group level - higher
get_pollinator_nestedness_contribution <- function (web, nsimul = 1000)
{
  if (is.null(colnames(web)))
    colnames(web) <- paste0("H", seq.int(ncol(web)))
  higher.out <- data.frame(row.names = colnames(web)) # higher - columns, pollinators
  higher.out$nestedcontribution <- NA
  if (any(dim(web) < 2)) {
    warning("Your web is too small for a meaningful computation of nested contr rank (and probably other indices)!")
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


# switch with group level - lower
get_species_features <- function(network_path, null_sim_features_path, null_dir, level="lower", nsim = 100)
{
  network <- process_network(network_path)
  network[is.na(network)] <- 0
  is_weighted = any(network > 1)

  if (level == "lower")
  {
    row_names = rownames(network)
  } else {
    row_names = colnames(network)
  }

  null_network_paths = list.files(null_dir)
  print(length(null_network_paths))
  null_species_features = list()
  null_species_features_dfs = list()
  for (i in 1:length(null_network_paths))
  {
    null_network = process_network(paste(null_dir, null_network_paths[i], sep=""))
    null_features = data.frame(specieslevel(web=null_network, level=level))
    null_species_features[[i]] = null_features
    null_features[["network_index"]] = i
    df = data.frame(null_features)
    row.names(df) = row_names
    colnames(df) = names(null_features)
    null_species_features_dfs[[i]] = df
  }
  null_species_features_df = do.call(rbind, null_species_features_dfs)
  null_species_features_df["observed_network"] = basename(network_path)
  write.csv(null_species_features_df, null_sim_features_path, row.names = TRUE)

  species_features = data.frame(specieslevel(web=network, level=level))
  row.names(species_features) = row_names
  if (level == "lower")
  {
    species_features["nestedness_contribution"] = get_plant_nestedness_contribution(web=network, nsimul=nsim)
  } else {
    species_features["nestedness_contribution"] = get_pollinator_nestedness_contribution(web=network, nsimul=nsim)
  }

  # delta transform using null features
  non_transformable_feature_names = list("degree",
                                         "d",
                                         "normalized.degree",
                                         "nestedness_contribution")
  features_names = names(species_features)
  for (i in 1:length(features_names))
  {
    if ((!features_names[i] %in% non_transformable_feature_names) & (!is.na(species_features[features_names[i]])))
    {
      null_values = get_feature_values(networks_features=null_species_features, feature_name=features_names[i])
      val = species_features[features_names[i]] - mean(null_values)
      species_features[paste("standardized_", features_names[i])] = val
    }
  }
  return(species_features)
}


get_interaction_features <- function(network_path, nsim=1000)
{
  network <- process_network(network_path)
  interaction_features = data.frame(matrix(ncol = 7, nrow = nrow(network)*ncol(network)))
  colnames(interaction_features) = c("network", "plant", "pollinator", "plant_dependence_on_pollinator", "plant_interactions_num", "pollinator_dependence_on_plant", "pollinator_interactions_num")
  plant_dependence_values = network/rowSums(network)
  pollinator_dependence_values = t(t(network)/colSums(network))
  interaction_features$network = basename(network_path)
  interaction_features$plant = rep(rownames(network), ncol(network))
  interaction_features$pollinator = rep(colnames(network), each=nrow(network))
  interaction_features$plant_dependence_on_pollinator = as.vector(plant_dependence_values)
  interaction_features$plant_interactions_num = rep(rowSums(network), ncol(network))
  interaction_features$pollinator_dependence_on_plant = as.vector(pollinator_dependence_values)
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


simulate_species_extinction <- function(network, nsim=1000,  survival_threshold = 0.5, level = "lower")
{
  species = rownames(network)
  if (level == "higher")
  {
    species = colnames(network)
  }
  Rstats_by_sp = list()

  for (s in (1:length(species)))
  {
    Rvalues<-NULL
    triggertally<-NULL
    cascadelength<-NULL
    c=1
    trigmat<-matrix(nrow=nsim, ncol=ncol(network))
    j=0
    plantsurvivorsC<-matrix(nrow=nsim, ncol=ncol(network))
    exttimesC<-matrix(nrow=nsim, ncol=ncol(network))

    for(k in 1:nsim){
      mymat<-network                                                      # make copy of original matrix to work on
      PLANT<-colSums(mymat)                                               # save plant degrees
      POL<-rowSums(mymat)                                                 # save pollinator degrees
      survivors<-NULL                                                     # create survivors to save pollinator counts
      plantdeaths<-NULL                                                   # create survivors to save plant counts
      pastplantdeaths<-NULL

      ntriggerhat<-NULL
      i=0
      j=j+1
      ttally=0

      while(sum(colSums(mymat))>0){

        if(length(ntriggerhat)>0){
          ntriggerhat<-sample(ntriggerhat)
          cascadelength[c]<-length(ntriggerhat)
          c=c+1
          for (m in 1:length(ntriggerhat)){
            ntrigger<-sp
            mymat[,ntrigger]<-0
            pastplantdeaths<-c(pastplantdeaths,ntrigger)
            i=i+1
            v.ext<-which(((POL-rowSums(mymat))/POL)>=survival_threshold)                   # check pollinators
            mymat[v.ext,]<-0                                                               # make pollinators extinct
            survivors[i]<-length(which(rowSums(mymat)>0))
            trigmat[j,i]<-1
            plantsurvivorsC[j,i]<-length(which(colSums(mymat)==0))
            exttimesC[k,i]<-ntrigger
          }
        }
        else{
          rtrigger<-s # the first primary extinction is of the species of interest
          if (level == "lower")
          {
            mymat[rtrigger,]<-0
          } else {
            mymat[,rtrigger]<-0
          }



        }

        p.ext<-which(((PLANT-colSums(mymat))/PLANT)>=survival_threshold)
        ntriggerhat<-c(setdiff(names(p.ext),pastplantdeaths),setdiff(pastplantdeaths,names(p.ext)))

      }

      triggertally[k]<-ttally
      survivorsx<-c(nrow(mymat),survivors)
      Rvalues[k]<-sum(survivorsx)/(ncol(mymat)*nrow(mymat))
    }

    Rstats = c("mean"=mean(Rvalues),
               "median"=median(Rvalues),
               "min"=min(Rvalues),
               "max"=max(Rvalues),
               "var"=var(Rvalues),
               "std"=sd(Rvalues))
    Rstats_by_sp[[sp]] = Rstats
  }
  return(Rstats_by_sp)
}


simulate_network_extinction <- function(network, # columns correspond to plants and rows to pollinators
                                        nsim = 1000, # number of simulations
                                        survival_threshold = 0.5) # probabilty for rewiring (similar to all species)
{
  Rvalues<-NULL
  triggertally<-NULL
  cascadelength<-NULL
  c=1
  trigmat<-matrix(nrow=nsim, ncol=ncol(network))
  j=0
  plantsurvivorsC<-matrix(nrow=nsim, ncol=ncol(network))
  exttimesC<-matrix(nrow=nsim, ncol=ncol(network))

  for(k in 1:nsim){
    mymat<-network                                                      # make copy of original matrix to work on
    PLANT<-colSums(mymat)                                               # save plant degrees
    POL<-rowSums(mymat)                                                 # save pollinator degrees
    survivors<-NULL                                                     # create survivors to save pollinator counts
    pastplantdeaths<-NULL

    ntriggerhat<-NULL
    i=0
    j=j+1
    ttally=0

    while(sum(colSums(mymat))>0){

      if(length(ntriggerhat)>0){
        ntriggerhat<-sample(ntriggerhat)
        cascadelength[c]<-length(ntriggerhat)
        c=c+1
        for (m in 1:length(ntriggerhat)){
          ntrigger<-ntriggerhat[m]
          mymat[,ntrigger]<-0
          pastplantdeaths<-c(pastplantdeaths,ntrigger)
          i=i+1
          v.ext<-which(((POL-rowSums(mymat))/POL)>=survival_threshold)                   # check pollinators
          mymat[v.ext,]<-0                                                               # make pollinators extinct
          survivors[i]<-length(which(rowSums(mymat)>0))
          trigmat[j,i]<-1
          plantsurvivorsC[j,i]<-length(which(colSums(mymat)==0))
          exttimesC[k,i]<-ntrigger
        }
      }
      else{
        rtrigger<-sample((names(colSums(mymat)[(colSums(mymat))>0])),1)
        mymat[,rtrigger]<-0
        pastplantdeaths<-c(pastplantdeaths,rtrigger)
        i=i+1
        v.ext<-which(((POL-rowSums(mymat))/POL)>=survival_threshold)                   # check pollinators
        mymat[v.ext,]<-0                                                               # make pollinators extinct
        survivors[i]<-length(which(rowSums(mymat)>0))                                  # save number of survivors
        ttally<-ttally+1
        trigmat[j,i]<-0
        plantsurvivorsC[j,i]<-length(which(colSums(mymat)==0))
        exttimesC[k,i]<-rtrigger
      }

      p.ext<-which(((PLANT-colSums(mymat))/PLANT)>=survival_threshold)
      ntriggerhat<-c(setdiff(names(p.ext),pastplantdeaths),setdiff(pastplantdeaths,names(p.ext)))

    }

    triggertally[k]<-ttally
    survivorsx<-c(nrow(mymat),survivors)
    Rvalues[k]<-sum(survivorsx)/(ncol(mymat)*nrow(mymat))
  }
  Rstats = c("robustness_mean"=mean(Rvalues),
             "robustness_median"=median(Rvalues),
             "robustness_min"=min(Rvalues),
             "robustness_max"=max(Rvalues),
             "robustness_var"=var(Rvalues),
             "robustness_std"=sd(Rvalues))
  return(Rstats)
}