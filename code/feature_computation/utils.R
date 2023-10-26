library(tidyverse)
library(bipartite)
library(matrixStats)


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

get_network_features <- function(network_path, null_sim_features_path, null_dir)
{
  network <- process_network(network_path)
  is_weighted = any(network != 1 & network != 0)
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

get_plant_nestedness_contribution <- function (web, nsim = 1000)
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
      nested.null <- sapply(1:nsim, function(x) {
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

# switch with group level - lower
get_species_features <- function(network_path, level="lower", nsim = 100)
{
  network <- process_network(network_path)
  network[is.na(network)] <- 0

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
    species_features["nestedness_contribution"] = get_plant_nestedness_contribution(web=network, nsim=nsim)
  } else {
    species_features["nestedness_contribution"] = get_pollinator_nestedness_contribution(web=network, nsim=nsim)
  }
  return(species_features)
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

# based on null model 2 of Bascompte 2003 (https://doi.org/10.1073/pnas.1633576100)
# generare_null_networks <- function(web, nsimul=100)
# {
#     web[web > 0] <- 1
# 	null_webs = replicate(nsimul, web, simplify = FALSE)
#     col_degs = colSums(web)/nrow(web)
#     row_degs = rowSums(web)/ncol(web)
#     for (s in 1:nsimul)
#     {
#       for (i in rownames(web))
#       {
#         for (j in colnames(web))
#         {
#           prob = 0.5*(col_degs[j]+row_degs[i])
#           null_webs[[s]][i, j] <- rbinom(1,1,prob)
#         }
#       }
#     }
# 	return (null_webs)
# }

generare_null_networks <- function(web, nsimul=100) # based on null model 2 of Bascompte 2003 (https://doi.org/10.1073/pnas.1633576100)
{
    web[web > 0] <- 1
	null_webs = replicate(nsimul, web, simplify = FALSE)
	for (i in rownames(web)) {
	  probs <- (rowSums(web)[i]/ncol(web) + colSums(web)/nrow(web))/2
	  for (s in 1:nsimul)
	  {
		null_webs[[s]][i, ] <- rbinom(ncol(web), 1, probs)
	  }
	}
	for (i in colnames(web)) {
	  probs <- (rowSums(web)/ncol(web) + colSums(web)[i]/nrow(web))/2
	  for (s in 1:nsimul)
	  {
		null_webs[[s]][, i] <- rbinom(nrow(web), 1, probs)
	  }
	}
	return (null_webs)
}

# implementation taken from Pires et al. 2020 https://doi.org/10.1002/ecy.3080
#Pires et al 2020
#adapted from Vieira and Almeida-Neto 2015; doi: 10.1111/ele.12394
#simulates extinction cascades by removing a target species
# The R arguments (R_rows, R_cols) can be a character string specifying an option (e.g. "normal", "VAlow", etc), or a numeric vector of length 1 or length == number of rows or columns. If a single number is given, this number is applied to all species in that guild.
# TargetGuild options: c("rows", "cols", "random_richness", "random_binary") (VA2015 used random_richness)
# TargetSpecies options: numeric (1:nrow or 1:ncol), "random_abundance", "random_binary"
#=========================================================================================================
netcascade_mod <- function(
					imatrix,
					R_val_rows,
					R_val_cols,
					TargetGuild,
					TargetSpecies,
					extinct_cols=NULL,
					extinct_rows=NULL,
					return.matrix =FALSE
					){

  #---------ARGUMENT CHECKS-----------------------------------

  if(class(imatrix)!="matrix" || (nrow(imatrix)+ncol(imatrix))<3){stop("
  	'imatrix' must be an object of class 'matrix', with at least three species
  	")}

   if(length(R_val_rows)!= nrow(imatrix)){stop("
   The length of vector'R_rows' must be equal to number of rows (i.e. species in guild) in 'imatrix'
   ")}

   if(length(R_val_cols)!= ncol(imatrix)){stop("
   The length of vector'R_cols' must be equal to number of columns in 'imatrix'
   ")}

  if((TargetGuild %in% c("rows","cols", "random_richness", "random_binary")) == FALSE){
	  stop('
  	Invalid target guild for primary extinction. Valid targets guilds are "rows", "cols", "random_richness", and "random_binary"
  	')}

  if(is.numeric(TargetSpecies)==FALSE && !(TargetSpecies %in% c("random_abundance", "random_binary"))){
	  stop('
  	Invalid value for the "TargetSpecies" argument. You may specify a single species by entering its row or column number or you may use a vector of relative probabilites for all species in the Target guild.
  	')}

  if(is.null(extinct_cols)==FALSE && class(extinct_cols)!= "integer"){stop("
  	extinct_cols must be either NULL an integer vector specifying the column numbers of species considered to be extinct on the original matrix
  	")}

  if(is.null(extinct_rows)==FALSE && class(extinct_rows)!= "integer"){stop("
  	extinct_rows must be either NULL or an integer vector specifying the row numbers of species considered to be extinct on the original matrix
  	")}

  #---------DEFINING SOME OBJECTS---------------------------
  nrows <- nrow(imatrix)
  ncols <- ncol(imatrix)
  cols <- 1:ncols
  rows <- 1:nrows
  colsNA <- 1:ncols
  rowsNA <- 1:nrows
  colsNA[extinct_cols] <- NA
  rowsNA[extinct_rows] <- NA

  degree_when_lost_cols <- c()
  degree_when_lost_rows <- c()



  #----------CALCULATING DEPENDENCE MATRICES-------------------
  interaction_strength <- array(0,dim=c(nrows,ncols,2))

  #matrix of columns' dependence on each row
  interaction_strength[,,1] <- t(t(imatrix)/rowSums(t(imatrix)))

  #matrix of rows' dependence on each column
  interaction_strength[,,2] <- imatrix/rowSums(imatrix)

  #-----------CHOOSING TARGET SPECIES FOR PRIMARY EXTINCTION---
  coext_rows <- c()
  coext_cols <- c()
  if(is.numeric(TargetSpecies) == TRUE){
	  if(length(TargetSpecies)==1){
		if(TargetGuild=="rows"){
		  if(TargetSpecies %in% extinct_rows){
			  stop('Specified target species for the primary extinction is already extinct')
		  }
		  coext_rows <- TargetSpecies
		  degree_when_lost_rows <- 1 #stores the degree of the extinction event of every row lost during the coextinction cascade.
		}
		if(TargetGuild=="cols"){
		  if(TargetSpecies %in% extinct_cols){
			  stop('Specified target species for the primary extinction is already extinct')
		  }
		  coext_cols <- TargetSpecies
		  degree_when_lost_cols <- 1
		}
	  }else{
		nspecies <- switch(TargetGuild, rows = nrows, cols = ncols)
		if(length(TargetSpecies)==nspecies){
		  if(TargetGuild =="rows"){
			alive_rows <- rows[is.na(rowsNA)==FALSE]
			coext_rows <- sample(c(alive_rows,0),1,prob = c(TargetSpecies[is.na(rowsNA)==FALSE],0))
			degree_when_lost_rows <- 1
		  }
		  if(TargetGuild =="cols"){
			alive_cols <- cols[is.na(colsNA)==FALSE]
			coext_cols <- sample(c(alive_cols,0),1,prob = c(TargetSpecies[is.na(colsNA)==FALSE],0))
			degree_when_lost_cols <- 1
		  }
		}else{
		  stop('Length of "TargetSpecies" must be 1 (specifying a single species within the Target guild) or else be equal to the number of species in the Target guild (specifying probabilities of primary extinction for each species in the Target guild)')
		}
	  }
  } else if(TargetSpecies == "random_abundance"){
  # TargetSpecies options: numeric (1:nrow or 1:ncol), "random_abundance", "random_binary"
	  if(TargetGuild=="rows"){
		alive_rows <- rows[is.na(rowsNA)==FALSE]
		coext_rows <- sample(
						x = alive_rows,
						size = 1,
						replace = FALSE,
						prob = rowSums(imatrix[alive_rows,])
			)
		degree_when_lost_rows <- 1 #stores the degree of the extinction event of every row lost during the coextinction cascade.
	  }
	  if(TargetGuild=="cols"){
		alive_cols <- cols[is.na(colsNA)==FALSE]
		coext_cols <- sample(
						x = alive_cols,
						size = 1,
						replace = FALSE,
						prob = colSums(imatrix[,alive_cols])
			)
		degree_when_lost_cols <- 1
	  }
  } else if(TargetSpecies == "random_binary"){
	  if(TargetGuild=="rows"){
		alive_rows <- rows[is.na(rowsNA)==FALSE]
		coext_rows <- sample(
						x = alive_rows,
						size = 1,
						replace = FALSE,
						prob = NULL
			)
		degree_when_lost_rows <- 1 #stores the degree of the extinction event of every row lost during the coextinction cascade.
	  }
	  if(TargetGuild=="cols"){
		alive_cols <- cols[is.na(colsNA)==FALSE]
		coext_cols <- sample(
						x = alive_cols,
						size = 1,
						replace = FALSE,
						prob = NULL
			)
		degree_when_lost_cols <- 1
	  }
  } else {
	  stop('Could not evaluate your non-numeric value for the "TargetSpecies" argument.')
  }


  imatrix[coext_rows,] <- 0
  imatrix[,coext_cols] <- 0

  lost_rows <- coext_rows #final list of rows which were "alive" in the original community but became extinct during this primary extinction + extinction cascade
  lost_cols <- coext_cols

  #-------------------CASCADE LOOP---------------------------
  equilibrium <- FALSE
  degree <- 1
  degree_table <- data.frame(degree,guild=factor(TargetGuild,levels=c("rows","cols")),n_extinctions=1)
  while(equilibrium == FALSE){
    extinct_rows <- coext_rows
    extinct_cols <- coext_cols
    colsNA[extinct_cols] <- NA
    rowsNA[extinct_rows] <- NA

    remaining_rows <- rows[is.na(rowsNA) == FALSE]
    remaining_cols <- cols[is.na(colsNA) == FALSE]

    # If one of the rows is extinct...
    if(length(extinct_rows)>0){
      for(i in 1:length(extinct_rows)){
        Target <- R_val_cols[remaining_cols]*interaction_strength[extinct_rows[i],remaining_cols,1] > runif(length(remaining_cols))
        coext_cols = c(coext_cols, remaining_cols[Target])
      }
      coext_rows <- c()
      coext_cols <- unique(coext_cols)
      colsNA[coext_cols] <- NA
      lost_cols <- c(lost_cols,coext_cols)

      # remove all of that column's interactions
      imatrix[,coext_cols] <- 0
      # make its interaction strengths zero; recalculate interaction strengths
      for(i in 1:ncols){
        if(sum(imatrix[,i])==0){
          interaction_strength[,i,1] <- 0
        }else{
          interaction_strength[,i,1] <- imatrix[,i]/sum(imatrix[,i])
        }
      }
      if(length(coext_cols)>0){
        degree <- degree + 1
        degree_when_lost_cols <- c(degree_when_lost_cols, rep(degree,length(coext_cols)))
        degree_table[degree,] <- data.frame(degree,"cols",length(coext_cols))
      }
    }else{
      for(i in 1:length(extinct_cols)){
        Target <- R_val_rows[remaining_rows]*interaction_strength[remaining_rows,extinct_cols[i],2] > runif(length(remaining_rows))
        coext_rows <- c(coext_rows, remaining_rows[Target])
      }
      coext_cols = c();
      coext_rows <- unique(coext_rows)
      lost_rows <- c(lost_rows,coext_rows)
      rowsNA[coext_rows] <- NA
      imatrix[coext_rows,] <- 0
      for(i in 1:nrows){
        if(sum(imatrix[i,])==0){
          interaction_strength[i,,2] <- 0
        }else{
          interaction_strength[i,,2] <- imatrix[i,]/sum(imatrix[i,])
        }
      }
      if(length(coext_rows)>0){
        degree <- degree + 1
        degree_when_lost_rows <- c(degree_when_lost_rows, rep(degree,length(coext_rows)))
        degree_table[degree,] <- data.frame(degree,"rows",length(coext_rows))
      }
    }
    equilibrium <- equilibrium + (length(coext_cols)+length(coext_rows))==0
  }

  #-------------------OUTPUT---------------------------

    if(length(lost_rows)>0){
      spp_data_rows <- data.frame(lost_rows = lost_rows, degree_of_extinction = degree_when_lost_rows)
    }else{
      spp_data_rows <- "No rows were lost"
    }
    if(length(lost_cols)>0){
      spp_data_cols <- data.frame(lost_col = lost_cols,degree_of_extinction=degree_when_lost_cols)
    }else{
      spp_data_cols <- "No columns were lost"
    }
  if(return.matrix==TRUE){
    return(
		list(
			interaction_matrix = imatrix,
			lost_rows = lost_rows,
			lost_cols = lost_cols,
			cascade_data = degree_table,
			rows_species_data = spp_data_rows,
			cols_species_data = spp_data_cols
			)
		)
  }else{
    return(
		list(
			cascade_data = degree_table,
			rows_species_data = spp_data_rows,
			cols_species_data = spp_data_cols
			)
		)
  }
}
