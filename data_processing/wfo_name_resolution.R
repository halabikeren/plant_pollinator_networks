#!/usr/bin/env Rscript
# DESCRIPTION -------------------------------------------------------------
#Author: Tyler Poppenwimer
#Date: Wed Mar 30 10:45:12 2022
#Description: Translates a set of names using the WFO package.  Afterwards it filters the names based on a set of criteria that are described below

#World Flora Online Resolution Steps
  #Step 1: Load set of names
  #Step 2: Conduct Name Resolution using WFO.Match
  #Step 3: Obtain ONE match for each submitted name using the WFO.one

#Filter Steps
  #Step 1: Change all NA's in Fuzzy.Dist to 100 (maximum value)
  #Step 2: Step 2: Extract good matches (based on fuzzy distance)
  #Step 3: Extract matches that maybe good (fuzzy distance greater than 2, but less than 100)
  #Step 4: Check the maybe good matches by checking if submitted names shows up in WFO matched name
  #Step 5: Combine the good data together and the bad data together
  #Step 6: Filter good matches by taxonomic rank (check the unique terms in the results to determine what you want to keep)
  #Step 7: Filter good matches by taxonomic status
  #Step 8: combine good and bad matches


# CLEANING ----------------------------------------------------------------
rm(list = ls())
cat("\014")


# LIBRARY LOADING ---------------------------------------------------------
library(WorldFlora)
library(stringr)
library("optparse")

WFO.remember()

#### ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~ ####
#### ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~ ####
#Step 1: Load set of names
 option_list = list(
  make_option(c("-f", "--input_path"), type="character", default=NULL,
              help="path to file with names to resolve", metavar="character"),
  make_option(c("-o", "--output_path"), type="character", default="name_resolution_world_flora_output.csv",
              help="output file name [default= %default]", metavar="character")
 );
 opt_parser = OptionParser(option_list=option_list);
 opt = parse_args(opt_parser);
 input_path = opt$input_path
 output_path = opt$output_path

#Load names to check
  names_df = data.frame(read.csv(input_path)[2])
  colnames(names_df)[1] <- "spec.name"
  names_df$spec.name <- tolower(names_df$spec.name)

  
# Conduct Name Resolution -------------------------------------------------
#Setp 2: Conduct Name Resolution using WFO.Match
  translation_data <- WFO.match(spec.data = names_df, WFO.data = WFO.data, Genus="New.Genus", Species="New.Species")

  
#Step 3: Obtain ONE match for each submitted name using the WFO.one
  translation_data <- WFO.one(translation_data)
  
  
  
#### ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~ ####
#### ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~ ####
# Filter the Resolution Results -------------------------------------------
#Step 1: Change all NA's in Fuzzy.Dist to 100 (maximum value)
  translation_data$Fuzzy.dist[is.na(translation_data$Fuzzy.dist)] = 100
  
  
#Step 2: Extract good matches (based on fuzzy distance)
  #Good matches whose fuzzy distances are 1 or 2
  good_match = translation_data[translation_data$Fuzzy.dist == 1 |
                                translation_data$Fuzzy.dist == 2, ]
  
  
#Step 3: Extract matches that maybe good (fuzzy distance greater than 2, but less than 100)
  maybe_good = translation_data[translation_data$Fuzzy.dist > 2 & translation_data$Fuzzy.dist < 100, ]
  
  
#Step 4: Check the maybe good matches by checking if submitted names shows up in WFO matched name
  #Set a holder for matches
  is_good = rep(NA, nrow(maybe_good))
  
  #go through each entry and check if the accepted name shows up in the given name
  for (i in 1:nrow(maybe_good)) 
  {
    
    #If there is at least 1 true, then there is a match
    is_good[i] = str_detect(maybe_good$spec.name.ORIG[i], tolower(maybe_good$scientificName[i]))
    
  }
  
  #Extract the maybe goods that are actually GOOD and the maybe goods that are actually BAD
  are_good = maybe_good[is_good, ]
  are_bad = maybe_good[!is_good, ]

  
#Step 5: Combine the good data together and the bad data together
  #Combine the good matches with those that we determined are actually GOOD
  good_match = rbind(good_match, 
                     are_good)
  
  #Combine the  data that had fuzzy dist value of 100 (no match) and those that are actually BAD
  bad_match = rbind(translation_data[translation_data$Fuzzy.dist == 100, ],
                    are_bad)
  
    #Replace the taxonRank, family, genus, specificEpithet, and scientificName to NA for the bad_match
    bad_match$taxonRank = NA
    bad_match$family = NA
    bad_match$genus = NA
    bad_match$specificEpithet = NA
    bad_match$scientificName = NA
    
    
#Step 6: Filter good matches by taxonomic rank (check the unique terms in the results to determine what you want to keep)
  #Remove data for those that DO NOT have the following rank
    #1) species
    #2) Species
    #3) SPECIES
    #4) VARIETY
    #5) SUBSPECIES
    if (!("taxonRank" %in% colnames(good_match)))
    {
      good_match$taxonRank = NA
    }
    good_match$taxonRank[is.na(good_match$taxonRank)] = "NAN"
    good_match[good_match$taxonRank != "species" &
               good_match$taxonRank != "Species" &
               good_match$taxonRank != "SPECIES" &
               good_match$taxonRank != "VARIETY" &
               good_match$taxonRank != "SUBSPECIES", 
                    c("taxonRank", "family", "genus", "specificEpithet", "scientificName")] = NA


#Step 7: Filter good matches by taxonomic status
  #Set all NA taxonomicStatus to "Unchecked"
    good_match$taxonomicStatus[is.na(good_match$taxonomicStatus)] = "Unchecked"

  #Remove the data for those that have a taxonomic status of Unchecked
  good_match[good_match$taxonomicStatus == "Unchecked", 
                  c("taxonRank", "family", "genus", "specificEpithet", "scientificName")] = NA
    
    
#Step 8: combine good and bad matches
  translation_data = rbind(good_match,
                           bad_match)

# Step 9: write translation data to csv
  write.csv(translation_data, output_path)
  
#CLEAN UP
  rm(are_bad, are_good, bad_match, good_match, 
     maybe_good, i, is_good)
  
  
#### ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~ ####
#### ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~ ####
