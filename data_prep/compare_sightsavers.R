##############################################
# Code author: erin r stearns
# Code objective: preparing dr. louise kelly-hope's literature-extracted entomological database to be a model input
# Date: 8.24.2022
#############################################

# PRE_REQs
#  - to run this script:
#        + make sure you change to your own sys env vars for setting data dirs
#        + define data_loc so makes sense in your environment         


rm(list = ls())

###################################################################################################################
# ----------------------------------------------- set up environment ---------------------------------------------
###################################################################################################################
#load packages
pacman::p_load(tidyverse,fs, readxl,
               ggplot2, viridis)

###################################################################################################################
# ----------------------------------------------- TO-DO ----------------------------------------------------------
###################################################################################################################
#set data directory
data_dir <- Sys.getenv("bmgf_proj_dir")

#LKH database
data_loc <- paste0(data_dir, "/Oncho_BlackFlyHabitatSuitability/data/")

###################################################################################################################
# ----------------------------------------------- define & load data ----------------------------------------------
###################################################################################################################
# sightsavers data
ssdata <-read.csv(paste0(data_loc, "BMGF NOEC FORMAT FOR COLLECTION OF BALCKFLY DATA FOR GOESPATIAL_Sightsavers_final.csv"),
                  stringsAsFactors = F)

# data from dr. makata
noecdata <- read.csv(paste0(data_loc, "BMGF NOEC FORMAT FOR COLLECTION OF BALCKFLY DATA FOR GOESPATIAL FINAL.csv"),
                     stringsAsFactors = F)

###################################################################################################################
# ----------------------------------------------- compare data ----------------------------------------------
###################################################################################################################
