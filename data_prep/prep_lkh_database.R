##############################################
# Code author: erin r stearns
# Code objective: preparing dr. louise kelly-hope's literature-extracted entomological database to be a model input
# Date: 2.25.2022
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
data_loc <- paste0(data_dir, "/Oncho_BlackFlyHabitatSuitability/data/LKH_NGA_LitExtraction/")


###################################################################################################################
# ----------------------------------------------- define & load data ------------------------------------------------------
###################################################################################################################
#config file
config <- read.csv(paste0(data_loc, "config_db_prep.csv"), stringsAsFactors = F)

#lkh database path
#db.path <- paste0(data_loc, "22_2.2.16_Nigeria_database.xlsx")
#updated db
db.path <- paste0(data_loc, "22_5.4_Blackfly historical_Database.xlsx")

###################################################################################################################
# ----------------------------------------------- prep data ------------------------------------------------------
###################################################################################################################
#make vector of sheet names
sheet.names <- db.path %>%
  excel_sheets()

# ---- concatenate all sheets, if relevant (new format from LKH is single tab) and transform into long format ----
if(length(sheet.names) == 1){
  #save sheet name
  df.name <- sheet.names[1]
  
  #get that data loaded in
  df_main <- read_excel(db.path, 1)
  
  # create col names sans leading/trailing white space
  new.names <- trimws(names(df_main))
  
  df_main <- df_main %>%
    rename_all(list(~new.names))
  
  #convert all cols to character bc will not row bind otherwise
  df_main[] <- lapply(df_main, as.character)  
  
} else { 
  for (i in 1:(length(sheet.names) - 1)) { #index ends at 2nd to last sheet bc final sheet is list of refs
    
    #save sheet name
    df.name <- sheet.names[i]
    
    #get that data loaded in
    df.temp <- read_excel(db.path, i) %>%
      #use tab name to create new field called pub_decade
      mutate(pub_decade = df.name) 
    
    # create col names sans leading/trailing white space
    new.names <- trimws(names(df.temp))
    
    df.temp <- df.temp %>%
      rename_all(list(~new.names))
    
    #convert all cols to character bc will not row bind otherwise
    df.temp[] <- lapply(df.temp, as.character)  
    
    #if first sheet, use to initiate main table
    if(i == 1) {
      df_main <- df.temp 
      
    } else {
      df_main <- df_main %>%
        bind_rows(df.temp) 
    }
  }
  
  }


# ---- subset to fields of interest --
#prepare vector of fields of interest
foi <- config %>%
  filter(include == 1) %>%
  select(model_input_field) #db_field

foi <- trimws(foi[["model_input_field"]]) #"db_field"

#prepare vector of new field names
new_foi <- config %>%
  filter(include == 1) %>%
  select(model_input_field)
  
new_foi <- trimws(new_foi[["model_input_field"]])

#subset data
df_main_sub <- df_main %>%
  select((foi)) #%>%
  #rename_all(list(~new_foi))

###################################################################################################################
# ----------------------------------------------- save ------------------------------------------------------
###################################################################################################################
write.csv(df_main_sub, file = paste0(data_loc, "22_5.4_nga_lit_extr.csv"))
