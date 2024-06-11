
###setup Maaslin2 packages etc into Rstudio:###
rm(list=ls())
getwd()
setwd("/Users/luke/Dropbox (Partners HealthCare)/WL-49-MS3-Luke/Fred/Fede-maaslin/Luke-troubleshooting/")

#dir.create("MaAsLin-pira") # Create a new directory
 # Change the current working directory 
getwd() #check if directory has been successfully changed

#Load MaAsLin2 package into the R environment
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("Maaslin2")

library(Maaslin2)
?Maaslin2

#Read in counts table
df_input_data = read.table("input.txt", header = TRUE, sep = "\t",
                           row.names = 1,
                           stringsAsFactors = FALSE)
df_input_data[1:5, 1:5]

#check metadata
lines <- readLines("metadata.txt")
df_input_metadata <- read.table("metadata.txt", header = TRUE, sep = "\t",
                                stringsAsFactors = FALSE, fill = TRUE)


#Read in metadata
df_input_metadata = read.table("metadata.txt", header = TRUE, sep = "\t",
                               row.names = 1,
                               stringsAsFactors = FALSE)
df_input_metadata[1:5, ]


#t_df_input_metadata <- t(df_input_metadata)
#rownames(df_input_metadata)

#####huttenhower commands with reference catagory#####

fit_data = Maaslin2(
  input_data = df_input_data ,
  input_metadata = df_input_metadata ,
  output = "Maaslin_L6_all_max" ,
  fixed_effects = c('age','edss','pirma'),
  reference = c('pirma', 0),
  random_effects = c('mergerepeats'),
  #max_significance = 0.572912628 ,
  standardize = FALSE,
  cores = 8)

