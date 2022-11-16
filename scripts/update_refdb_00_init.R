#VARIABELEN
#-----------

input_source <- "1NidznDq9EVHN4_wfrhv0W5gZ0Q-S8wEn5jvXbb9rf8k"
user_name <- "pieter.verschelde@inbo.be"
googlesheets4::gs4_auth(user_name)
db_name <- "refdb_2022-11-11"
if (!dir.exists("input")) 
  dir.create("input")
if (!dir.exists("database")) 
  dir.create("database") 
if (!dir.exists(file.path("database", db_name))) 
  dir.create(file.path("database", db_name))
if (!dir.exists("output"))
  dir.create("output")
if (!dir.exists(file.path("output", db_name)))
  dir.create(file.path("output", db_name))
taxdump_name <- "taxdump_2022-11-09"
fasta_name <- "input.fasta"
fasta_name_obi3 <- "input_obi3.fasta"

#drive_ls("https://drive.google.com/drive/folders/1l7Wrv0Np7tt1I6qTC6uyOwmFoTy22g8e")

### locatie van input fasta (drive_download werkt niet op dit moment)
fasta_input_folder <- "https://drive.google.com/drive/folders/1l7Wrv0Np7tt1I6qTC6uyOwmFoTy22g8e"

fasta_input_path <- 
   file.path("G:/.shortcut-targets-by-id/0B1XJuciaZSENZG55ZnlDQ0FvT0E", 
             "PRJ_eDNA/PRJ_eDNA_Refdb_2021",
            "input_seqs/Riaz")

fasta_input_path_local <- "input"

### Refdb name


if (!dir.exists("input")) dir.create("input")

#PAKKETTEN EN FUNCTIES
#-----------------------

library(tidyverse)
library(googledrive)
library(googlesheets4)
library(RSQLite)
library(refdb) #https://github.com/fkeck/refdb
source("scripts/_functions_fasta.R")
source("scripts/_functions_refdb.R")
