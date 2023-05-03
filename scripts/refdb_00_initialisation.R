#VARIABELEN
#-----------
db_name <- "refdb_2023-06-01"
user_name <- "pieter.verschelde@inbo.be"
docker_path <- "docker_image_obitools3"
docker_container_name <- "obitools3testcontainer" 
root_gdrive <- file.path("G:", 
                         ".shortcut-targets-by-id",
                         "0B1XJuciaZSENZG55ZnlDQ0FvT0E", 
                         "PRJ_eDNA",
                         "PRJ_eDNA_Refdb_2023" )
fasta_inputs <- file.path(root_gdrive, 'input_seqs', 'Riaz', 'import')
root_gdrive_key <- "16R6hn8ov_B3eFHJ1Obhyn8BxuhK-rZnP"
metadata_gdrive_key <- "1NidznDq9EVHN4_wfrhv0W5gZ0Q-S8wEn5jvXbb9rf8k"

taxdump_name <- "new_taxdump_2022-05-02"
fasta_name <- "input.fasta"

googlesheets4::gs4_auth(user_name)

#folders
#--------

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
