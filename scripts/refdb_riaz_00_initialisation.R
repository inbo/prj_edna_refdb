
#R omgeving
#-----------
db_name <- "refdb_riaz_2023-09-01"
taxdump_name <- "taxonomy/taxdump_2023-07-24.tar.gz"
user_name <- "pieter.verschelde@inbo.be"
googlesheets4::gs4_auth(user_name) #googlesheets4::gs4_deauth()
library(tidyverse)
library(googledrive)
library(googlesheets4)
library(RSQLite)
library(refdb) #https://github.com/fkeck/refdb
source("scripts/_functions_fasta.R")
source("scripts/_functions_refdb.R")
source("scripts/_functions_postprocessing.R")

#VARIABELEN
#-----------

#variabelen projectstructuur

refdb_location <- paste0("database/", db_name)
output_path <- paste0("database/", db_name, "/output")
input_fasta <- paste0(refdb_location, "/", "input.fasta")
if (!dir.exists("input")) 
  dir.create("input")
if (!dir.exists("database")) 
  dir.create("database") 
if (!dir.exists(file.path("database", db_name))) {
  dir.create(file.path("database", db_name))  
  dir.create(file.path("database", db_name, "output"))
}



#verbindingen met googld drive
root_gdrive <- file.path("G:", 
                         ".shortcut-targets-by-id",
                         "0B1XJuciaZSENZG55ZnlDQ0FvT0E", 
                         "PRJ_eDNA",
                         "PRJ_eDNA_Refdb_2023" )
root_gdrive_key <- "16R6hn8ov_B3eFHJ1Obhyn8BxuhK-rZnP"
metadata_gdrive_key <- "1NidznDq9EVHN4_wfrhv0W5gZ0Q-S8wEn5jvXbb9rf8k"

#files waarmee gewerkt wordt
fasta_name <- "input.fasta"
fasta_inputs <- file.path(root_gdrive, 'input_seqs', 'import')


#docker stuurvariabelen
obi_script <- "refdb_riaz.sh"
docker_win <- "docker_obitools3_Riaz"
docker_container_name <- "obitools3container" 
image_name <- "obitools3pv"
docker_path <- paste0(docker_container_name, ":/app/")
docker_taxonomy <- paste0(docker_path, "taxdump.tar.gz")
docker_script <- paste0('/app/', obi_script)
docker_taxonomy <- paste0(docker_path, "taxdump.tar.gz")
docker_input_fasta <- paste0(docker_path, "input.fasta")

#powershell commando's om docker aan te sturen
pscommand_build <- paste("docker", "build", "-t", image_name, paste0("./", docker_win, "/")) 
pscommand_delete <- paste("docker", "rm", docker_container_name)
pscommand_run <- paste("docker run --rm -i -d --name", 
                       docker_container_name, image_name)# -i for interactive mode, -d for detached mode (run in background) 



