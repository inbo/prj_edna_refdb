
#R omgeving
#-----------
db_name <- "refdb_teleo_2023-09-01"
user_name <- "pieter.verschelde@inbo.be"
googlesheets4::gs4_auth(user_name) #googlesheets4::gs4_deauth()
library(tidyverse)
library(googledrive)
library(googlesheets4)
library(RSQLite)
library(refdb) #https://github.com/fkeck/refdb (nog niet echt gebruikt)
source("scripts/_functions_fasta.R")
source("scripts/_functions_refdb.R")
source("scripts/_functions_postprocessing.R")

#VARIABELEN
#-----------

#variabelen projectstructuur

refdb_location <- paste0("database/", db_name)
output_path <- paste0("database/", db_name, "/output")

#google drive
fasta_inputs <- "G:\\.shortcut-targets-by-id\\0B1XJuciaZSENZG55ZnlDQ0FvT0E\\PRJ_eDNA\\PRJ_eDNA_Refdb_2023\\input_seqs\\import"
root_gdrive_key <- "16R6hn8ov_B3eFHJ1Obhyn8BxuhK-rZnP"
metadata_gdrive_key <- "1NidznDq9EVHN4_wfrhv0W5gZ0Q-S8wEn5jvXbb9rf8k"

#files waarmee gewerkt wordt

taxdump_name <- "taxonomy/taxdump_2023-07-24.tar.gz"
obi_script <- "refdb_teleo.sh"
fasta_name <- "input.fasta"

#docker stuurvariabelen
docker_win <- "docker_obitools3_teleo"
docker_container_name <- "teleocontainer" 
image_name <- "obitools3_teleo"

#namen in docker container zelf
docker_path <- paste0(docker_container_name, ":/app/")
docker_script <- paste0('/app/', obi_script)
docker_taxonomy <- paste0(docker_path, "taxdump.tar.gz")
docker_input_fasta <- paste0(docker_path, "input.fasta")

#powershell commando's om docker aan te sturen
pscommand_build <- paste("docker", "build", "-t", image_name, paste0("./", docker_win, "/")) 
pscommand_delete <- paste("docker", "rm", docker_container_name)
pscommand_run <- paste("docker run --rm -i -d --name", 
                       docker_container_name, image_name)# -i for interactive mode, -d for detached mode (run in background) 



