
#R omgeving
#-----------
user_name <- "pieter.verschelde@inbo.be"
googlesheets4::gs4_auth(user_name) #googlesheets4::gs4_deauth()
library(tidyverse)
library(googledrive)
library(googlesheets4)
library(RSQLite)
library(refdb) #https://github.com/fkeck/refdb
source("scripts/_functions_fasta.R")
source("scripts/_functions_refdb.R")

#STRUCTUUR
#---------

#workdir (basisdirectory in R)
  #database
    #refdb_name
      #input.fasta
      #outputs nieuwe refdb
  #taxonomy
    #new_taxdump_YYYY-MM-DD.tar.gz (https://ftp.ncbi.nih.gov/pub/taxonomy/new_taxdump/)
  #docker_obitools3
    #Dockerfile (bestandje die bepaalt hoe de docker image eruit ziet)
  #scripts
    #refdb_00_initialization (altijd laten lopen, configuratie R sessie)
    #refdb_01_preprocess_inputs (de inputs verzamelen tot 1 input.fasta)
    #refdb_02_docker_processing (de refdb laten genereren via dofker)
    #refdb_03_postprocessing (zoeken naar problemen in de gecreëerde refdb)
    #_functions_fasta (ondersteunende functies voor fasta parsing)
    #_functions_refdb (ondersteundende functies voor referentiedb)
    #_functions_problems (ondersteunende functies voor opsporen problemen in refdb)


#VARIABELEN
#-----------

#variabelen projectstructuur
db_name <- "refdb_2023-06-01"
refdb_location <- paste0("database/", db_name)
input_fasta <- paste0(refdb_location, "/", "input.fasta")
if (!dir.exists("input")) 
  dir.create("input")
if (!dir.exists("database")) 
  dir.create("database") 
if (!dir.exists(file.path("database", db_name))) 
  dir.create(file.path("database", db_name))


#verbindingen met googld drive
root_gdrive <- file.path("G:", 
                         ".shortcut-targets-by-id",
                         "0B1XJuciaZSENZG55ZnlDQ0FvT0E", 
                         "PRJ_eDNA",
                         "PRJ_eDNA_Refdb_2023" )
root_gdrive_key <- "16R6hn8ov_B3eFHJ1Obhyn8BxuhK-rZnP"
metadata_gdrive_key <- "1NidznDq9EVHN4_wfrhv0W5gZ0Q-S8wEn5jvXbb9rf8k"

#files waarmee gewerkt wordt
fasta_inputs <- file.path(root_gdrive, 'input_seqs', 'Riaz', 'import')
taxdump_name <- "taxonomy/new_taxdump_2023-05-02.tar.gz"
obi_script <- "generate_refdb.sh"

#docker stuurvariabelen
docker_win <- "docker_obitools3"
docker_container_name <- "obitools3container" 
image_name <- "obitools3pv"
container_name <- "obitools3refdbcontainer"

#namen in docker container zelf
docker_path <- paste0(container_name, ":/app/")
docker_script <- paste0('/app/', obi_script)
docker_taxonomy <- paste0(docker_path, "taxdump.tar.gz")
docker_input_fasta <- paste0(docker_path, "input.fasta")

#powershell commando's om docker aan te sturen
pscommand_build <- paste("docker", "build", "-t", image_name, paste0("./", docker_win, "/")) 
pscommand_delete <- paste("docker", "rm", container_name)
pscommand_run <- paste("docker run --rm -i -d --name", 
                       container_name, image_name)# -i for interactive mode, -d for detached mode (run in background) 



