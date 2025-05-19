# Adapted From: https://github.com/inbo/prj_edna_refdb

# library(RSQLite)
library(tidyverse)
library(googlesheets4)

# Load source
VSC_DATA_dir = Sys.getenv("VSC_DATA")

### User Input variables
### Google Drive Authentication
###-------------------------
#google drive gebruikersnaam
user_name <- "nick.dillen@inbo.be"
googlesheets4::gs4_auth(user_name) #googlesheets4::gs4_deauth()

###-------------------------
# INPUT on HPC
## Could be made into cmd arguments
output_folder = "/staging/leuven/stg_00184/genetic_diversity/reference_databases/PRJ_eDNA_Refdb_2025"
taxdump_name <- file.path(output_folder, "taxonomy/2025-05-16-taxdump.tar.gz")

# INPUT on GDRIVE
#zelf te kiezen namen en locaties
## This is not really necessary here, maybe split teleo and riaz so we can update and execute seperatly
## Now we do everything x2
db_name_teleo <- "refdb_teleo_XXX"
db_name_riaz  <- "refdb_riaz_XXX"

### Globale variabelen
###----------------------
## HARDCODED
#shortcut naar google drive path van het project
## MyDrive + shortcut to PROJECTEN... in a MOUNTED folder Rclone

system("rclone version") # make sure rclone module is loaded in the pre-run scriplet! (module load rclone)
system("rclone mount gdrive: ~/G_MOUNT -vv --daemon")

root_gdrive <- file.path("~/G_MOUNT", 
                         "GDRIVE_shortcuts",
                         "PRJ_eDNA_Refdb_2023" )

fasta_inputs_location <- file.path(root_gdrive, 'input_seqs', 'import')

fasta_inputs_location_id = "1Q0OgyA_WSHGWdZjy3LvjAzP9MGiqBmCr"

#google key van de sheet waar de soortenlijst, multihit, ... bewaard wordt
metadata_gdrive_key <- "1NidznDq9EVHN4_wfrhv0W5gZ0Q-S8wEn5jvXbb9rf8k"

cleaned_input_fasta <- file.path("database", "input", "cleaned_inputs.fasta") 
all_input_fasta <- file.path("database", "input", "all_inputs.fasta") 
fasta_name <- 'input.fasta'

refdb_location_riaz  <- file.path("database", db_name_riaz)
refdb_location_teleo <- file.path("database", db_name_teleo)

output_path_riaz  <- file.path("database", db_name_riaz, "output")
output_path_teleo <- file.path("database", db_name_teleo, "output")

### Environment klaarzetten
###----------------------------
source(file.path(VSC_DATA_dir, "prj_edna_refdb","scripts/_functions_fasta.R"))
source(file.path(VSC_DATA_dir, "prj_edna_refdb","scripts/_functions_refdb.R"))
source(file.path(VSC_DATA_dir, "prj_edna_refdb","scripts/_functions_postprocessing.R"))

# Create output directories in target folder
setwd(output_folder)
if (!dir.exists("database")) dir.create("database")
if (!dir.exists(refdb_location_riaz)) dir.create(file.path("database", db_name_riaz))
if (!dir.exists(refdb_location_teleo)) dir.create(file.path("database", db_name_teleo))

### Docker variabelen
###---------------------

# docker_image_name <- "obitools3pv"
obi_script_riaz  <- "refdb_riaz.sh"
obi_script_teleo <- "refdb_teleo.sh"

# docker_win_riaz  <- "docker_obitools3_riaz"
# docker_win_teleo <- "docker_obitools3_teleo"

# docker_container_name <- "obitools3container" 
# docker_path           <- paste0(docker_container_name, ":/app/")
# docker_taxonomy       <- paste0(docker_path, "taxdump.tar.gz")
# docker_input_fasta    <- paste0(docker_path, "input.fasta")

# docker_script_riaz  <- paste0('/app/', obi_script_riaz)
# docker_script_teleo <- paste0('/app/', obi_script_teleo)

# pscommand_build_riaz  <- paste("docker", "build", "-t", 
#                                docker_image_name, paste0("./", docker_win_riaz, "/"))
# pscommand_build_teleo <- paste("docker", "build", "-t", 
#                                docker_image_name, paste0("./", docker_win_teleo, "/"))
# 
# pscommand_delete <- paste("docker", "rm", docker_container_name)
# pscommand_run <- paste("docker run --rm -i -d --name", 
#                        docker_container_name, docker_image_name) 
# # -i for interactive mode, -d for detached mode (run in background) 
