# Adapted From: https://github.com/inbo/prj_edna_refdb

# Load source
VSC_DATA_dir = Sys.getenv("VSC_DATA")

### User Input variables
### Google Drive Authentication
###-------------------------
#google drive gebruikersnaam
user_name <- "nick.dillen@inbo.be"

#shortcut naar google drive path van het project
root_gdrive <- file.path("G:", 
                         ".shortcut-targets-by-id",
                         "0B1XJuciaZSENZG55ZnlDQ0FvT0E", 
                         "PRJ_eDNA",
                         "PRJ_eDNA_Refdb_2023" )
fasta_inputs_location <- file.path(root_gdrive, 'input_seqs', 'import')

#google key van de sheet waar de soortenlijst, multihit, ... bewaard wordt
metadata_gdrive_key <- "1NidznDq9EVHN4_wfrhv0W5gZ0Q-S8wEn5jvXbb9rf8k"

googlesheets4::gs4_auth(user_name) #googlesheets4::gs4_deauth()

###-------------------------

# INPUT
data_folder = "/staging/leuven/stg_00184/genetic_diversity/reference_databases/PRJ_eDNA_Refdb_2025"
taxdump_name <- file.path(data_folder, "taxonomy/taxdump_2023-07-24.tar.gz")
fasta_inputs_location <- file.path(data_folder, 'input_seqs', 'import')

setwd(data_folder)

#zelf te kiezen namen en locaties
db_name_teleo <- "refdb_teleo_XXX"
db_name_riaz  <- "refdb_riaz_XXX"

### Globale variabelen
###----------------------

cleaned_input_fasta <- file.path("database", "input", "cleaned_inputs.fasta") 
all_input_fasta <- file.path("database", "input", "all_inputs.fasta") 
fasta_name <- 'input.fasta'

refdb_location_riaz  <- file.path("database", db_name_riaz)
refdb_location_teleo <- file.path("database", db_name_teleo)

output_path_riaz  <- file.path("database", db_name_riaz, "output")
output_path_teleo <- file.path("database", db_name_teleo, "output")

### Environment klaarzetten
###----------------------------

library(tidyverse)
# library(googlesheets4)
# library(RSQLite)
source(file.path(VSC_DATA_dir, "prj_edna_refdb","scripts/_functions_fasta.R"))
source(file.path(VSC_DATA_dir, "prj_edna_refdb","scripts/_functions_refdb.R"))
source(file.path(VSC_DATA_dir, "prj_edna_refdb","scripts/_functions_postprocessing.R"))


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
