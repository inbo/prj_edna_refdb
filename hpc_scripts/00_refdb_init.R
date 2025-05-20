###################
### REF DB INIT ###
###################

# RUN INTERACTIVELY #

# Adapted From: https://github.com/inbo/prj_edna_refdb

## total INPUT:
# - user_name
# - output_folder
# - taxdump_name
# - db_name_teleo
# - db_name_riaz
# - root_gdrive

### ------------------------- ###
### ----- INPUT by USER ----- ###
### ------------------------- ###

#directory where the results will be stored
USER_OUTPUT_DIR = "/staging/leuven/stg_00184/genetic_diversity/reference_databases/PRJ_eDNA_Refdb_2025_TEST"

#google drive gebruikersnaam
USER_NAME <- "nick.dillen@inbo.be"

##################

library(tidyverse)
library(googlesheets4)

# Get VSC_DATA dir -> script source location
VSC_DATA_dir = Sys.getenv("VSC_DATA")

### User Input variables
### Google Drive Authentication
###-------------------------
googlesheets4::gs4_auth(USER_NAME) #googlesheets4::gs4_deauth()

###-------------------------
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

#system("module load rclone") # Load Rclone module on HPC -THIS DOES NOT WORK. NEED TO ADD TO Rstudio server FORM
system("rclone version") # make sure rclone module is loaded in the pre-run scriplet! (module load rclone)
system("rclone mount gdrive: ~/G_MOUNT -vv --daemon") # loads some minutes
system("ls ~/G_MOUNT") # should list your gdrive: content

# construct path to desired folder, USER SPECIFIC HARDCODE!
root_gdrive <- file.path("~/G_MOUNT", 
                         "GDRIVE_shortcuts",
                         "PRJ_eDNA_Refdb_2023" )

# folder structure ~ DMP
fasta_inputs_location <- file.path(root_gdrive, 'input_seqs', 'import')

#google key van de sheet waar de soortenlijst, multihit, ... bewaard wordt
metadata_gdrive_key <- "1NidznDq9EVHN4_wfrhv0W5gZ0Q-S8wEn5jvXbb9rf8k"

## HARDOCED OUTPUT NAMES
my_TS = format(Sys.time(), "%Y%m%d")

# For processing of INPUT from GDrive
import_dir_name = paste0("input_reference_sequences")
all_input_fasta <- file.path("database", import_dir_name, paste0(my_TS, "_", basename(root_gdrive),"_reference_sequences.fasta"))
cleaned_input_fasta <- file.path("database", import_dir_name, paste0(my_TS,"_", basename(root_gdrive),"_reference_sequences_cleaned.fasta"))

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
setwd(USER_OUTPUT_DIR)

# DATABASE output
dir.create(file.path("database", import_dir_name), recursive = T)
dir.create(file.path("database", db_name_riaz))
dir.create(file.path("database", db_name_teleo))
