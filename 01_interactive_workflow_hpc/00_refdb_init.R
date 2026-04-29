###################
### REF DB INIT ###
###################

# !RUN INTERACTIVELY! #
# In interactive Rstudio-server session on HPC
# See: https://ondemand.hpc.kuleuven.be/pun/sys/dashboard/batch_connect/sys/rstudio/session_contexts/new

# Adapted From: https://github.com/inbo/prj_edna_refdb

## total INPUT:
# - user_name
# - output_folder
# - root_gdrive

### ------------------------- ###
### ----- INPUT by USER ----- ###
### ------------------------- ###

# Working directory where the results will be stored/read
USER_OUTPUT_DIR = "/staging/leuven/stg_00184/genetic_diversity/reference_databases/12S-Riaz-Teleo/PRJ_eDNA_Refdb/PRJ_eDNA_Refdb_2026_apr/"

#google drive gebruikersnaam -> assume this is set up in git config, and is the same as work-email
USER_NAME <- system("git config --global user.email", intern = T)

# Mountpoint for PRJ_eDNA_Refdb:
RCLONE_REFDB_MOUNTPOINT="~/rclone_mnt/PRJ_eDNA_Refdb"

##################

library(tidyverse)
library(googlesheets4)

# Get VSC_DATA dir -> script source location
VSC_DATA_dir = Sys.getenv("VSC_DATA")

### User Input variables
### Google Drive Authentication
###-------------------------
googlesheets4::gs4_auth(USER_NAME) #googlesheets4::gs4_deauth()

### Globale variabelen
###----------------------
# Mount PRJ_eDNA_Refdb: on user defined mountpoint
system("rclone version") # make sure rclone module is loaded in the pre-run scriplet! (module load rclone)
system(paste0("rclone mount PRJ_eDNA_Refdb: ", RCLONE_REFDB_MOUNTPOINT, " -vv --daemon")) # loads some minutes
system(paste0("ls ", RCLONE_REFDB_MOUNTPOINT)) # should list your gdrive: content

# construct path to desired folder, USER SPECIFIC HARDCODE!
root_gdrive <- file.path(RCLONE_REFDB_MOUNTPOINT,
                         "PRJ_eDNA_Refdb_Water" )

# folder structure ~ DMP
fasta_inputs_location <- file.path(root_gdrive, 'input_data' ,'input_seqs', 'import')

#google key van de sheet waar de soortenlijst, multihit, ... bewaard wordt
metadata_gdrive_key <- "1NidznDq9EVHN4_wfrhv0W5gZ0Q-S8wEn5jvXbb9rf8k"

### Environment klaarzetten
###----------------------------
source(file.path(VSC_DATA_dir, "prj_edna_refdb","r_functions/_functions_fasta.R"))
source(file.path(VSC_DATA_dir, "prj_edna_refdb","r_functions/_functions_refdb.R"))
source(file.path(VSC_DATA_dir, "prj_edna_refdb","r_functions/_functions_postprocessing.R"))

# Create output directories in target folder
dir.create(USER_OUTPUT_DIR)
setwd(USER_OUTPUT_DIR)
