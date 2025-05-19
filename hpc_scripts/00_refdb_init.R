###################
### REF DB INIT ###
###################

# Adapted From: https://github.com/inbo/prj_edna_refdb

## INPUT:
# - user_name
# - output_folder
# - taxdump_name
# - db_name_teleo
# - db_name_riaz
# - root_gdrive

##################

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
output_folder = "/staging/leuven/stg_00184/genetic_diversity/reference_databases/PRJ_eDNA_Refdb_2025_TEST"

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
system("ls ~/G_MOUNT") # should list your gdrive: content

# construct path to desired folder, USER SPECIFIC HARDCODE!
root_gdrive <- file.path("~/G_MOUNT", 
                         "GDRIVE_shortcuts",
                         "PRJ_eDNA_Refdb_2023" )

# folder structure ~ DMP
fasta_inputs_location <- file.path(root_gdrive, 'input_seqs', 'import')

# Not used?
# fasta_inputs_location_id = "1Q0OgyA_WSHGWdZjy3LvjAzP9MGiqBmCr"

#google key van de sheet waar de soortenlijst, multihit, ... bewaard wordt
metadata_gdrive_key <- "1NidznDq9EVHN4_wfrhv0W5gZ0Q-S8wEn5jvXbb9rf8k"

## HARDOCED OUTPUT NAMES
my_TS = format(Sys.time(), "%Y%m%d")

cleaned_input_fasta <- file.path("database", "input", paste0(my_TS,"_input_seqs_cleaned.fasta"))
all_input_fasta <- file.path("database", "input", paste0(my_TS,"_input_seqs.fasta"))

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

# DATABASE output
if (!dir.exists("database")) dir.create("database")
if (!dir.exists("database/input")) dir.create("database/input")
if (!dir.exists(refdb_location_riaz)) dir.create(file.path("database", db_name_riaz))
if (!dir.exists(refdb_location_teleo)) dir.create(file.path("database", db_name_teleo))
