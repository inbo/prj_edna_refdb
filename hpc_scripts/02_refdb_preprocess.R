### ------------------------- ###
### ----- INPUT by USER ----- ###
### ------------------------- ###

PRIMER_NAME="riaz"
# PRIMER_NAME="teleo"

# Assert that inputs are still defined! These INPUT vars are assumed to exist
stopifnot(file.exists(cleaned_input_fasta),
          file.exists(all_input_fasta),
          exists("metadata_gdrive_key"))

# ---------------------------
### HARDOCED OUTPUT NAMES ###

my_TS = format(Sys.time(), "%Y%m%d")
# set names
db_name <- paste0( my_TS, "_refdb_", PRIMER_NAME)
refdb_location  <- file.path("database", db_name)
output_path <- file.path(refdb_location, "output")

# output FASTA name base on input
refseq_ID = stringr::str_split(basename(cleaned_input_fasta), pattern="_ref")[[1]][1]
fasta_name = file.path(refdb_location, paste0(refseq_ID, "_", PRIMER_NAME, "_input.fasta"))

##############################
### INPUTS OP ORDE STELLEN ### 
##############################

#!inputs inlezen uit vorige stap 01_refdb_import.R (nadat de sequentiefouten eruit zijn)
df_inputs_cleaned <- read_rds(str_replace(cleaned_input_fasta, '.fasta', '.RDS'))
df_inputs_raw <- read_rds(str_replace(all_input_fasta, '.fasta', '.RDS'))

# Read input Google Sheet
#! opnieuw toegestaan
if (PRIMER_NAME == "riaz"){
  df_passlist <- read_sheet(metadata_gdrive_key, "Passlist_Riaz")
  df_multihit <- read_sheet(metadata_gdrive_key, "Multihitlist_Riaz", range = "A:H", 
                             col_types = 'ccccccnn')
  df_allowed_merges <- read_sheet(metadata_gdrive_key, "Toegelaten_merges_Riaz")
} else if (PRIMER_NAME == "teleo"){
  df_passlist <- read_sheet(metadata_gdrive_key, "Passlist_Teleo")
  df_multihit <- read_sheet(metadata_gdrive_key, "Multihitlist_Teleo", range = "A:H", 
                             col_types = 'ccccccnn')
  df_allowed_merges <- read_sheet(metadata_gdrive_key, "Toegelaten_merges_Teleo")
}

df_passlist <- df_inputs_raw %>% 
  filter(genbank_id %in% df_passlist$ENTRY_ID) %>% 
  mutate(taxid = as.numeric(taxid))
nrow(df_passlist)

#! multihitsoorten
df_multihit = df_multihit %>% rename_all(., .funs = tolower)

#! Pas de multihitlijst toe
#-----------------------------

df_inputs <- df_inputs_cleaned %>% 
  filter(!taxid %in% (df_multihit %>% 
                        filter(taxid != pref_taxid) %>% 
                        pull(taxid)))

#!  Herintroduceer sequenties uit de passlist
#---------------------------

df_inputs <- df_inputs %>% bind_rows(df_passlist) 

# remove full duplicates
df_inputs <- df_inputs %>% distinct(genbank_id, taxid, dna_hash, 
                                    .keep_all = TRUE)

## create fasta input for ecopcr (obitools2)
#---------------------------------------------

# Create output dir
dir.create(refdb_location)

# Write FASTA file
create_input_fasta(file = fasta_name, 
                   lowercase = TRUE, 
                   data = df_inputs)

#################################
### GENEREREN DB
#################################

# Launch job script
# Example for riaz:
# sbatch $VSC_DATA/prj_edna_refdb/hpc_scripts/03_refdb_ecopcr.slurm \
#   -i 20250520_PRJ_eDNA_Refdb_2023_riaz_input.fasta \
#   -t ../../taxonomy/2025-05-16-taxdump.tar.gz \
#   -p riaz
cat(paste0("In your HPC-terminal:\n\n#Go to the directory with input sequences\ncd ", USER_OUTPUT_DIR, "/",refdb_location,
           "\n\n# Launch the OBITools3-ecopcr JOB to perform in-silico PCR and create a REFDB :\n",
  "sbatch",
  " $VSC_DATA/prj_edna_refdb/hpc_scripts/03_refdb_ecopcr.slurm",
             " -i ", basename(fasta_name),
             " -t ", "../../taxonomy/taxdump.tar.gz",
             " -p ", PRIMER_NAME, "\n\nWait untill the job has completed before you continue!\nYou can track this with:\n\nsqueue -M all"))
