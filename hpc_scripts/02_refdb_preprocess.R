### ------------------------- ###
### ----- INPUT by USER ----- ###
### ------------------------- ###

PRIMER_NAME="teleo"

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

#!inputs inlezen (nadat de sequentiefouten eruit zijn)
df_inputs_cleaned <- read_rds(str_replace(cleaned_input_fasta, '.fasta', '.RDS'))
df_inputs_raw <- read_rds(str_replace(all_input_fasta, '.fasta', '.RDS'))

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

#!toegelaten merges op hoger niveau
df_allowed_merges = df_allowed_merges %>% rename_all(., .funs = tolower)

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

# Launch job script 03_refdb_ecopcr_riaz.slurm with -i IN_FASTA and -t taxdump.gz
