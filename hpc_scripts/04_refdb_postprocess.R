### ------------------------- ###
### ----- INPUT by USER ----- ###
### ------------------------- ###
getwd()

PRIMER_NAME="riaz"
# PRIMER_NAME="teleo"

# output of the job script 03_refdb_ecopcr.slurm
# Most likely its in this list:
grep(".obidms$", x = list.dirs(".", full.names = T, recursive = T), value = T)

detected_ecopcr_files = list.files(".", full.names = T, recursive = T, pattern = "ecopcr.fasta")
OBITOOLS_OUTPUT = dirname(detected_ecopcr_files[max(grep(pattern = PRIMER_NAME, detected_ecopcr_files))])
# OBITOOLS_OUTPUT=file.path("./database/20251021_refdb_riaz/2025-10-23-obitools3-refdb-riaz/")

# Assert that inputs are still defined! These INPUT vars are assumed to exist
stopifnot(file.exists(cleaned_input_fasta),
          exists("metadata_gdrive_key"))

# ---------------------------
### HARDOCED OUTPUT NAMES ###

refdb_location  <- file.path(OBITOOLS_OUTPUT)

obi_input_fasta_path = list.files(refdb_location, pattern = "-input.fasta", full.names = T)
obi_ecopcr_fasta_path = list.files(refdb_location, pattern = "-ecopcr.fasta", full.names = T)
obi_refdb_fasta_path = list.files(refdb_location, pattern = "-ecopcr_final_0.995.fasta", full.names = T)

stopifnot(file.exists(obi_input_fasta_path),
          file.exists(obi_ecopcr_fasta_path),
          file.exists(obi_refdb_fasta_path))

##################
### READ INPUT ###
##################

# INPUT from GSHEET
#! allowed merges per primer
if (PRIMER_NAME == "riaz"){
  df_allowed_merges <- read_sheet(metadata_gdrive_key, "Toegelaten_merges_Riaz")
  df_multihit <- read_sheet(metadata_gdrive_key, "Multihitlist_Riaz", range = "A:H", 
                            col_types = 'ccccccnn')
} else if (PRIMER_NAME == "teleo"){
  df_allowed_merges <- read_sheet(metadata_gdrive_key, "Toegelaten_merges_Teleo")
  df_multihit <- read_sheet(metadata_gdrive_key, "Multihitlist_Teleo", range = "A:H", 
                            col_types = 'ccccccnn')
}

df_multihit = df_multihit %>% rename_all(., .funs = tolower)
df_allowed_merges = df_allowed_merges %>% rename_all(., .funs = tolower)

df_soortenlijst <- read_sheet(metadata_gdrive_key, "Soortenlijst") %>% 
  mutate(taxid = Taxid, priority = Priority, rank = Rank) %>% 
  filter(priority != 9)

# INPUT form OBITools3 DMS
obitools_input_data <- parse_refdb_fasta(obi_input_fasta_path)

ecopcr_data <- parse_refdb_fasta(obi_ecopcr_fasta_path)

merged_data <- parse_refdb_fasta(obi_refdb_fasta_path)  %>% 
  mutate(obi_taxid = str_replace(substring(lca_taxid, 2), "]", ""),
         taxid = str_replace(taxid, ";", ""))

# INPUT from pre-processed reference sequences GDrive
input_data_before_multihit <- read_rds(str_replace(cleaned_input_fasta, '.fasta', '.RDS'))

######################
### POSTPROCESSING ###
######################

ecopcr_combined <- combine_ecpocr_with_merged(ecopcr_data, merged_data) 


df_conflicts <- genereer_conflicten(ecopcr_combined, df_soortenlijst,
                                    df_allowed_merges)

df_soortenevaluatie <- genereer_soortenevaluatie(ecopcr_combined, 
                                                 merged_data,
                                                 input_data_before_multihit,
                                                 df_soortenlijst, 
                                                 df_multihit, 
                                                 df_allowed_merges, 
                                                 df_conflicts)
### ---------------
### write output LOCAL ###
write_excel_csv2(df_conflicts, 
                 file = file.path(refdb_location,"niet_op_soort_gebracht.csv"))

write_excel_csv2(df_soortenevaluatie, 
                 file = file.path(refdb_location, paste0("soortenevaluatie_", PRIMER_NAME,".csv")))


### write output GDrive ###
googlesheets4::write_sheet(df_soortenevaluatie, 
                           sheet = basename(OBITOOLS_OUTPUT),
                           ss = "1rS0ZsP4AO1wrgCnDPmAfTUbxlzmHNtR9bA7rqgVI2xE")

googlesheets4::write_sheet(df_conflicts, 
                           sheet = basename(OBITOOLS_OUTPUT),
                           ss = "1kTHcNz0HktqbFfVy8atN6lR5MvCCQxNUSrLQFph9pVY")

