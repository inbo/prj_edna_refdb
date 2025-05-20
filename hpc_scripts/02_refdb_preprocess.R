##############################
### INPUTS OP ORDE STELLEN ### 
##############################

#!inputs inlezen (nadat de sequentiefouten eruit zijn)
df_inputs_cleaned <- read_rds(str_replace(cleaned_input_fasta, '.fasta', '.RDS'))
df_inputs_raw <- read_rds(str_replace(all_input_fasta, '.fasta', '.RDS'))

#! opnieuw toegestaan
df_passlist <- read_sheet(metadata_gdrive_key, "Passlist_Riaz")
df_passlist <- df_inputs_raw %>% 
  filter(genbank_id %in% df_passlist$ENTRY_ID) %>% 
  mutate(taxid = as.numeric(taxid))
nrow(df_passlist)


#! multihitsoorten
df_multihit <-  read_sheet(metadata_gdrive_key, "Multihitlist_Riaz", range = "A:H", 
                           col_types = 'ccccccnn') %>% 
  rename_all(., .funs = tolower)

#!toegelaten merges op hoger niveau
df_allowed_merges <- read_sheet(metadata_gdrive_key, "Toegelaten_merges_Riaz") %>% 
  rename_all(., .funs = tolower)

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

create_input_fasta(file = file.path("database", db_name_riaz, "input.fasta"), 
                   lowercase = TRUE, 
                   data = df_inputs)

#################################
### GENEREREN DB
#################################

# Launch job script 03_refdb_ecopcr_riaz.slurm with -i IN_FASTA and -t taxdump.gz
