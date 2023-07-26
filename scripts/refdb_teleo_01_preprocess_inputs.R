
## GOOGLESHEET REFDB

df_species_all <- read_sheet(metadata_gdrive_key, "Soortenlijst") %>% 
  mutate(ROWNR = 1:n())

df_species <- df_species_all %>% 
  filter(Taxid > 0, !is.na(Priority))

df_priority <- df_species %>% 
  select(Taxid, Priority)

df_passlist <- read_sheet(metadata_gdrive_key, "Passlist_Teleo")

df_seq_errors <- read_sheet(metadata_gdrive_key, "Sequentiefouten")

df_multihit <-  read_sheet(metadata_gdrive_key, "Multihitlist_Teleo")

df_allowed_merges <- read_sheet(metadata_gdrive_key, "Toegelaten_merges")


### INPUTS

#get the source files
files <- sort(list.files(fasta_inputs, pattern = ".fasta"))

#read the source files
df_inputs_all <- NULL
for (file in files) {
  cat("\n\nINLEZEN VAN ", file, "\n--------------------------------\n")
  parsed <- parse_refdb_fasta(file.path(fasta_inputs, file))
  df_inputs_all <- df_inputs_all %>% 
    bind_rows(parsed)
}
whidup <- which(duplicated(df_inputs_all %>%  select(genlab_id, taxid, dna_sequence)))
df_inputs_all <- df_inputs_all[-whidup, ]
df_inputs_all$taxid <- trimws(df_inputs_all$taxid)

#find duplicate ids with different species
dup_ids <- df_inputs_all %>% group_by(genlab_id) %>% 
  summarise(aantal_seqs = n(), 
            sources = paste(source, collapse = ","),
            aantal_taxa = n_distinct(taxid), 
            taxa = paste(unique(taxid), collapse = ",")) %>% 
  filter(aantal_seqs > 1 | aantal_taxa > 1)
write_excel_csv2(dup_ids, paste0(output_path, '/', "duplicate_sequenties_verschillend_taxid.csv"))

## CURATE INPUTS

#let op: De sequentiefouten van Riaz worden hier al verwijderd in tegenstelling tot 2019
df_inputs <- df_inputs_all %>% 
  mutate(taxid = as.numeric(taxid)) %>% 
  filter(!(genlab_id %in% (df_seq_errors %>% pull(ENTRY_ID)))) %>% 
  filter(!taxid %in% (df_seq_errors %>% 
                        filter(TYPE == "taxid") %>% pull(TAXID)))

### remove alle sequenties die in multihit taxa staan en waar taxid != pref_taxid
df_inputs <- df_inputs %>% 
  filter(!taxid %in% (df_multihit %>% 
                        filter(TAXID != PREF_TAXID) %>% 
                        pull(TAXID)))

### add passlist sequences
df_inputs <- df_inputs %>% bind_rows(df_passlist) 

### remove full duplicates
df_inputs <- df_inputs %>% distinct(genlab_id, taxid, AMPLICON_HASH, 
                                    .keep_all = TRUE)

## create fasta input for ecopcr (obitools2)

#create_input_fasta(file = file.path("database", db_name, fasta_name), data = df_inputs)

#voor obitools3
create_input_fasta(file = file.path("database", db_name, fasta_name), 
                   lowercase = TRUE, 
                   data = df_inputs)

saveRDS(df_inputs, file=file.path(output_path, "df_inputs.RDS"))

