###############
#### INPUT ####
###############
# Adapted by Nick Dillen

#! soortenlijst
df_species_all <- read_sheet(metadata_gdrive_key, "Soortenlijst") %>% 
  mutate(ROWNR = 1:n())

df_species <- df_species_all %>% 
  filter(Taxid > 0, !is.na(Priority), Priority != 9)

#! prioriteiten
df_priority <- df_species %>% 
  select(Taxid, Priority)

#! input sequenties
## Make sure to have a working rclone mount-point!
print(fasta_inputs_location)

input_fasta_files <- sort(list.files(fasta_inputs_location, pattern = ".fasta", full.names = T))
df_inputs_orig <- NULL

for (f in input_fasta_files) {
  cat("\n\nINLEZEN VAN ", f, "\n--------------------------------\n")
  f_parsed <- parse_refdb_fasta(f)
  df_inputs_orig <- df_inputs_orig %>% 
    bind_rows(f_parsed)
}

#! sequentiefouten
df_seq_errors <- read_sheet(metadata_gdrive_key, "Sequentiefouten")


#! nakijken of er duplicaten zijn
#----------------------------------
# Clean whitespace (why 2 =/= functions?)
df_inputs_all <- df_inputs_orig %>% 
  mutate(genbank_id = stringr::str_trim(genbank_id),
         taxid = trimws(df_inputs_orig$taxid))

# Select duplicated records
whidup <- which(duplicated(df_inputs_all %>%  select(genbank_id)))
dupids <- df_inputs_all %>% slice(whidup) %>% pull(genbank_id)
#df_inputs_all %>% filter(genbank_id %in% dupids) %>% arrange(genbank_id) %>% view()

# remove (dangerous in interactive, run exactly ONCE!)
if (length(whidup)) df_inputs_all <- df_inputs_all[-whidup, ]

saveRDS(df_inputs_all, file.path(str_replace(all_input_fasta, ".fasta", '.RDS' )))
create_input_fasta(file = all_input_fasta, 
                   lowercase = TRUE, 
                   data = df_inputs_all)


#Kijk of er nog prioriteit-9 (soorten die niet weerhouden worden wegens bvb hybriden) aan seq_errors moet toegevoegd worden

#! Haal foutieve sequenties uit 
#---------------------------------

df_inputs <- df_inputs_all %>% 
  mutate(taxid = as.numeric(taxid)) %>% 
  filter(!(genbank_id %in% (df_seq_errors %>% pull(ENTRY_ID)))) %>% 
  filter(!taxid %in% (df_seq_errors %>% 
                        filter(TYPE == "taxid") %>% pull(TAXID)))

#check of geen prioriteit-9 soorten overblijven (anders toevoegen aan foutieve sequenties)
df_inputs %>% filter(taxid %in% (df_species %>% filter(Priority == 9) %>% pull(Taxid)))

#find duplicate ids with different species
dup_ids <- df_inputs %>% group_by(genbank_id) %>% 
  summarise(aantal_seqs = n(), 
            sources = paste(source, collapse = ","),
            aantal_taxa = n_distinct(taxid), 
            taxa = paste(unique(taxid), collapse = ",")) %>% 
  filter(aantal_seqs > 1 & aantal_taxa > 1)

#! Bewaar de gecleande data
#-------------------------------
saveRDS(df_inputs,
        file.path(str_replace(cleaned_input_fasta, ".fasta", '.RDS' )))

create_input_fasta(file = cleaned_input_fasta, 
                   lowercase = TRUE, 
                   data = df_inputs)
