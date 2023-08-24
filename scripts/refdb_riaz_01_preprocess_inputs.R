
## GOOGLESHEET REFDB

df_species_all <- read_sheet(metadata_gdrive_key, "Soortenlijst") %>% 
  mutate(ROWNR = 1:n())

df_species <- df_species_all %>% 
  filter(Taxid > 0, !is.na(Priority), Priority != 9)

df_priority <- df_species %>% 
  select(Taxid, Priority)

df_passlist <- read_sheet(metadata_gdrive_key, "Passlist_Riaz")

df_seq_errors <- read_sheet(metadata_gdrive_key, "Sequentiefouten")

df_multihit <-  read_sheet(metadata_gdrive_key, "Multihitlist_Riaz")

df_allowed_merges <- read_sheet(metadata_gdrive_key, "Toegelaten_merges_Riaz")


### INPUTS

#get the source files
files <- sort(list.files(fasta_inputs, pattern = ".fasta"))

#read the source files
df_inputs_all <- df_inputs_orig <- NULL
for (file in files) {
  cat("\n\nINLEZEN VAN ", file, "\n--------------------------------\n")
  parsed <- parse_refdb_fasta(file.path(fasta_inputs, file))
  df_inputs_all <- df_inputs_all %>% 
    bind_rows(parsed)
}
df_inputs_orig <- df_inputs_all

#nakijken of er duplicaten zijn, en zeker tegenstrijdigheden opsporen
df_inputs_all <- mutate(df_inputs_all, 
                        genbank_id = substring(stringr::str_trim(genbank_id),1,19),
                        taxid = as.numeric(trimws(df_inputs_all$taxid)))



whidup <- which(duplicated(df_inputs_all %>%  select(genbank_id, taxid, dna_sequence)))
df_inputs_all <- df_inputs_all[-whidup, ]

dup_id_diff_dna <- df_inputs_all$genbank_id[duplicated(df_inputs_all$genbank_id)]
dup_id_diff_dna <- df_inputs_all %>% filter(genbank_id %in% dup_id_diff_dna) %>% arrange(genbank_id, source)
write_excel_csv2(dup_id_diff_dna, file = file.path(output_path, 'duplicate_genlabid_verschillend_dna.csv'))

#Kijk of er nog prioriteit-9 (soorten die niet weerhouden worden wegens bvb hybriden) aan seq_errors moet toegevoegd worden
df_inputs_all %>% filter(taxid %in% (df_species %>% filter(Priority == 9) %>% pull(Taxid)))

#let op: De sequentiefouten van Riaz worden hier al verwijderd in tegenstelling tot 2019
df_inputs_all <- df_inputs_all %>% 
  mutate(taxid = as.numeric(taxid)) %>% 
  filter(!(genbank_id %in% (df_seq_errors %>% pull(ENTRY_ID)))) %>% 
  filter(!taxid %in% (df_seq_errors %>% 
                        filter(TYPE == "taxid") %>% pull(TAXID)))

#find duplicate ids with different species
dup_ids <- df_inputs_all %>% group_by(genbank_id) %>% 
  summarise(aantal_seqs = n(), 
            sources = paste(source, collapse = ","),
            aantal_taxa = n_distinct(taxid), 
            taxa = paste(unique(taxid), collapse = ",")) %>% 
  filter(aantal_seqs > 1 & aantal_taxa > 1)
write_excel_csv2(dup_ids, paste0(output_path, '/', "duplicate_sequenties_verschillend_taxid.csv"))

## CURATE INPUTS

saveRDS(df_inputs_all, file.path("database", db_name, "inputs_before_multihit.RDS"))

### remove alle sequenties die in multihit taxa staan en waar taxid != pref_taxid
#multihitlist van teleo gebruiken is OK, want multihits van riaz hoeven geen multihits te zijn bij teleo
df_inputs <- df_inputs_all %>% 
  filter(!taxid %in% (df_multihit %>% 
                        filter(TAXID != PREF_TAXID) %>% 
                        pull(TAXID)))

### add passlist sequences
df_inputs <- df_inputs %>% bind_rows(df_passlist) 

### remove full duplicates
df_inputs <- df_inputs %>% distinct(genbank_id, taxid, DNA_HASH, 
                                    .keep_all = TRUE)

## create fasta input for ecopcr (obitools2)

#create_input_fasta(file = file.path("database", db_name, fasta_name), data = df_inputs)

#voor obitools3
create_input_fasta(file = file.path("database", db_name, fasta_name), 
                   lowercase = TRUE, 
                   data = df_inputs)

saveRDS(df_inputs, file=file.path(output_path, "df_inputs.RDS"))

