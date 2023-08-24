

source("scripts/refdb_teleo_00_initialisation.R")

## SOORTENLIJST
## --------------

df_species_all <- read_sheet(metadata_gdrive_key, "Soortenlijst") %>% 
  mutate(ROWNR = 1:n())

df_species <- df_species_all %>% 
  filter(Taxid > 0, !is.na(Priority), Priority != 9)

df_priority <- df_species %>% 
  select(Taxid, Priority)

## INPUTS
## --------

files <- sort(list.files(fasta_inputs, pattern = ".fasta"))

df_inputs_orig <- NULL
for (file in files) {
  cat("\n\nINLEZEN VAN ", file, "\n--------------------------------\n")
  parsed <- parse_refdb_fasta(file.path(fasta_inputs, file))
  df_inputs_orig <- df_inputs_orig %>% 
    bind_rows(parsed)
}

## AFGELEIDE FILES
## ---------------

#nakijken of de sequenties wel gevonden worden (voorlopig 21/08/23 is dat niet het geval)
df_passlist <- read_sheet(metadata_gdrive_key, "Passlist_Teleo")
df_passlist <- df_inputs_orig %>% 
  filter(genbank_id %in% df_passlist$ENTRY_ID) %>% 
  mutate(taxid = as.numeric(taxid))

df_seq_errors <- read_sheet(metadata_gdrive_key, "Sequentiefouten")

df_multihit <-  read_sheet(metadata_gdrive_key, "Multihitlist_Teleo", range = "A:H", col_types = 'ccccccnn')

df_allowed_merges <- read_sheet(metadata_gdrive_key, "Toegelaten_merges_Teleo")


### INPUTS

#get the source files

#read the source files

#nakijken of er duplicaten zijn, en zeker tegenstrijdigheden opsporen
df_inputs_all <- df_inputs_orig %>% 
  mutate(genbank_id = stringr::str_trim(genbank_id),
         taxid = trimws(df_inputs_orig$taxid))

whidup <- which(duplicated(df_inputs_all %>%  select(genbank_id, taxid, dna_sequence)))
if (length(whidup)) df_inputs_all <- df_inputs_all[-whidup, ]


## CURATE INPUTS

#Kijk of er nog prioriteit-9 (soorten die niet weerhouden worden wegens bvb hybriden) aan seq_errors moet toegevoegd worden

#Haal foutieve sequenties uit 
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
write_excel_csv2(dup_ids, paste0(output_path, '/', "duplicate_sequenties_verschillend_taxid.csv"))

saveRDS(df_inputs, paste0(output_path, '/', "inputs_before_multihit.RDS"))

### remove alle sequenties die in multihit taxa staan en waar taxid != pref_taxid
#multihitlist van teleo gebruiken is OK, want multihits van riaz hoeven geen multihits te zijn bij teleo
df_inputs <- df_inputs%>% 
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

