
#SESSION

source("scripts/update_refdb_00_init.R")

######################################################################

# INPUT GSHEET

df_species <- read_sheet(input_source, "Soortenlijst") %>% 
  mutate(ROWNR = 1:n()) %>% 
  filter(Taxid > 0, !is.na(Priority))

df_priority <- df_species %>% 
  select(Taxid, Priority)

df_multihit <-  read_sheet(input_source, "Multihitlist_riaz")
df_seq_errors <- read_sheet(input_source, "Sequentie_fouten")
df_allowed_merges <- read_sheet(input_source, "Toegelaten_merges")
df_passlist <- read_sheet(input_source, "Pass_list_Riaz")

## INPUTSEQS 

### LEES ALLE INPUTS,
#VOEG MULTIHIT_PREF_TAXID SEQS TOE en VOEG PASS_LIST SEQS TOE
files <- list.files(fasta_input_path, pattern = ".fasta")
if (!length(files)) {
  warning("inputbestanden niet gevonden: gedeelde map van gendiv in GDrive. Lokale kopie gebruikt")
  fasta_input_path <- 'input'
  files <- list.files(fasta_input_path, pattern = ".fasta")
  print(files)
}

df_inputs_all <- NULL
for (file in files) {
  fullpath <- file.path(fasta_input_path, file)
  cat("\n\nINLEZEN VAN ", file, "\n--------------------------------\n")
  parsed <- parse_refdb_fasta(fullpath)
  df_inputs_all <- df_inputs_all %>% 
    bind_rows(parsed)
  rm(fullpath)
}

## CURATE INPUTS

passed <- df_inputs_all %>% 
  filter(genlab_id %in% (df_passlist %>% pull(ENTRY_ID)))

### remove seq errors
df_inputs <- df_inputs_all %>% 
  filter(!(genlab_id %in% (df_seq_errors %>% pull(ENTRY_ID)))) %>% 
  filter(!taxid %in% (df_seq_errors %>% 
                        filter(TYPE == "taxid") %>% pull(TAXID)))

### remove alle sequenties die in multihit taxa staan en waar taxid != pref_taxid
df_inputs <- df_inputs %>% 
  filter(!taxid %in% (df_multihit %>% 
                        filter(TAXID != PREF_TAXID) %>% 
                        pull(TAXID)))

### add passlist sequences
df_inputs <- df_inputs %>% bind_rows(passed) 

### remove full duplicates
df_inputs <- df_inputs %>% distinct(genlab_id, taxid, AMPLICON_HASH, 
                                           .keep_all = TRUE)

## create fasta input for ecopcr

create_input_fasta(file = file.path("database", db_name, fasta_name), 
                   data = df_inputs)

#voor obitools3
create_input_fasta(file = file.path("database", db_name, fasta_name_obi3), 
                   lowercase = TRUE, 
                   data = df_inputs)

