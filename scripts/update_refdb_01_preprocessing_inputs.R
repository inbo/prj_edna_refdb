#FILE DOES NOT ALTER FILES OR SAVES ANYTHING
#IT JUST CHECKS INPUTS FOR PROBLEMS

source("scripts/update_refdb_00_init.R")

## SPECIES
df_species_all <- read_sheet(input_source, "Soortenlijst") %>% 
  mutate(ROWNR = 1:n())

df_species <- df_species_all %>% 
  filter(Taxid > 0, !is.na(Priority))

df_priority <- df_species %>% 
  select(Taxid, Priority)

### MULTIHIT

df_multihit <-  read_sheet(input_source, "Multihitlist_riaz")

### SEQ_FOUTEN

df_seq_errors <- read_sheet(input_source, "Sequentie_fouten")

### INPUT SEQS


files <- list.files(fasta_input_path, pattern = ".fasta")
if (!length(files)) {
  warning("using local input files")
  files <- list.files(fasta_input_path_local, pattern = ".fasta")
  rootpath <- file.path(fasta_input_path_local, file)
} else {
  rootpath <- file.path(fasta_input_path, file) 
}

df_inputs_all <- NULL
for (file in files) {
  cat("\n\nINLEZEN VAN ", file, "\n--------------------------------\n")
  parsed <- parse_refdb_fasta(rootpath)
  df_inputs_all <- df_inputs_all %>% 
    bind_rows(parsed)
  rm(fullpath)
}
#remove full duplicates (same name, amplicon and taxid)
df_inputs <- df_inputs_all %>% 
  group_by(genlab_id, taxid, dna_sequence) %>% 
  summarise(amplicon_hash = AMPLICON_HASH[1])

#duplicate input sequence names with different dna_sequence or taxid
dup_ids <- df_inputs %>% group_by(genlab_id) %>% 
  summarise(aantal_seqs = n(), 
            sources = paste(source, collapse = ","),
            aantal_taxa = n_distinct(taxid), 
            taxa = paste(unique(taxid), collapse = ",")) %>% 
  filter(aantal_seqs > 1 | aantal_taxa > 1)
  

dup_ids <- df_inputs %>% 
  filter(duplicated(genlab_id) | str_detect(genlab_id, ".PV")) %>% pull(genlab_id)
df_inputs %>% filter(genlab_id %in% dup_ids) %>% 
  select(source, genlab_id, taxid, AMPLICON_HASH, dna_sequence) %>% 
  #filter(!duplicated(.)) %>% 
  arrange(genlab_id) %>% view()

#duplicate input DNA seqs
dup_hash <- df_inputs %>% filter(duplicated(AMPLICON_HASH)) %>% pull(AMPLICON_HASH)
dupli_seqs <- df_inputs %>% 
  filter(AMPLICON_HASH %in% dup_hash) %>% 
  select(genlab_id, taxid, AMPLICON_HASH, source, dna_sequence) %>% 
  arrange(AMPLICON_HASH, genlab_id)

#all duplicate sequences
dupli_seqsg <- dupli_seqs %>% 
  group_by(AMPLICON_HASH) %>% 
  summarize(aantal_inputs = n(),
            distinct_inputs = n_distinct(genlab_id),
            aantal_taxa = n_distinct(taxid), 
            taxa = paste(unique(taxid), collapse = ","),
            inputs = paste(unique(genlab_id), collapse = ",")) %>% 
  filter(distinct_inputs > 1 | aantal_taxa > 1) %>% 
  arrange(aantal_taxa, distinct_inputs)

#a priori detected sequences
dupli_seqs_filt <- dupli_seqs %>% 
  filter(taxid %in% df_multihit$TAXID) %>% 
  filter(!(genlab_id %in% df_seq_errors$ENTRY_ID)) %>% 
  filter(!(taxid %in% (df_seq_errors %>% 
                        filter(TYPE == "taxid") %>% 
                        pull(TAXID)))) %>% 
  group_by(AMPLICON_HASH) %>% 
  summarize(aantal_inputs = n(),
            distinct_inputs = n_distinct(genlab_id),
            aantal_taxa = n_distinct(taxid), 
            taxa = paste(unique(taxid), collapse = ","),
            inputs = paste(unique(genlab_id), collapse = ",")) %>% 
  filter(aantal_taxa > 1) %>% 
  arrange(aantal_taxa, distinct_inputs)



  