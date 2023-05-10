refdb_name <- 'refdb_2023-06-01'
ecopcr_riaz_file <- "amplified_clean.fasta"
merged_riaz_file <- "amplified_clean_uniq.fasta"
input_fasta_file <- "input.fasta"
path <- "output"


library(tidyverse)
source("scripts/_functions_fasta.R")

## ECOPCR INHOUD

input_data <- 
  parse_refdb_fasta(file.path("database", refdb_name, input_fasta_file), 
                    is_merged_file = FALSE)
ecopcr_data <- 
  parse_refdb_fasta(file.path("database", refdb_name, ecopcr_riaz_file), 
                    is_merged_file = FALSE) %>% 
  rename(taxid = TAXID,
         amplicon_hash = AMPLICON_HASH,
         coubt = COUNT)
merged_data <- 
  parse_refdb_fasta(file.path("database", refdb_name, merged_riaz_file), 
                    is_merged_file = TRUE) %>% 
  rename(merged_count = COUNT, 
         amplicon_hash = AMPLICON_HASH,
         merged_taxid = MERGED_TAXID,
         taxid = TAXID
         )

ecopcr_combined <- combine_ecpocr_with_merged(ecopcr_data, merged_data) 
table(ecopcr_combined$is_merged)
table(ecopcr_combined$obi_rank)
table(ecopcr_combined$obi_rank, ecopcr_combined$is_merged)

### overzicht van soorten en families die in de ecopcr zitten (origineel, nog niet gemerged)

ecopcr_data %>% 
  group_by(family_name) %>% 
  summarise(verschillende_genera = n_distinct(genus), 
            verschillende_taxa = n_distinct(taxid)) %>% 
  arrange(desc(verschillende_taxa)) %>% 
  write_excel_csv2(file = paste0(path, "taxa_in_ecopcr.csv"))


### Overzicht sequenties die in input maar niet in de output zitten
input_data %>% anti_join(ecopcr_data, by = "genlab_id") %>% view()


###############################################################################

df_soortenlijst <- read_sheet(metadata_gdrive_key, sheet = "Soortenlijst") %>% 
  filter(Priority %in% c(1,2), Taxid > 0) %>% 
  rename(taxid = Taxid, priority = Priority)

df_ok_merges <- read_sheet(metadata_gdrive_key, sheet = "Toegelaten_merges") 

#overzicht van records die niet op soort gebracht konden worden
ecopcr_filtered_sp <- ecopcr_combined %>% 
  inner_join(df_soortenlijst %>% select(taxid, priority)) %>% 
  select(genlab_id, taxid, rank, priority, species_name,
         amplicon_hash, obi_rank, obi_taxid, genus_name, family_name,
         merged_overview, obi_count) %>% 
  filter(!(obi_rank %in% c("subspecies", "species"))) %>% 
  group_by(across(-genlab_id)) %>% 
  summarise(genlab_id = paste(genlab_id, collapse = ";")) %>% 
  arrange(obi_taxid, taxid, amplicon_hash) %>% 
  dplyr::filter(!(obi_rank == "genus" & 
            obi_taxid %in% (df_ok_merges %>% 
                              filter(RANK == "genus") %>% pull(TAXID))))

# Kijk of er soorten mergebaar zijn (prioriteit 1)
ecopcr_filtered_sp %>% group_by(obi_taxid) %>% 
  do(check_multihits2(.)) %>% view()


check_multihits2 <- function(data) {
  np1 <- sum(data$priority == 1)
  np2 <- sum(data$priority == 2)
  np3 <- sum(data$priority == 3)
  nrows <- nrow(data)
  data$EVAL <- NA
  data$NEWID <- NA
  data$NEWRANK <- NA
  if (np1 + np2 + np3 == 1) {
    data$EVAL <- "OK, unique"
    data$NEWID <- data$taxid
    data$NEWRANK <- data$rank   
  }
  if (np1 == 1 & np2 + np3 > 0) {
    whi1 <- which(data$priority == 1)
    ntx <- data$taxid[whi1]
    data$EVAL[whi1] <- "OK, Pref species"
    data$EVAL[-whi1] <- "multihit_with_p1"
    data$NEWID <- ntx
    data$NEWRANK[whi1] <- data$rank[whi1]
  }
  if (np1 > 1) {
    data$EVAL <- "CLASH P1"
    data$NEWID <- NA
    data$NEWRANK <- NA
  }
  if (np1 == 0 & np2 + np3 == 1 ){
    data$EVAL <- "OK, in this db"
    data$NEWID <- data$taxid
    data$NEWRANK <- data$rank
  }
  if (np1 == 0 & np2 + np3 > 1) {
    #to be implemented
  }
  data
}
