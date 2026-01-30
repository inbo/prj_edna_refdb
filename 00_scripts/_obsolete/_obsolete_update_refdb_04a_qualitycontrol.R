
refdb_name <- 'refdb_2022-11-11'
ecopcr_riaz_file <- "amplified_clean.fasta"
merged_riaz_file <- "amplified_clean_uniq.fasta"
input_fasta_file <- "input.fasta"
path <- "output"


library(tidyverse)
source("scripts/_functions_fasta.R")


#(INSPIRATIE UIT 02_Ecopcr_refdb_20220517.R)

## ECOPCR INHOUD

input_data <- 
  parse_refdb_fasta(file.path("database", refdb_name, input_fasta_file), 
                    is_merged_file = FALSE)
ecopcr_data <- 
  parse_refdb_fasta(file.path("database", refdb_name, ecopcr_riaz_file), 
                    is_merged_file = FALSE)
merged_data <- 
  parse_refdb_fasta(file.path("database", refdb_name, merged_riaz_file), 
                    is_merged_file = TRUE) 

ecopcr_riaz_db <- combine_ecpocr_with_merged(ecopcr_data, merged_data) 
table(ecopcr_riaz_tmp$IS_MERGED)
table(ecopcr_riaz_tmp$OBI_RANK)

### overzicht van soorten en families die in de ecopcr zitten (origineel, nog niet gemerged)

ecopcr_data %>% 
  group_by(family_name) %>% 
  summarise(verschillende_genera = n_distinct(genus), 
            verschillende_taxa = n_distinct(taxid)) %>% 
  arrange(desc(verschillende_taxa)) %>% 
  write_excel_csv2(file = paste0(path, "taxa_in_ecopcr.csv"))

### soorten uit de soortenlijst die niet in de ecopcr zitten



### inputs die niet in de ecopcr zitten (dit is enkel bij het ecopcr proces)

#ofwel niet geamplificeerd (instelling error rate e, l, L, ...)
#ofwel geen consistente naamgeving in NCBI taxdump

input_data %>% anti_join(ecopcr_data, by = "genlab_id") %>% view()

### inputs die in ecopcr zitten, maar niet in merged_data 


## NIEUWE MULTIHITS

### Vindt nieuwe multihit soorten

find_multihit_taxa

### Hou rekening met allowed_merges

# SAMENVATTENDE TABEL

## tabel

## nieuwe bloksoorten, ...

## lijst van soorten in refdb (na curatie)

#VOEG TOE AAN SQLITE DB (Dit is voor het volgende script)
