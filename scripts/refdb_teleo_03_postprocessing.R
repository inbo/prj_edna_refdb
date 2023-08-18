####TELEO

#db_name uit initialisation
before_ecopcr_file <- "amplified.fasta"
ecopcr_file <- "amplified_clean.fasta"
merged_file <- "final_db_0.99.fasta"
input_fasta_file <- "input.fasta"


library(tidyverse)
library(googlesheets4)
source("scripts/_functions_fasta.R")


df_soortenlijst_all <- read_sheet(metadata_gdrive_key, sheet = "Soortenlijst") %>% 
  filter(Taxid > 0) %>% 
  rename(taxid = Taxid, priority = Priority)

df_soortenlijst <- df_soortenlijst_all %>% 
  filter(priority %in% c(1:7))

df_ok_merges <- read_sheet(metadata_gdrive_key, sheet = "Toegelaten_merges_Teleo") 

## ECOPCR INHOUD

input_data <- 
  parse_refdb_fasta(file.path("database", db_name, "kept_input.fasta"), 
                    is_merged_file = FALSE)

input_data_before_multihit <- readRDS(file.path(output_path, 'inputs_before_multihit.RDS'))

ecopcr_data <- 
  parse_refdb_fasta(file.path("database", db_name, ecopcr_file), 
                    is_merged_file = FALSE) %>% 
  rename(taxid = TAXID,
         amplicon_hash = AMPLICON_HASH,
         count = COUNT)
merged_data <- 
  parse_refdb_fasta(file.path("database", db_name, merged_file), 
                    is_merged_file = TRUE) %>% 
  rename(amplicon_hash = AMPLICON_HASH,
         merged_taxid = MERGED_TAXID,
         taxid = TAXID
         )
ecopcr_combined <- combine_ecpocr_with_merged(ecopcr_data, merged_data) 

write_rds(ecopcr_data, "ecopcr_data.RDS")
write_rds(merged_data, "merged_data.RDS")
write_rds(input_data, "input_data.RDS")
write_rds(ecopcr_combined, "ecopcr_combined.RDS")

table(ecopcr_combined$is_merged)
table(ecopcr_combined$obi_rank)
table(ecopcr_combined$obi_rank, ecopcr_combined$is_merged)


### Tijdelijke details amplicons
# 
# probleemsoorten1 <- ecopcr_combined %>% 
#   filter(taxid %in% c(8047,8049,8052,8056,8058,8060,44932,80720,80721,185735,185739,1042646)) %>% 
#   arrange(amplicon_hash) %>% 
#   select(genlab_id, amplicon, amplicon_hash, rank, taxid, species_name, merged_overview, obi_rank, obi_taxid)
# write_excel_csv2(probleemsoorten1, "multihit_1.csv")
# 
# probleemsoorten2 <- ecopcr_combined %>%
#   filter(taxid %in% c(8056,44932,185735)) %>% 
#   arrange(amplicon_hash) %>% 
#   select(genlab_id, amplicon, amplicon_hash, rank, taxid, species_name, merged_overview, obi_rank, obi_taxid)
# write_excel_csv2(probleemsoorten2, "multihit_2.csv")

### overzicht van soorten en families die in de ecopcr zitten (origineel, nog niet gemerged)

ecopcr_data %>% 
  group_by(family_name) %>% 
  summarise(verschillende_genera = n_distinct(genus), 
            verschillende_taxa = n_distinct(taxid)) %>% 
  arrange(desc(verschillende_taxa)) %>% 
  write_excel_csv2(file = paste0(output_path,"/", "taxa_in_ecopcr.csv"))


### Overzicht sequenties die in input maar niet in de output zitten
niet_geamplificeerd <- (input_data %>% anti_join(ecopcr_data, by = "genlab_id") %>% 
  mutate(taxid = as.numeric(TAXID)) %>% 
  left_join(df_soortenlijst %>% 
              select(taxid, priority, NameScientific, NameEnglish))) %>%
  write_excel_csv2(file = paste0(output_path,"/", "sequenties_niet_geamplificeerd.csv"))


mtch <- NULL
for (i in seq_len(nrow(niet_geamplificeerd)))
  mtch[i] <- find_matching_seq(niet_geamplificeerd$dna_sequence[i])
table(mtch)
geweigerd_ondanks_match <- niet_geamplificeerd %>% 
  mutate(gematched = mtch) %>% 
  filter(substring(mtch, 1,3) == "1_1")
view(geweigerd_ondanks_match)
write_excel_csv2(geweigerd_ondanks_match, 
                 file = file.path(output_path, "geweigerd_ondanks_perfecte_match.csv"))

#tijdelijke check (mag weg)
before_data <- 
  parse_refdb_fasta(file.path("database", db_name, before_ecopcr_file), 
                    is_merged_file = FALSE) %>% 
  rename(taxid = TAXID,
         amplicon_hash = AMPLICON_HASH,
         count = COUNT)
#species, genus en family  moet gekend zijn voor de ecopcr, soms is species wel gekend, maar genus of familie niet, en die gevallen worden verwijderd
lost_in_cleaning <- before_data %>% filter(!(genlab_id %in% ecopcr_data$genlab_id))
write_excel_csv2(lost_in_cleaning, 
                 file = file.path(output_path, "geweigerd_onvolledige_taxdump_hierarchie.csv"))






###############################################################################
### CLEANING
###############################################################################

#overzicht van records die niet op soort gebracht konden worden
#indien 0 rijen dan is alles op soort gebracht of een toegelaten genusmerge of toegelaten familymerge
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
                  obi_taxid %in% (df_ok_merges %>% filter(RANK == "genus") %>% pull(TAXID))), 
                !(obi_rank == "family" & 
                  obi_taxid %in% (df_ok_merges %>% filter(RANK == "family") %>% pull(TAXID))  ))

# Kijk of er soorten mergebaar zijn (prioriteit 1)
niet_op_soort <- (ecopcr_filtered_sp %>% group_by(obi_taxid) %>% 
  do(check_multihits2(.))) %>% 
  write_excel_csv2(file = paste0(output_path,"/", "niet_op_soort_gebracht.csv"))

## ANDERE MOGELIJKHEID
#overzicht soorten
df_spm <- ecopcr_filtered_sp %>%
  arrange(obi_taxid) %>%
  group_by(obi_taxid) %>%
 do ({
   get_taxa_from_merged(., df_soortenlijst)
 }
 )

df_beoordeeld <- df_spm %>% group_by(obi_taxid) %>%
  do ({
    judge_species(.)
  })
write_excel_csv2(df_beoordeeld, file = paste0(output_path,"/", "niet_op_soort_gebracht_alternatief.csv"))

###################################################################
### SOORTENTABEL
###################################################################



df_soortentabel <- df_soortenlijst %>% 
  group_by(Group, NameScientific, NameEnglish, NameDutch, taxid, priority) %>% 
  summarise(n_input_orig = sum(input_data_before_multihit$taxid == .data$taxid[1]),
            #n_input_seqs = sum(as.numeric(input_data$TAXID) == .data$taxid[1]),
            n_input_seqs = sum(as.numeric(input_data$TAXID == .data$taxid[1]) > 0, na.rm = TRUE),
            n_amplicons = sum(ecopcr_combined$taxid == .data$taxid[1]),
            n_determined = sum(ecopcr_combined$obi_taxid == .data$taxid[1]))

ecopcr_soorten <- ecopcr_combined %>% 
  right_join(df_soortenlijst %>% select(taxid)) %>% 
  group_by(taxid) %>% 
  summarise(merged = sum(is_merged)>0, 
         obi_taxid = paste(unique(obi_taxid), collapse = '|'),
         obi_rank = paste(unique(obi_rank), collapse = "|")) %>% 
  mutate(merged = merged & taxid != obi_taxid)

beoordeling <- df_beoordeeld %>% group_by(taxid, pref_taxid) %>%
  summarise(oordeel = paste(unique(oordeel), collapse = "|"))


df_multihits <-
  df_beoordeeld %>% group_by(taxid) %>%
  do({
    taxid <- .$taxid[1]
    obitaxids <- (unique(.$obi_taxid))
    taxids <- df_beoordeeld$taxid[df_beoordeeld$obi_taxid %in% obitaxids]
    data.frame(multihit_taxa = paste("|", sort(unique(taxids)), collapse = '', sep = ''))
  })

df_soortenevaluatie <- df_soortentabel %>% 
  left_join(ecopcr_soorten) %>% 
  left_join(beoordeling) %>% 
  left_join(df_multihits) %>% 
  mutate(pref_taxid = ifelse(merged == FALSE, taxid, pref_taxid)) %>% #klopt niet, merged is fout
  left_join(df_soortenlijst %>% select(taxid, NameEnglish, NameScientific), by = c("pref_taxid" = "taxid"))

write_excel_csv2(df_soortenevaluatie, file = paste0(output_path,"/", "soortenevaluatie_ruw.csv"))





















