####TELEO

#db_name uit initialisation
before_ecopcr_file <- "amplified.fasta"
ecopcr_file <- "amplified_clean.fasta"
merged_file <- "final_db_0.99.fasta"
input_fasta_file <- "input.fasta"


library(tidyverse)
library(googlesheets4)
source("scripts/_functions_fasta.R")
source("scripts/refdb_teleo_00_initialisation.R")


df_soortenlijst_all <- read_sheet(metadata_gdrive_key, sheet = "Soortenlijst") %>% 
  filter(Taxid > 0) %>% 
  rename(taxid = Taxid, priority = Priority)

df_soortenlijst <- df_soortenlijst_all %>% 
  filter(priority %in% c(1:7))

df_ok_merges <- read_sheet(metadata_gdrive_key, sheet = "Toegelaten_merges_Teleo") 

df_multihits <- read_sheet(metadata_gdrive_key, sheet = "Multihitlist_Teleo") 

## ECOPCR INHOUD

input_data <- 
  parse_refdb_fasta(file.path("database", db_name, "kept_input.fasta"), 
                    is_merged_file = FALSE)

input_data_before_multihit <- readRDS(file.path("database", db_name, 'inputs_before_multihit.RDS'))

ecopcr_data <- 
  parse_refdb_fasta(file.path("database", db_name, ecopcr_file), 
                    is_merged_file = FALSE) %>% 
  rename(taxid = TAXID,
         dna_hash = DNA_HASH,
         count = COUNT)
merged_data <- 
  parse_refdb_fasta(file.path("database", db_name, merged_file), 
                    is_merged_file = TRUE) %>% 
  rename(dna_hash = DNA_HASH,
         merged_taxid = MERGED_TAXID,
         taxid = TAXID
         )
ecopcr_combined <- combine_ecpocr_with_merged(ecopcr_data, merged_data) 

write_rds(df_soortenlijst, file.path("database", db_name, "df_soortenlijst.RDS"))
write_rds(df_ok_merges, file.path("database", db_name, "df_ok_merges.RDS"))
write_rds(ecopcr_data, file.path("database", db_name, "ecopcr_data.RDS"))
write_rds(merged_data, file.path("database", db_name, "merged_data.RDS"))
write_rds(input_data, file.path("database", db_name, "input_data.RDS"))
write_rds(ecopcr_combined, file.path("database", db_name, "ecopcr_combined.RDS"))
write_rds(df_multihits, file.path("database", db_name, "df_multihtitlist.RDS"))

table(ecopcr_combined$is_merged)
table(ecopcr_combined$obi_rank)
table(ecopcr_combined$obi_rank, ecopcr_combined$is_merged)


### Tijdelijke details amplicons
# 
# probleemsoorten1 <- ecopcr_combined %>% 
#   filter(taxid %in% c(8047,8049,8052,8056,8058,8060,44932,80720,80721,185735,185739,1042646)) %>% 
#   arrange(amplicon_hash) %>% 
#   select(genbank_id, amplicon, amplicon_hash, rank, taxid, species_name, merged_overview, obi_rank, obi_taxid)
# write_excel_csv2(probleemsoorten1, "multihit_1.csv")
# 
# probleemsoorten2 <- ecopcr_combined %>%
#   filter(taxid %in% c(8056,44932,185735)) %>% 
#   arrange(amplicon_hash) %>% 
#   select(genbank_id, amplicon, amplicon_hash, rank, taxid, species_name, merged_overview, obi_rank, obi_taxid)
# write_excel_csv2(probleemsoorten2, "multihit_2.csv")

### overzicht van soorten en families die in de ecopcr zitten (origineel, nog niet gemerged)

ecopcr_data %>% 
  group_by(family_name) %>% 
  summarise(verschillende_genera = n_distinct(genus), 
            verschillende_taxa = n_distinct(taxid)) %>% 
  arrange(desc(verschillende_taxa)) %>% 
  write_excel_csv2(file = paste0(output_path,"/", "taxa_in_ecopcr.csv"))




#### Overzicht sequenties die in input maar niet in de output zitten
niet_geamplificeerd <- (input_data %>% anti_join(ecopcr_data, by = "genbank_id") %>% 
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
# before_data <- 
#   parse_refdb_fasta(file.path("database", db_name, before_ecopcr_file), 
#                     is_merged_file = FALSE) %>% 
#   rename(taxid = TAXID,
#          dna_hash = DNA_HASH,
#          count = COUNT)
# #species, genus en family  moet gekend zijn voor de ecopcr, soms is species wel gekend, maar genus of familie niet, en die gevallen worden verwijderd
# lost_in_cleaning <- before_data %>% filter(!(genbank_id %in% ecopcr_data$genbank_id))
# write_excel_csv2(lost_in_cleaning, 
#                  file = file.path(output_path, "geweigerd_onvolledige_taxdump_hierarchie.csv"))



###############################################################################
### CONFLICTEN OPSPOREN
###############################################################################

from_rds <- FALSE
if (from_rds) {
  library(tidyverse)
  source("scripts/_functions_fasta.R")
  source("scripts/refdb_teleo_00_initialisation.R")
  df_soortenlijst <- read_rds(file.path("database", db_name, "df_soortenlijst.RDS"))
  df_ok_merges <- read_rds(file.path("database", db_name, "df_ok_merges.RDS"))
  #ecopcr_data <- read_rds(file.path("database", db_name, "ecopcr_data.RDS"))
  #merged_data <- read_rds(file.path("database", db_name, "merged_data.RDS"))
  #input_data <- read_rds(file.path("database", db_name, "input_data.RDS"))
  ecopcr_combined <- read_rds(file.path("database", db_name, "ecopcr_combined.RDS")) 
  input_data_before_multihit <- read_rds(file.path("database", db_name,  'inputs_before_multihit.RDS'))
  df_multihits <- read_rds(file.path("database", db_name, "df_multihtitlist.RDS"))
  
}


#overzicht van records die niet op soort gebracht konden worden
#indien 0 rijen dan is alles op soort gebracht of een toegelaten genusmerge of toegelaten familymerge
ecopcr_filtered_sp <- ecopcr_combined %>% 
  inner_join(df_soortenlijst %>% select(taxid, priority)) %>% 
  select(genbank_id, taxid, rank, priority, species_name,
         dna_hash, obi_rank, obi_taxid, genus_name, family_name,
         merged_overview, obi_count) %>% 
  filter(!(obi_rank %in% c("subspecies", "species"))) %>% 
  group_by(across(-genbank_id)) %>% 
  summarise(genbank_id = paste(genbank_id, collapse = ";"), .groups = "drop") %>% 
  arrange(obi_taxid, taxid, dna_hash) %>% 
  dplyr::filter(!(obi_rank == "genus" & 
                  obi_taxid %in% (df_ok_merges %>% filter(RANK == "genus") %>% pull(TAXID))), 
                !(obi_rank == "family" & 
                  obi_taxid %in% (df_ok_merges %>% filter(RANK == "family") %>% pull(TAXID))  ))


#overzicht soorten
df_beoordeeld <- 
  ecopcr_filtered_sp %>%
  arrange(obi_taxid) %>%
  group_by(obi_taxid) %>%
 do ({
   get_taxa_from_merged(., df_soortenlijst)
 }) %>% 
  group_by(obi_taxid) %>% 
  do({
    judge_species(.)    
  })
  
write_excel_csv2(df_beoordeeld, file = paste0(output_path,"/", "niet_op_soort_gebracht_alternatief.csv"))

###################################################################
### SOORTENEVALUATIE
###################################################################

#soortentabel linken aan ecopcr output
ecopcr_soorten <- ecopcr_combined %>% 
  right_join(df_soortenlijst %>% select(taxid)) %>% 
  group_by(taxid) %>% 
  summarise(merged = sum(is_merged)>0, 
         obi_taxid = paste(unique(obi_taxid), collapse = '|'),
         obi_rank = paste(unique(obi_rank), collapse = "|")) %>%
  mutate(obi_taxid = ifelse(obi_taxid == "NA", NA, obi_taxid),
         obi_rank = ifelse(obi_rank == "NA", NA, obi_rank)) %>% 
  mutate(merged = taxid != obi_taxid)


#multihit-conflicten opsporen
df_conflicts <-
  df_beoordeeld %>% group_by(taxid) %>%
  do({
    taxid <- .$taxid[1]
    obitaxids <- (unique(.$obi_taxid))
    taxids <- df_beoordeeld$taxid[df_beoordeeld$obi_taxid %in% obitaxids]
    data.frame(has_conflicts = TRUE, 
               conflict_taxa = paste(sort(unique(taxids)), collapse = '|'),
               conflict_oordeel = paste(unique(.$oordeel), collapse = "|"))
  })

#soorten-multihit
multihit_species <- df_multihits %>% 
  select(taxid = TAXID, pref_taxid = PREF_TAXID) %>% 
  inner_join(df_soortenlijst %>% select(taxid)) %>% 
  mutate(is_multihit = TRUE,
         is_chosen = taxid == pref_taxid) 

multihit_species <- multihit_species %>% 
  group_by(taxid) %>% 
  do({
    pref = .$pref_taxid
    taxids <- sort(multihit_species$taxid[pref == multihit_species$pref_taxid])
    whi <- which(taxids == pref )
    if(length(whi)) {
      taxids = c(taxids[whi], taxids[-whi])
    }
    cbind(., hitlist = paste(taxids, collapse = " | "))
  })
  
#samenvattende tabel

df_soortenevaluatie <- df_soortenlijst %>% 
  group_by(Group, NameScientific, NameEnglish, NameDutch, taxid, priority) %>% 
  summarise(n_input_orig = sum(input_data_before_multihit$taxid == .data$taxid[1]),
            n_haplotypes = sum(ecopcr_combined$taxid == .data$taxid[1])) %>% 
  left_join(ecopcr_soorten) %>% 
  left_join(multihit_species) %>% 
  mutate(pref_taxid = ifelse(is.na(pref_taxid) & (obi_taxid == taxid | obi_rank == "species"), 
                             taxid, 
                             pref_taxid)) %>% 
  left_join(df_soortenlijst %>% select(taxid, 
                                       Pref_NameEnglish = NameEnglish, 
                                       Pref_NameScientific = NameScientific), 
            by = c("pref_taxid" = "taxid")) %>% 
  left_join(df_conflicts)

write_excel_csv2(df_soortenevaluatie, file = paste0(output_path,"/", "soortenevaluatie_ruw.csv"))





















