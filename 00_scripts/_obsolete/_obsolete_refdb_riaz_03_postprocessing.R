####RIAZ

#db_name uit initialisation
before_ecopcr_file <- "amplified.fasta"
ecopcr_file <- "amplified_clean.fasta"
merged_file <- "final_db_0.99.fasta"
input_fasta_file <- "input.fasta"


library(tidyverse)
library(googlesheets4)
source("scripts/_functions_fasta.R")
source("scripts/refdb_riaz_00_initialisation.R")


df_soortenlijst_all <- read_sheet(metadata_gdrive_key, sheet = "Soortenlijst") %>% 
  filter(Taxid > 0) %>% 
  rename(taxid = Taxid, priority = Priority)

df_soortenlijst <- df_soortenlijst_all %>% 
  filter(priority %in% c(1:7))

df_ok_merges <- read_sheet(metadata_gdrive_key, sheet = "Toegelaten_merges_Riaz") %>% 
  rename(taxid = TAXID, rank = RANK)

df_multihits <- read_sheet(metadata_gdrive_key, sheet = "Multihitlist_Riaz") %>% 
  rename(taxid = TAXID, pref_taxid = PREF_TAXID)

## ECOPCR INHOUD

input_data <- #vermoedelijk ook gewoon de RDS halen uit de db
  parse_refdb_fasta(file.path("database", db_name, "kept_input.fasta"), 
                    is_merged_file = FALSE)

input_data_before_multihit <- readRDS(file.path("database", db_name, 'inputs_before_multihit.RDS'))

ecopcr_data <- 
  parse_refdb_fasta(file.path("database", db_name, ecopcr_file), 
                    is_merged_file = FALSE) %>% 
  rename(taxid = TAXID,
         dna_hash = DNA_HASH,
         count = COUNT) %>% 
  mutate(taxid = str_replace(taxid, ";", ""))
merged_data <- 
  parse_refdb_fasta(file.path("database", db_name, merged_file), 
                    is_merged_file = TRUE) %>% 
  rename(dna_hash = DNA_HASH,
         merged_taxid = MERGED_TAXID,
         taxid = TAXID) %>% 
  mutate(obi_taxid = str_replace(substring(LCA_TAXID, 2), "]", ""),
         taxid = str_replace(taxid, ";", ""))

ecopcr_combined <- combine_ecpocr_with_merged(ecopcr_data, merged_data) 

write_rds(df_soortenlijst, file.path("database", db_name, "df_soortenlijst.RDS"))
write_rds(df_ok_merges, file.path("database", db_name, "df_ok_merges.RDS"))
write_rds(df_multihits, file.path("database", db_name, "df_multihtitlist.RDS"))
write_rds(ecopcr_data, file.path("database", db_name, "ecopcr_data.RDS"))
write_rds(merged_data, file.path("database", db_name, "merged_data.RDS"))
write_rds(input_data, file.path("database", db_name, "input_data.RDS"))
write_rds(ecopcr_combined, file.path("database", db_name, "ecopcr_combined.RDS"))

table(ecopcr_combined$is_merged)
table(ecopcr_combined$obi_rank, ecopcr_combined$is_merged)


### overzicht van soorten en families die in de ecopcr zitten (origineel, nog niet gemerged)

ecopcr_data %>% 
  group_by(family_name) %>% 
  summarise(verschillende_genera = n_distinct(genus), 
            verschillende_taxa = n_distinct(taxid)) %>% 
  arrange(desc(verschillende_taxa)) %>% 
  write_excel_csv2(file = paste0(output_path,"/", "taxa_in_ecopcr.csv"))

#### Overzicht sequenties die in input maar niet in de output zitten
niet_geamplificeerd <- input_data %>% anti_join(ecopcr_data, by = "genbank_id") %>% 
  mutate(taxid = as.numeric(TAXID)) %>% 
  left_join(df_soortenlijst %>% 
              select(taxid, priority, NameScientific, NameEnglish), 
            by = "taxid")
mtch <- NULL
for (i in seq_len(nrow(niet_geamplificeerd)))
  mtch[i] <- find_matching_seq(niet_geamplificeerd$dna_sequence[i])[4]
table(mtch)
geweigerd_ondanks_match <- niet_geamplificeerd %>% 
  mutate(gematched = mtch) %>% 
  filter(mtch == 11)
view(geweigerd_ondanks_match)
write_excel_csv2(geweigerd_ondanks_match, 
                 file = file.path(output_path, "geweigerd_ondanks_perfecte_match.csv"))

###############################################################################
### CONFLICTEN OPSPOREN
###############################################################################

#importeren data

from_rds <- FALSE
if (from_rds) {
  library(tidyverse)
  source("scripts/_functions_fasta.R")
  source("scripts/refdb_riaz_00_initialisation.R")
  source("scripts/_functions_postprocessing.R")
  df_soortenlijst <- read_rds(file.path("database", db_name, "df_soortenlijst.RDS"))
  df_ok_merges <- read_rds(file.path("database", db_name, "df_ok_merges.RDS"))
  ecopcr_data <- read_rds(file.path("database", db_name, "ecopcr_data.RDS"))
  merged_data <- read_rds(file.path("database", db_name, "merged_data.RDS"))
  #input_data <- read_rds(file.path("database", db_name, "input_data.RDS"))
  ecopcr_combined <- read_rds(file.path("database", db_name, "ecopcr_combined.RDS")) 
  input_data_before_multihit <- read_rds(file.path("database", db_name,  'inputs_before_multihit.RDS'))
  df_multihits <- read_rds(file.path("database", db_name, "df_multihtitlist.RDS"))
  inputs_used <- read_rds(file.path("database", db_name, 'inputs_sent_to_obitools.RDS'))
}

input_data_before_multihit <- read_rds(file.path("database", db_name,  'inputs_before_multihit.RDS'))
df_conflicts <- genereer_conflicten(ecopcr_combined, df_soortenlijst, df_ok_merges)
write_excel_csv2(df_conflicts, file = paste0(output_path,"/", "niet_op_soort_gebracht.csv"))

df_soortenevaluatie <- genereer_soortenevaluatie(ecopcr_combined, 
                                                 merged_data,
                                                 input_data_before_multihit,
                                                 df_soortenlijst, 
                                                 df_multihits, 
                                                 df_ok_merges, 
                                                 df_conflicts)
write_excel_csv2(df_soortenevaluatie, file = paste0(output_path,"/", "soortenevaluatie_riaz.csv"))



################################################################################
# 
# 
# inputs_used <- read_rds(file.path("database", db_name, 'inputs_sent_to_obitools.RDS'))
# test <- inputs_used %>% 
#   mutate(taxid = as.numeric(taxid)) %>% 
#   left_join(ecopcr_combined %>% mutate(is_amplified = TRUE),
#             by = c('genbank_id', 'taxid')) %>% 
#   mutate(is_amplified = ifelse(is.na(is_amplified), FALSE, is_amplified))
# inputampli <- test %>% 
#   group_by(source) %>% 
#   summarise(geamplificeerd = sum(is_amplified), niet_geamplificeerd = sum(!is_amplified))
# write_excel_csv2(inputampli, file.path("database", db_name, "output", "amplified_inputs_overview.csv"))
# 
# 
# 
# 
# #overzicht van records die niet op soort gebracht konden worden
# #indien 0 rijen dan is alles op soort gebracht of een toegelaten genusmerge of toegelaten familymerge
# ecopcr_filtered_sp <- ecopcr_combined %>% 
#   inner_join(df_soortenlijst %>% select(taxid, priority)) %>% 
#   select(genbank_id, taxid, rank, priority, species_name,
#          dna_hash, obi_rank, obi_taxid, genus_name, family_name,
#          merged_overview, obi_count) %>% 
#   filter(!(obi_rank %in% c("subspecies", "species", "species;"))) %>% 
#   group_by(across(-genbank_id)) %>% 
#   summarise(genbank_id = paste(genbank_id, collapse = ";"), .groups = "drop") %>% 
#   arrange(obi_taxid, taxid, dna_hash) %>% 
#   dplyr::filter(!(substring(obi_rank,1,5) == "genus" & 
#                     obi_taxid %in% (df_ok_merges %>% filter(RANK == "genus") %>% pull(TAXID))), 
#                 !(substring(obi_rank,1,6) == "family" & 
#                     obi_taxid %in% (df_ok_merges %>% filter(RANK == "family") %>% pull(TAXID))  ))
# write_excel_csv2(ecopcr_filtered_sp, file = paste0(output_path,"/", "conflicten_ruw.csv"))
# 
# #overzicht soorten
# df_beoordeeld <- 
#   ecopcr_filtered_sp %>%
#   arrange(obi_taxid) %>%
#   group_by(obi_taxid) %>%
#   do ({
#     df1 <- get_taxa_from_merged(., df_soortenlijst)
#     #genbank_ids <- paste(.$genbank_id, collapse = ";")
#     #bind_cols(df1, genbank_id = genbank_ids)
#     df1
#   }) %>% 
#   group_by(obi_taxid) %>% 
#   do({
#     judge_species(.)    
#   })
# df_conflicts <- df_beoordeeld %>% 
#   group_by(taxid) %>%
#   do({
#     taxid <- .$taxid[1]
#     obitaxids <- (unique(.$obi_taxid))
#     taxids <- unique(df_beoordeeld$taxid[df_beoordeeld$obi_taxid %in% obitaxids])
#     taxid_order <- order(taxids) 
#     scinams <- unique(df_beoordeeld$NameScientific[df_beoordeeld$obi_taxid %in% obitaxids])
#     data.frame(has_conflicts = TRUE, 
#                conflict_taxa = paste(taxids[taxid_order], collapse = '|'),
#                conflict_names = paste(scinams[taxid_order], collapse = '|'),
#                conflict_oordeel = paste(unique(.$oordeel), collapse = "|"))
#   }) %>% 
#   arrange(conflict_taxa)
# 
# write_excel_csv2(df_conflicts, file = paste0(output_path,"/", "niet_op_soort_gebracht.csv"))
# 
# ###################################################################
# ### SOORTENEVALUATIE
# ###################################################################
# 
# #soortentabel linken aan ecopcr output
# ecopcr_soorten <- ecopcr_combined %>% 
#   right_join(df_soortenlijst %>% select(taxid), by = "taxid") %>% 
#   group_by(taxid) %>% 
#   do({
#     merged = sum(.$is_merged)>0
#     obi_ranks = factor(.$obi_rank, levels = c("subspecies", "species", "subgenus", "genus", "subfamily", "family"))
#     whimaxrank = which.max(obi_ranks)
#     if (length(whimaxrank)) {
#       corresptaxid = as.numeric(.$obi_taxid[whimaxrank])
#       corresprank = as.character(.$obi_rank[whimaxrank])      
#     } else {
#       corresptaxid = as.numeric(.$obi_taxid[1])
#       corresprank = as.character(.$obi_rank[1])       
#     }
#     merged = .$taxid[1] != corresptaxid
#     rv <- data.frame(merged = merged, obi_taxid = corresptaxid, obi_rank = corresprank)
#     rv
#   })
# 
# #soorten-multihit
# multihit_species <- df_multihits %>% 
#   select(taxid = TAXID, pref_taxid = PREF_TAXID) %>% 
#   inner_join(df_soortenlijst %>% select(taxid)) %>% 
#   mutate(is_multihit = TRUE,
#          is_chosen = taxid == pref_taxid) 
# 
# multihit_species <- multihit_species %>% 
#   group_by(taxid) %>% 
#   do({
#     pref = .$pref_taxid[1]
#     taxids <- sort(multihit_species$taxid[pref == multihit_species$pref_taxid])
#     whi <- which(taxids == pref )
#     if(length(whi)) {
#       taxids = c(taxids[whi], taxids[-whi])
#     }
#     cbind(., hitlist = paste(taxids, collapse = " | "))
#   })
# 
# #ok_merges (namen komen dan hieruit ipv multihitlijst)
# df_merges <- df_ok_merges %>% 
#   left_join(ecopcr_combined %>%  select(obi_taxid, merged_overview),
#             by = c("TAXID" = "obi_taxid")) %>% 
#   distinct() %>% 
#   group_by(TAXID) %>% 
#   do({
#     taxids <- paste(.$merged_overview, collapse = ",")
#     taxids <- str_replace_all(taxids, "\\{", "")
#     taxids <- str_replace_all(taxids, "\\}", "")
#     taxids <- str_replace_all(taxids, "\\'", "")
#     separated <- unlist(str_split(taxids, ", "))
#     separated <- substring(separated, 1, regexpr(":", separated)-1)
#     combined <- paste(sort(as.numeric(unique(separated))), collapse = " | ")
#     data.frame(taxid_okm = .$TAXID[1], rank = .$RANK[1], name_okm = .$SCI_NAME[1], hitlist_okm = combined)
#   })
# 
# #voorkomen van soorten in de refdb als merged species zoek naar '#####' met ##### = taxid
# n_in_amplified <- df_soortenlijst %>% 
#   group_by(taxid) %>% 
#   do({
#     aantal <- sum(str_detect(merged_data$merged_taxid, paste0("'", .$taxid[1], "'")))
#     data.frame(n_in_amplicons = aantal)
#   })
# 
# #samenvattende tabel
# 
# df_soortenevaluatie <- df_soortenlijst %>% 
#   group_by(Group, NameScientific, NameEnglish, NameDutch, taxid, rank = Rank, priority) %>% 
#   summarise(n_input_orig = sum(input_data_before_multihit$taxid == .data$taxid[1]),
#             n_haplotypes = sum(merged_data$taxid == .data$taxid[1]),
#             .groups = "drop") %>% 
#   
#   left_join(n_in_amplified, by = "taxid") %>% 
#   
#   left_join(ecopcr_soorten, by = "taxid") %>% 
#   
#   left_join(multihit_species, by = "taxid") %>% 
#   mutate(pref_taxid = ifelse(is.na(pref_taxid) & (obi_taxid == taxid | obi_rank == "species"), 
#                              taxid, 
#                              pref_taxid)) %>%
#   mutate(is_multihit = ifelse(!(obi_rank %in% c("species", "subspecies")), TRUE, is_multihit)) %>% 
#   
#   left_join(df_soortenlijst %>% select(taxid, 
#                                        Pref_NameEnglish = NameEnglish, 
#                                        Pref_NameScientific = NameScientific), 
#             by = c("pref_taxid" = "taxid")) %>% 
#   
#   left_join(df_merges %>% transmute(taxid = TAXID, taxid_okm, name_okm, hitlist_okm), 
#             by = c("obi_taxid" = "taxid")) %>% 
#   mutate(Pref_NameScientific = ifelse(is.na(Pref_NameScientific), name_okm, Pref_NameScientific),
#          pref_taxid = ifelse(is.na(hitlist), taxid_okm, pref_taxid),
#          hitlist = ifelse(is.na(hitlist), hitlist_okm, hitlist)
#   ) %>% 
#   
#   left_join(df_conflicts, by = "taxid") %>% 
#   
#   select(Group, priority, n_input_orig, n_merged_taxid = n_in_amplicons, n_obi_taxid = n_haplotypes,
#          NameScientific, NameEnglish, NameDutch,
#          is_multihit, is_chosen, Pref_NameScientific, 
#          rank_obi_taxid = obi_rank, obi_taxid, rank, taxid, pref_taxid,
#          hitlist, Pref_NameEnglish, has_conflicts, conflict_taxa, conflict_oordeel, conflict_names) %>% 
#   arrange(Group, pref_taxid, desc(is_chosen))
# 
# write_excel_csv2(df_soortenevaluatie, file = paste0(output_path,"/", "soortenevaluatie_riaz.csv"))
# 
# 
#