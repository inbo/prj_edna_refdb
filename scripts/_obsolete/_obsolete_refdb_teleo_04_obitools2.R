#ZORG DAT DOCKER DESKTOP RUNT en OBITOOLS geladen is
#run de container obitools_pv 

#################
## INIT
#################

source(file.path("scripts", "refdb_teleo_00_initialisation.R"))

docker_container <- "obitools_pv"
fasta_name <- "input.fasta"
db_namepath <- paste0(db_name, "/", db_name) 

e <- 5
l <- 20
L <- 100
p1 <- select_primer_tags("TELEO")["oligo1"]
p2 <- select_primer_tags("TELEO")['oligo2']

#######################
### DOCKER PROCESSING
#######################

## Maak structuur aan en kopieer input.fasta
system2("docker", args = c("exec", docker_container, "mkdir", db_name))
system2("docker", args = c("cp", file.path("database", db_name, fasta_name), 
                           paste0(docker_container,":/",db_name)))

#KOPIEREN TAXONOMY indien nog niet gebeurd
system2("docker", args = c("exec", docker_container, "rm", "-r", "taxdump"))
system2("docker", args = c("exec", docker_container, "mkdir", "taxdump"))
system2("docker", args = c("cp", file.path(taxdump_name), 
                           paste0(docker_container,":/","taxdump/",  "taxdump.tar.gz")))
system2("docker", args = c("exec", docker_container, "tar", "-xzvf", "taxdump/taxdump.tar.gz", "-C", "taxdump"))

## converteer input+taxonomie naar voor obitools bruikbaar formaat

obiconvertstring <- paste("obiconvert", "--fasta", paste0("--ecopcrdb-output=",db_namepath),
                          "-t", "taxdump", 
                          paste0(db_name, "/", fasta_name))
system2('docker', args = c('exec', docker_container, 'sh', '-c', paste0('"', obiconvertstring, '"')))

## voer de in silico PCR uit
#als ecoPCR niet gekend is ./home/ecopcr/src/ecoPCR
ecopcrstring <- paste("/home/ecopcr/src/ecoPCR", "-d", db_namepath, 
                      "-e", e, "-l", l, "-L", L, p1, p2, 
                      ">", paste0(db_name, "/", "amplified2.ecopcr"))
system2('docker', args = c('exec', docker_container, 'sh', '-c', paste0('"', ecopcrstring, '"')))

## clean:(hou enkel die over die een parent rank hebben die soort en familie zijn)

#cleaning stap origineel
# grepstring <- paste("obigrep", "-d", db_namepath, 
#                     "--require-rank=species",
#                     "--require-rank=family",
#                     paste0(db_name, "/", "amplified2.ecopcr"),
#                     ">", 
#                     paste0(db_name, "/", "amplified_clean2.fasta"))
# system2('docker', args = c('exec', docker_container, 'sh', '-c', paste0('"', grepstring, '"')))

#non cleaning stap:DUMMY CLEANING
grepstring <- paste("obigrep", "-d", db_namepath,
                    paste0(db_name, "/", "amplified2.ecopcr"),
                    ">",
                    paste0(db_name, "/", "amplified_clean2.fasta"))
system2('docker', args = c('exec', docker_container, 'sh', '-c', paste0('"', grepstring, '"')))

## Keep unique sequences (merge identical amplicon sequences)

#UNIQ INDIEN GECLEAND OF DUMMY GECLEAND
uniqstring <- paste("obiuniq", "-d", db_namepath,
                    paste0(db_name, "/", "amplified_clean2.fasta"),
                    ">",
                    paste0(db_name, "/", "amplified_clean_uniq2.fasta"))
system2('docker', args = c('exec', docker_container, 'sh', '-c', paste0('"', uniqstring, '"')))

## Copy back to windows
system2("docker", args = c("cp", paste0(docker_container, ':', db_name),file.path("database")))


#############################
### POSTPROCESSING
#############################
library(tidyverse)
source(file.path("scripts", "refdb_teleo_00_initialisation.R"))
source(file.path("scripts", "_functions_postprocessing.R"))

ecopcr_file <- "amplified_clean2.fasta" #zelfde filename ook indien niet gecleand
merged_file <- "amplified_clean_uniq2.fasta"
input_fasta_file <- "input.fasta"

df_soortenlijst_all <- read_sheet(metadata_gdrive_key, sheet = "Soortenlijst") %>% 
  filter(Taxid > 0) %>% 
  rename(taxid = Taxid, priority = Priority)
df_soortenlijst <- df_soortenlijst_all %>% 
  filter(priority %in% c(1,2,3,4))
df_ok_merges <- read_sheet(metadata_gdrive_key, sheet = "Toegelaten_merges_Teleo") %>% 
  rename(taxid = TAXID, rank = RANK)
df_multihits <- read_sheet(metadata_gdrive_key, sheet = "Multihitlist_Teleo") %>% 
  rename(taxid = TAXID, pref_taxid = PREF_TAXID)

## ECOPCR INHOUD

input_data2 <- 
  parse_refdb_fasta(file.path("database", db_name, "kept_input.fasta"))

input_data_before_multihit <- readRDS(file.path("database", db_name, 'inputs_before_multihit.RDS'))

ecopcr_data2 <- 
  parse_refdb_fasta(file.path("database", db_name, ecopcr_file))

merged_data2 <- merged_data2_backup <- 
  parse_refdb_fasta(file.path("database", db_name, merged_file))
merged_data2 <- merged_data2 %>% 
  mutate(merged_taxid = str_replace(merged_taxid, "\\{", "{'"),
         merged_taxid = str_replace_all(merged_taxid, ":", "':"),
         merged_taxid = str_replace_all(merged_taxid, ", ", ", '"))

ecopcr_combined2 <- combine_ecpocr_with_merged(ecopcr_data2, merged_data2) 

table(ecopcr_combined2$is_merged)
table(ecopcr_combined2$obi_rank)
table(ecopcr_combined2$obi_rank, ecopcr_combined2$is_merged)


df_conflicts2 <- genereer_conflicten(ecopcr_combined2, df_soortenlijst, df_ok_merges)
write_excel_csv2(df_conflicts2, file = paste0(output_path,"/", "niet_op_soort_gebracht_obi2.csv"))

df_soortenevaluatie2 <- genereer_soortenevaluatie(ecopcr_combined2, 
                                                  merged_data2,
                                                  input_data_before_multihit,
                                                  df_soortenlijst, 
                                                  df_multihits, 
                                                  df_ok_merges, 
                                                  df_conflicts2)
write_excel_csv2(df_soortenevaluatie2, file = paste0(output_path,"/", "soortenevaluatie_teleo_obi2.csv"))

### CONFLICTEN EN SOORTENEVALUATIE
# 
# ecopcr_filtered_sp2 <- ecopcr_combined2 %>% 
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
# write_excel_csv2(ecopcr_filtered_sp2, file = paste0(output_path,"/", "conflicten_ruw2.csv"))
# 
# #overzicht soorten
# df_beoordeeld2 <- 
#   ecopcr_filtered_sp2 %>%
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
# df_conflicts2 <- df_beoordeeld2 %>% 
#   group_by(taxid) %>%
#   do({
#     taxid <- .$taxid[1]
#     obitaxids <- (unique(.$obi_taxid))
#     taxids <- unique(df_beoordeeld2$taxid[df_beoordeeld2$obi_taxid %in% obitaxids])
#     taxid_order <- order(taxids) 
#     scinams <- unique(df_beoordeeld2$NameScientific[df_beoordeeld2$obi_taxid %in% obitaxids])
#     data.frame(has_conflicts = TRUE, 
#                conflict_taxa = paste(taxids[taxid_order], collapse = '|'),
#                conflict_names = paste(scinams[taxid_order], collapse = '|'),
#                conflict_oordeel = paste(unique(.$oordeel), collapse = "|"))
#   }) %>% 
#   arrange(conflict_taxa)
# 
# write_excel_csv2(df_conflicts2, file = paste0(output_path,"/", "niet_op_soort_gebracht_obi2.csv"))
# 
# ###################################################################
# ### SOORTENEVALUATIE
# ###################################################################
# 
# #soortentabel linken aan ecopcr output
# ecopcr_soorten2 <- ecopcr_combined2 %>% 
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
# multihit_species2 <- df_multihits %>% 
#   select(taxid = TAXID, pref_taxid = PREF_TAXID) %>% 
#   inner_join(df_soortenlijst %>% select(taxid)) %>% 
#   mutate(is_multihit = TRUE,
#          is_chosen = taxid == pref_taxid) 
# 
# multihit_species2 <- multihit_species2 %>% 
#   group_by(taxid) %>% 
#   do({
#     pref = .$pref_taxid[1]
#     taxids <- sort(multihit_species2$taxid[pref == multihit_species2$pref_taxid])
#     whi <- which(taxids == pref )
#     if(length(whi)) {
#       taxids = c(taxids[whi], taxids[-whi])
#     }
#     cbind(., hitlist = paste(taxids, collapse = " | "))
#   })
# 
# #ok_merges (namen komen dan hieruit ipv multihitlijst)
# df_merges2 <- df_ok_merges %>% 
#   left_join(ecopcr_combined2 %>%  select(obi_taxid, merged_overview),
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
# n_in_amplified2 <- df_soortenlijst %>% 
#   group_by(taxid) %>% 
#   do({
#     aantal <- sum(str_detect(merged_data2$merged_taxid, paste0("'", .$taxid[1], "'")))
#     data.frame(n_in_amplicons = aantal)
#   })
# 
# #samenvattende tabel
# 
# df_soortenevaluatie2 <- df_soortenlijst %>% 
#   group_by(Group, NameScientific, NameEnglish, NameDutch, taxid, rank = Rank, priority) %>% 
#   summarise(n_input_orig = sum(input_data_before_multihit$taxid == .data$taxid[1]),
#             n_haplotypes = sum(merged_data2$taxid == .data$taxid[1]),
#             .groups = "drop") %>% 
#   
#   left_join(n_in_amplified2, by = "taxid") %>% 
#   
#   left_join(ecopcr_soorten2, by = "taxid") %>% 
#   
#   left_join(multihit_species2, by = "taxid") %>% 
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
#   left_join(df_merges2 %>% transmute(taxid = TAXID, taxid_okm, name_okm, hitlist_okm), 
#             by = c("obi_taxid" = "taxid")) %>% 
#   mutate(Pref_NameScientific = ifelse(is.na(Pref_NameScientific), name_okm, Pref_NameScientific),
#          pref_taxid = ifelse(is.na(hitlist), taxid_okm, pref_taxid),
#          hitlist = ifelse(is.na(hitlist), hitlist_okm, hitlist)
#   ) %>% 
#   
#   left_join(df_conflicts2, by = "taxid") %>% 
#   
#   select(Group, priority, n_input_orig, n_merged_taxid = n_in_amplicons, n_obi_taxid = n_haplotypes,
#          NameScientific, NameEnglish, NameDutch,
#          is_multihit, is_chosen, Pref_NameScientific, 
#          rank_obi_taxid = obi_rank, obi_taxid, rank, taxid, pref_taxid,
#          hitlist, Pref_NameEnglish, has_conflicts, conflict_taxa, conflict_oordeel, conflict_names) %>% 
#   arrange(Group, pref_taxid, desc(is_chosen))
# 
# write_excel_csv2(df_soortenevaluatie2, file = paste0(output_path,"/", "soortenevaluatie_teleo_obi2.csv"))



# ecopcr_filtered_sp2 <- ecopcr_combined2 %>% 
#   inner_join(df_soortenlijst %>% select(taxid, priority)) %>% 
#   select(genbank_id, taxid, rank, priority, species_name,
#          amplicon_hash, obi_rank, obi_taxid, genus_name, family_name,
#          merged_overview, obi_count) %>% 
#   filter(!(obi_rank %in% c("subspecies", "species"))) %>% 
#   group_by(across(-genbank_id)) %>% 
#   summarise(genbank_id = paste(genbank_id, collapse = ";")) %>% 
#   arrange(obi_taxid, taxid, amplicon_hash) %>% 
#   dplyr::filter(!(obi_rank == "genus" & 
#                     obi_taxid %in% (df_ok_merges %>% filter(RANK == "genus") %>% pull(TAXID))), 
#                 !(obi_rank == "family" & 
#                     obi_taxid %in% (df_ok_merges %>% filter(RANK == "family") %>% pull(TAXID))  ))
# 
# df_spm2 <- ecopcr_filtered_sp2 %>%
#   arrange(obi_taxid) %>%
#   group_by(obi_taxid) %>%
#   do ({
#     get_taxa_from_merged(., df_soortenlijst)
#   }
#   )
# 
# df_beoordeeld2 <- df_spm2 %>% group_by(obi_taxid) %>%
#   do ({
#     judge_species(.)
#   })
# 
# beoordeling2 <- df_beoordeeld2 %>% group_by(taxid, pref_taxid) %>%
#   summarise(oordeel = paste(unique(oordeel), collapse = "|"))
# 
# df_soortentabel2 <- df_soortenlijst %>% 
#   group_by(Group, NameScientific, NameEnglish, NameDutch, taxid, priority) %>% 
#   summarise(n_input_orig = sum(input_data_before_multihit$taxid == .data$taxid[1]),
#             #n_input_seqs = sum(as.numeric(input_data$TAXID) == .data$taxid[1]),
#             n_input_seqs = sum(as.numeric(input_data2$TAXID == .data$taxid[1]) > 0, na.rm = TRUE),
#             n_amplicons = sum(ecopcr_combined2$taxid == .data$taxid[1]),
#             n_determined = sum(ecopcr_combined2$obi_taxid == .data$taxid[1]))
# 
# ecopcr_soorten2 <- ecopcr_combined2 %>% 
#   right_join(df_soortenlijst %>% select(taxid)) %>% 
#   group_by(taxid) %>% 
#   summarise(merged = sum(is_merged)>0, 
#             obi_taxid = paste(unique(obi_taxid), collapse = '|'),
#             obi_rank = paste(unique(obi_rank), collapse = "|")) %>% 
#   mutate(merged = merged & taxid != obi_taxid)
# 
# df_multihits2 <-
#   df_beoordeeld2 %>% group_by(taxid) %>%
#   do({
#     taxid <- .$taxid[1]
#     obitaxids <- (unique(.$obi_taxid))
#     taxids <- df_beoordeeld2$taxid[df_beoordeeld2$obi_taxid %in% obitaxids]
#     data.frame(multihit_taxa = paste("|", sort(unique(taxids)), collapse = '', sep = ''))
#   })
# 
# df_soortenevaluatie2 <- df_soortentabel2 %>% 
#   left_join(ecopcr_soorten2) %>% 
#   left_join(beoordeling2) %>% 
#   left_join(df_multihits2) %>% 
#   mutate(pref_taxid = ifelse(merged == FALSE, taxid, pref_taxid)) %>% #klopt niet, merged is fout
#   left_join(df_soortenlijst %>% select(taxid, NameEnglish, NameScientific), by = c("pref_taxid" = "taxid"))
# 
# write_excel_csv2(df_soortenevaluatie2, file = paste0(output_path,"/", "soortenevaluatie_obi2_ruw.csv"))
# 
# 
