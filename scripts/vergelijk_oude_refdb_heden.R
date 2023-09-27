source("scripts/refdb_00_initialisation.R")

soortenlijst <- read_sheet(metadata_gdrive_key, "Soortenlijst") %>% 
  filter(Priority != 9) %>% 
  rename(taxid = Taxid, priority = Priority, rank = Rank)

##########################################
### RIAZ
##########################################

uniq_riaz_2022 <- "database\\refdb_riaz_2022-05-17\\amplified_clean_uniq.fasta"
uniq_riaz_2023 <- "database\\refdb_riaz_2023-09-06\\amplified_uniq2.fasta"

uniq22 <- parse_refdb_fasta(uniq_riaz_2022) %>% 
  select(genbank_id, merged_taxid, rank, taxid, dna_hash) %>% 
  mutate(merged_taxid = str_replace(merged_taxid, "\\{", "{'"),
         merged_taxid = str_replace_all(merged_taxid, ":", "':"),
         merged_taxid = str_replace_all(merged_taxid, ", ", ", '"))

uniq23 <- parse_refdb_fasta(uniq_riaz_2023) %>% 
  mutate(merged_taxid = str_replace(merged_taxid, "\\{", "{'"),
         merged_taxid = str_replace_all(merged_taxid, ":", "':"),
         merged_taxid = str_replace_all(merged_taxid, ", ", ", '"))

multihits <- read_sheet(metadata_gdrive_key, "Multihitlist_Riaz") %>% 
  rename_all(., .funs = tolower)

outputtabel <- soortenlijst %>% 
  group_by(NameScientific, NameEnglish, NameDutch, taxid, priority, rank) %>% 
  do({
    
    found22 <- sum(uniq22$taxid == .$taxid[1])
    found23 <- sum(uniq23$taxid == .$taxid[1])
    
    mh_chosen <- sum(.$taxid[1] %in% multihits$pref_taxid )
    mh_not_chosen <- sum(.$taxid[1] %in% multihits$taxid & !(.$taxid[1] %in% multihits$pref_taxid ))
    mh_status <- ifelse(mh_chosen, 
                        "multihit_chosen", 
                        ifelse(mh_not_chosen, 
                               "multihit_not_chosen", 
                               "no_multihit"))
    
  if (mh_status %in% c("multihit_chosen", "multihit_not_chosen")) {
    status22 <- status23 <- mh_status
  } else {
    status22 <- ifelse(found22, "determined", "not in db") 
    status23 <- ifelse(found23, "determined", "not in db") 
  }
    haplo22 <- sum(str_detect(uniq22$merged_taxid, paste0("'", .$taxid[1], "'")))
    haplo23 <- sum(str_detect(uniq23$merged_taxid, paste0("'", .$taxid[1], "'")))
    verschil <- haplo23 - haplo22
    data.frame(status22 = status22, status23 = status23, 
               n_merged_taxid_22 = haplo22, n_merged_taxid_23 = haplo23, 
               verschil_n_merged = verschil)
  })

table(outputtabel$status23, outputtabel$status22, outputtabel$priority)  
#verdwenen soorten (genus --> species)
outputtabel %>% 
  mutate(taxid = as.character(taxid)) %>% 
  filter(status23 == "not in db", status22 == "determined") %>% 
  arrange(priority) %>% view()

#gewonnen soorten
outputtabel %>% 
  mutate(taxid = as.character(taxid)) %>% 
  filter(status23 == "determined", status22 == "not in db") %>% 
  arrange(priority) %>% view()

##########################################
### TELEO
##########################################

uniq_teleo_2019 <- "database\\refdb_teleo_2019-08_Teleo_Aug2019\\Teleo_RefDB_aug2019_amplified_uniq.fasta"
uniq_teleo_2023 <- "database\\refdb_teleo_2023-09-06\\amplified_uniq2.fasta"

uniqteleo19 <- parse_refdb_fasta(uniq_teleo_2019) %>% 
  select(genbank_id, merged_taxid, rank, taxid, dna_hash) %>% 
  mutate(merged_taxid = str_replace(merged_taxid, "\\{", "{'"),
         merged_taxid = str_replace_all(merged_taxid, ":", "':"),
         merged_taxid = str_replace_all(merged_taxid, ", ", ", '"))

uniqteleo23 <- parse_refdb_fasta(uniq_teleo_2023) %>% 
  mutate(merged_taxid = str_replace(merged_taxid, "\\{", "{'"),
         merged_taxid = str_replace_all(merged_taxid, ":", "':"),
         merged_taxid = str_replace_all(merged_taxid, ", ", ", '"))

multihits_teleo <- read_sheet(metadata_gdrive_key, "Multihitlist_Teleo") %>% 
  rename_all(., .funs = tolower)

outputteleo <- soortenlijst %>% 
  group_by(NameScientific, NameEnglish, NameDutch, taxid, priority, rank) %>% 
  do({
    
    found19 <- sum(uniqteleo19$taxid == .$taxid[1])
    found23 <- sum(uniqteleo23$taxid == .$taxid[1])
    
    mh_chosen <- sum(.$taxid[1] %in% multihits_teleo$pref_taxid )
    mh_not_chosen <- sum(.$taxid[1] %in% multihits_teleo$taxid & 
                           !(.$taxid[1] %in% multihits_teleo$pref_taxid ))
    mh_status <- ifelse(mh_chosen, 
                        "multihit_chosen", 
                        ifelse(mh_not_chosen, 
                               "multihit_not_chosen", 
                               "no_multihit"))
    
    if (mh_status %in% c("multihit_chosen", "multihit_not_chosen")) {
      status19 <- status23 <- mh_status
    } else {
      status19 <- ifelse(found19, "determined", "not in db") 
      status23 <- ifelse(found23, "determined", "not in db") 
    }
    haplo19 <- sum(str_detect(uniqteleo19$merged_taxid, paste0("'", .$taxid[1], "'")))
    haplo23 <- sum(str_detect(uniqteleo23$merged_taxid, paste0("'", .$taxid[1], "'")))
    verschil <- haplo23 - haplo19
    data.frame(status19 = status19, status23 = status23, 
               n_merged_taxid_19 = haplo19, n_merged_taxid_23 = haplo23, 
               verschil_n_merged = verschil)
  })

table(outputteleo$status23, outputteleo$status19, outputteleo$priority)  
#verdwenen soorten
outputteleo %>% 
  mutate(taxid = as.character(taxid)) %>% 
  filter(status23 == "not in db", status19 == "determined") %>% 
  arrange(priority) %>% view()

#gewonnen soorten
outputteleo %>% 
  mutate(taxid = as.character(taxid)) %>% 
  filter(status23 == "determined", status19 == "not in db") %>% 
  arrange(priority) %>% view()


