
path <- "database/refdb_teleo_2023-09-01"

# amplificaties

ampli2 <- read_delim(file.path(path, "amplified2.ecopcr"), skip = 13, col_names = FALSE, delim = " | ") %>% 
  select(c(1,3,4))
colnames(ampli2) <- c("genbank_id", "taxid", "rank") 
ampli2 <- ampli2 %>% transmute(genbank_id = trimws(genbank_id))

ampli3 <- parse_refdb_fasta(file.path(path, "amplified.fasta")) %>% 
  select(genbank_id, taxid = TAXID, rank)

table(substring(ampli2$genbank_id, 1, 20) %in% substring(ampli3$genbank_id, 1, 20))
table(substring(ampli3$genbank_id, 1, 20) %in% substring(ampli2$genbank_id, 1, 20))

# cleaning

ecopcr3 <- parse_refdb_fasta(file.path(path, "amplified_clean.fasta"), is_merged_file = FALSE)
ecopcr2 <- parse_refdb_fasta(file.path(path, "amplified_clean2.fasta"), is_merged_file = FALSE)

table(substring(ecopcr2$genbank_id, 1, 20) %in% substring(ecopcr3$genbank_id, 1, 20))
table(substring(ecopcr3$genbank_id, 1, 20) %in% substring(ecopcr2$genbank_id, 1, 20))

# uniques

ampli3r <- parse_refdb_fasta(file.path(path, "final_db_0.99.fasta"), is_merged_file = TRUE)
ampli3 <- parse_refdb_fasta(file.path(path, "amplified_clean_uniq.fasta"), is_merged_file = TRUE)
ampli2 <- parse_refdb_fasta(file.path(path, "amplified_clean_uniq2.fasta"), is_merged_file = TRUE)

df_compare <- data.frame(DNA_HASH = unique(c(ampli3$DNA_HASH, ampli2$DNA_HASH)))
df_compare <- df_compare %>% 
  left_join(ampli2 %>% select(DNA_HASH, taxid, merged_taxid, rank)) %>% 
  left_join(ampli3 %>% transmute(DNA_HASH, 
                                 taxid3 = as.numeric(str_replace(TAXID, ";", "")),
                                 merged_taxid3 = MERGED_TAXID, 
                                 rank3 = rank)) %>% 
  mutate(ranks_equal = rank == rank3,
         taxid_equal = taxid == taxid3)

table(df_compare$ranks_equal)
table(df_compare$taxid_equal)


table(ampli2$DNA_HASH %in% ampli3$DNA_HASH)
table(ampli3$DNA_HASH %in% ampli2$DNA_HASH)


not_in_ecopcr3 <- ecopcr2 %>% 
  slice(which(!(ecopcr2$genlab_id %in% ecopcr3$genlab_id))) %>% 
  select(genlab_id, taxid, family_name, species_name, rank) %>% 
  mutate("niet in obitools3")

not_in_ecopcr2 <- ecopcr3 %>% 
  slice(which(!(ecopcr3$genlab_id %in% ecopcr2$genlab_id))) %>% 
  select(genlab_id, taxid = TAXID, family_name, species_name, rank) %>% 
  mutate("niet in obitools2")

bind_rows(not_in_ecopcr2, not_in_ecopcr3) %>% 
  write_excel_csv2(file.path(output_path, "verschillen_obitools2_en_3.csv"))




which(!(ampli2$genlab_id %in% ampli3$genlab_id))
which(!(ampli3$genlab_id %in% ampli2$genlab_id))

comparedata <- data.frame(genlab_id = unique(c(ampli3$genlab_id, ampli2$genlab_id)))
comparedata <- comparedata %>% 
  left_join(ampli2 %>% select(genlab_id, 
                              count2 = "merged_count",
                              merged_taxid2 = "merged_taxid",
                              obi_rank2 = "rank",
                              obi_taxid2 = "taxid") %>% 
              distinct()) %>% 
  left_join(ampli3 %>% select(genlab_id, 
                              count3 = "COUNT",
                              merged_taxid3 = "MERGED_TAXID",
                              obi_rank3 = "rank",
                              obi_taxid3 = "TAXID") %>% 
              distinct()) %>% 
  left_join(ampli3r %>% select(genlab_id, 
                               count3r = "COUNT",
                               merged_taxid3r = "MERGED_TAXID",
                               obi_rank3r = "rank",
                               obi_TXID = "TAXID",
                               obi_TXLCA = "LCA_TAXID")) %>% 
  mutate(merged_taxid3 = stringr::str_replace_all(merged_taxid3, "'", "") ,
         obi_TXLCA = stringr::str_replace(obi_TXLCA, "\\[", ""),
         obi_TXLCA = stringr::str_replace(obi_TXLCA, "\\]", "")) %>% 
  mutate(same_count = count2 == count3, 
         same_merge = merged_taxid2 == merged_taxid3,
         same_rank = obi_rank2 == obi_rank3,
         same_taxid = obi_taxid2 == obi_taxid3,
         count3r = count3 == count3r,
         rankr =  obi_rank3 == obi_rank3r,
         taxid33 = obi_TXID == obi_TXLCA,
         taxid23 = obi_taxid2 == obi_TXLCA, 
         rank2lca = obi_rank2 == obi_rank3r)

te_bekijken_seqs <- comparedata %>% filter(taxid33 == FALSE) 
write_excel_csv2(te_bekijken_seqs, file.path(output_path,"te_controleren_seq.csv"))


################################################################

path <- "database/refdb_teleo_2023-09-06"
soortenevaluatie2 <- read_csv2(file.path(path, "output", "soortenevaluatie_teleo_obi2.csv"))
soortenevaluatie3 <- read_csv2(file.path(path, "output", "soortenevaluatie_teleo.csv"))

soortenevaluatiecomp <- soortenevaluatie3 %>% distinct() %>% 
  left_join(soortenevaluatie2 %>% distinct() %>% 
              select(taxid, 
                     n_merged_taxid2 = n_merged_taxid, 
                     n_obi_taxid2 = n_obi_taxid)) %>% 
  mutate(merged_comp = n_merged_taxid - n_merged_taxid2,
         obi_comp = n_obi_taxid - n_obi_taxid2)






