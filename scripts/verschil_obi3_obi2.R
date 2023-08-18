
path <- "database/refdb_teleo_2023-08-10"

ecopcr3 <- parse_refdb_fasta(file.path(path, "amplified_clean.fasta"), is_merged_file = FALSE)
ecopcr2 <- parse_refdb_fasta(file.path(path, "amplified_clean2.fasta"), is_merged_file = FALSE)

ampli3r <- parse_refdb_fasta(file.path(path, "final_db_0.99.fasta"), is_merged_file = TRUE)
ampli3 <- parse_refdb_fasta(file.path(path, "amplified_clean_uniq.fasta"), is_merged_file = TRUE)
ampli2 <- parse_refdb_fasta(file.path(path, "amplified_clean_uniq2.fasta"), is_merged_file = TRUE)

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










