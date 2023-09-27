

df_soortenlijst_all <- read_sheet("1NidznDq9EVHN4_wfrhv0W5gZ0Q-S8wEn5jvXbb9rf8k", sheet = "Soortenlijst") %>% 
  filter(Taxid > 0) %>% 
  rename(taxid = Taxid, priority = Priority)

df_soortenlijst <- df_soortenlijst_all %>% 
  filter(priority %in% c(1:7))

merged2019orig <- 
  parse_refdb_fasta("input/Teleo_RefDB_aug2019_amplified_uniq.fasta", 
                    is_merged_file = TRUE)

merges <- merged2019orig$merged_taxid
merges <- str_replace(merges, "\\{", "{'")
merges <- str_replace(merges, ":", "':")
obitaxids <- paste0("'", merged2019orig$taxid, "'")


df_taxids <- df_soortenlijst %>% 
  select(taxid) %>%
  group_by(taxid) %>% 
  do({
    taxid = paste0("'", .$taxid[1], "'")
    count_obi <- sum(str_detect(obitaxids, taxid))
    count_amp <- sum(str_detect(merges, taxid))
    data.frame(n_merged_taxid_2019 = count_amp, n_obi_taxid_2019 = count_obi)
  })


df_soortenevaluatie_2019 <- read_csv2("database\\refdb_teleo_2023-09-01\\output\\soortenevaluatie_teleo.csv") %>% 
  left_join(df_taxids, by = "taxid")

write_excel_csv2(df_soortenevaluatie_2019, "output/soortenevaluatie_teleo_compared_2019.csv")

plot(jitter(df_soortenevaluatie_2019$n_merged_taxid), 
     jitter(df_soortenevaluatie_2019$n_merged_taxid_2019))
abline(0,1)

plot(jitter(df_soortenevaluatie_2019$n_obi_taxid), 
     jitter(df_soortenevaluatie_2019$n_obi_taxid_2019))
abline(0,1)

ggplot(df_soortenevaluatie_2019 %>% filter(priority == 1), 
       aes(x = jitter(n_merged_taxid), y = jitter(n_merged_taxid_2019))) +
  geom_point(pch = 1) + geom_abline(slope = 1)


