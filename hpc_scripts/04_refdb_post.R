#################################
### POSTPROCESSING
#################################

obitools_input_data <- parse_refdb_fasta(file.path("database", db_name_riaz, "kept_input.fasta"))
input_data_before_multihit <- read_rds(str_replace(cleaned_input_fasta, '.fasta', '.RDS'))

ecopcr_data <- 
  parse_refdb_fasta(file.path("database", db_name_riaz, "amplified.fasta"))

merged_data <- 
  parse_refdb_fasta(file.path("database", db_name_riaz, "final_db_0.99.fasta"))  %>% 
  mutate(obi_taxid = str_replace(substring(lca_taxid, 2), "]", ""),
         taxid = str_replace(taxid, ";", ""))

ecopcr_combined <- combine_ecpocr_with_merged(ecopcr_data, merged_data) 

df_soortenlijst <-  read_sheet(metadata_gdrive_key, "Soortenlijst") %>% 
  mutate(taxid = Taxid, priority = Priority, rank = Rank) %>% 
  filter(priority != 9)

###

df_conflicts <- genereer_conflicten(ecopcr_combined, df_soortenlijst, df_allowed_merges)
write_excel_csv2(df_conflicts, 
                 file = file.path('database', db_name_riaz,"niet_op_soort_gebracht.csv"))

df_soortenevaluatie <- genereer_soortenevaluatie(ecopcr_combined, 
                                                 merged_data,
                                                 input_data_before_multihit,
                                                 df_soortenlijst, 
                                                 df_multihit, 
                                                 df_allowed_merges, 
                                                 df_conflicts)
write_excel_csv2(df_soortenevaluatie, 
                 file = file.path('database', db_name_riaz, "soortenevaluatie_riaz.csv"))
