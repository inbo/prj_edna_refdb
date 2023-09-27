################################
### INPUTS OP ORDE STELLEN
################################

source("scripts/refdb_00_initialisation.R")

#!inputs inlezen (nadat de sequentiefouten eruit zijn)
df_inputs_cleaned <- read_rds(str_replace(cleaned_input_fasta, '.fasta', '.RDS'))
df_inputs_raw <- read_rds(str_replace(all_input_fasta, '.fasta', '.RDS'))

#! opnieuw toegestaan
df_passlist <- read_sheet(metadata_gdrive_key, "Passlist_Teleo")
df_passlist <- df_inputs_raw %>% 
  filter(genbank_id %in% df_passlist$ENTRY_ID) %>% 
  mutate(taxid = as.numeric(taxid))
nrow(df_passlist)


#! multihitsoorten
df_multihit <-  read_sheet(metadata_gdrive_key, "Multihitlist_Teleo", range = "A:H", 
                           col_types = 'ccccccnn') %>% 
  rename_all(., .funs = tolower)

#!toegelaten merges op hoger niveau
df_allowed_merges <- read_sheet(metadata_gdrive_key, "Toegelaten_merges_Teleo") %>% 
  rename_all(., .funs = tolower)

#! Pas de multihitlijst toe
#-----------------------------

df_inputs <- df_inputs_cleaned %>% 
  filter(!taxid %in% (df_multihit %>% 
                        filter(taxid != pref_taxid) %>% 
                        pull(taxid)))

#!  Herintroduceer sequenties uit de passlist
#---------------------------
df_inputs <- df_inputs %>% bind_rows(df_passlist) 

# remove full duplicates
df_inputs <- df_inputs %>% distinct(genbank_id, taxid, dna_hash, 
                                    .keep_all = TRUE)


## create fasta input for ecopcr (obitools2)
#---------------------------------------------

create_input_fasta(file = file.path("database", db_name_teleo, "input.fasta"), 
                   lowercase = TRUE, 
                   data = df_inputs)


#################################
### GENEREREN DB
#################################

   #
  # #
 # ! #   #Zorg dat docker desktop voor windows aan het draaien is 
#######   


#Bouw de image ("docker build -t obitools3pv ./docker_image_obitools3/")
system2("powershell", args = pscommand_build_teleo)#Run de container van de imag
#delete eerst docker rm obitools3testcontainer
#run: "docker run --rm -i -d --name obitools3testcontainer obitools3pv")
system2("powershell", args = pscommand_delete)
system2("powershell", args = pscommand_run)

##copy to container
system2('docker',  c('cp', 
                     file.path("database", db_name_teleo, obi_script_teleo), 
                     paste0(docker_container_name, ':', docker_script_teleo)))
system2('docker',  c('cp', taxdump_name, docker_taxonomy))
system2('docker',  c('cp', file.path("database", db_name_teleo, fasta_name), docker_input_fasta)) 

##execute main script
system2('docker', c("exec", docker_container_name, "chmod", "777", docker_script_teleo)) #???
system2('docker', c('exec', docker_container_name, 'bash', '-c', docker_script_teleo))

##return results from container to windows
system2('docker', c("cp", paste0(docker_path, "logfile.txt"), refdb_location_teleo))
system2('docker', c("cp", paste0(docker_path, obi_script_teleo), refdb_location_teleo))
system2('docker', c("cp", paste0(docker_path, "kept_input.fasta"), refdb_location_teleo))
system2('docker', c("cp", paste0(docker_path, "amplified.fasta"), refdb_location_teleo))
#system2('docker', c("cp", paste0(docker_path, "amplified_clean.fasta"), refdb_location_teleo))
system2('docker', c("cp", paste0(docker_path, "amplified_uniq.fasta"), refdb_location_teleo))
#system2('docker', c("cp", paste0(docker_path, "amplified_clean_uniq.fasta"), refdb_location_teleo))
system2('docker', c("cp", paste0(docker_path, "final_db_0.97.fasta"), refdb_location_teleo))
system2('docker', c("cp", paste0(docker_path, "final_db_0.99.fasta"), refdb_location_teleo))
system2('docker', c("cp", paste0(docker_path, "refdb.obidms"), refdb_location_teleo))        

##exit container
system2("powershell", paste("docker stop ", docker_container_name))       


#################################
### POSTPROCESSING
#################################

obitools_input_data <- parse_refdb_fasta(file.path("database", db_name_teleo, "kept_input.fasta"))
input_data_before_multihit <- read_rds(str_replace(cleaned_input_fasta, '.fasta', '.RDS'))

ecopcr_data <- 
  parse_refdb_fasta(file.path("database", db_name_teleo, "amplified.fasta"))

merged_data <- 
  parse_refdb_fasta(file.path("database", db_name_teleo, "final_db_0.99.fasta"))  %>% 
  mutate(obi_taxid = str_replace(substring(lca_taxid, 2), "]", ""),
         taxid = str_replace(taxid, ";", ""))

ecopcr_combined <- combine_ecpocr_with_merged(ecopcr_data, merged_data) 

df_soortenlijst <-  read_sheet(metadata_gdrive_key, "Soortenlijst") %>% 
  mutate(taxid = Taxid, priority = Priority, rank = Rank) %>% 
  filter(priority != 9)

###

df_conflicts <- genereer_conflicten(ecopcr_combined, df_soortenlijst, df_allowed_merges)
write_excel_csv2(df_conflicts, 
                 file = file.path('database', db_name_teleo,"niet_op_soort_gebracht.csv"))

df_soortenevaluatie <- genereer_soortenevaluatie(ecopcr_combined, 
                                                 merged_data,
                                                 input_data_before_multihit,
                                                 df_soortenlijst, 
                                                 df_multihit, 
                                                 df_allowed_merges, 
                                                 df_conflicts)
write_excel_csv2(df_soortenevaluatie, 
                 file = file.path('database', db_name_teleo, "soortenevaluatie_teleo.csv"))


