   #
  ###
 # ! #
#######

#ZORG DAT DOCKER DESKTOP RUNT en OBITOOLS geladen is
#run de container obitools_pv (image romunov/obitools:1.2.13)

#################
## INIT
#################

source(file.path("scripts", "refdb_00_initialisation.R"))

docker_container <- "obitools_pv"
fasta_name <- "input.fasta"
db_namepath_riaz <- paste0(db_name_riaz, "/", db_name_riaz) 

e <- 4
l <- 50
L <- 160
p1 <- select_primer_tags("RIAZ")["oligo1"]
p2 <- select_primer_tags("RIAZ")['oligo2']

#######################
### DOCKER PROCESSING
#######################

## Maak structuur aan en kopieer input.fasta
system2("docker", args = c("exec", docker_container, "mkdir", db_name_riaz))
system2("docker", args = c("cp", file.path("database", db_name_riaz, fasta_name), 
                           paste0(docker_container,":/",db_name_riaz)))

#KOPIEREN TAXONOMY indien nog niet gebeurd
system2("docker", args = c("exec", docker_container, "rm", "-r", "taxdump"))
system2("docker", args = c("exec", docker_container, "mkdir", "taxdump"))
system2("docker", args = c("cp", file.path(taxdump_name), 
                           paste0(docker_container,":/","taxdump/",  "taxdump.tar.gz")))
system2("docker", args = c("exec", docker_container, "tar", "-xzvf", "taxdump/taxdump.tar.gz", "-C", "taxdump"))

## ecopcrdb: converteer input+taxonomie naar voor obitools bruikbaar formaat
obiconvertstring <- paste("obiconvert", "--fasta", paste0("--ecopcrdb-output=",db_namepath_riaz),
                          "-t", "taxdump", 
                          paste0(db_name_riaz, "/", fasta_name))
system2('docker', args = c('exec', docker_container, 'sh', '-c', paste0('"', obiconvertstring, '"')))

## voer de in silico PCR uit( ./home/ecopcr/src/ecoPCR)
ecopcrstring <- paste("/home/ecopcr/src/ecoPCR", "-d", db_namepath_riaz, 
                      "-e", e, "-l", l, "-L", L, p1, p2, 
                      ">", paste0(db_name_riaz, "/", "amplified2.ecopcr"))
system2('docker', args = c('exec', docker_container, 'sh', '-c', paste0('"', ecopcrstring, '"')))

## convert ecopcr naar fasta
amplifiedfastastring <- paste("obiconvert","--ecopcr", "--fasta-output", 
                              paste0(db_name_riaz, "/","amplified2.ecopcr"),
                              ">", paste0(db_name_riaz, "/","amplified2.fasta"))
system2('docker', args = c('exec', docker_container, 'sh', '-c', paste0('"', amplifiedfastastring, '"')))

## extra annotatiestap (nodig voor Annelies om foutmelding te vermijden)
extraannotatestring <- paste("obiannotate", "--uniq-id", 
                             paste0(db_name_riaz, "/","amplified2.fasta"), 
                             ">", paste0(db_name_riaz, "/","amplified2_extrastep.fasta"))
system2('docker', args = c('exec', docker_container, 'sh', '-c', paste0('"', extraannotatestring, '"')))


## unique sequences (merge identical amplicon sequences)

uniqstring <- paste("obiuniq", "-d", db_namepath_riaz, 
                    paste0(db_name_riaz, "/", "amplified2_extrastep.fasta"),
                    ">", paste0(db_name_riaz, "/", "amplified_uniq2.fasta"))
system2('docker', args = c('exec', docker_container, 'sh', '-c', paste0('"', uniqstring, '"')))

## Copy back to windows
system2("docker", args = c("cp", paste0(docker_container, ':', db_name_riaz),file.path("database")))



#############################
### POSTPROCESSING
#############################


obitools_input_data <- parse_refdb_fasta(file.path("database", db_name_riaz, "kept_input.fasta"))
input_data_before_multihit <- read_rds(str_replace(cleaned_input_fasta, '.fasta', '.RDS'))

ecopcr_data2 <- 
  parse_refdb_fasta(file.path("database", db_name_riaz, "amplified2.fasta"))

merged_data2 <- merged_data2_backup <- 
  parse_refdb_fasta(file.path("database", db_name_riaz, "amplified_uniq2.fasta"))
merged_data2 <- merged_data2 %>% 
  mutate(merged_taxid = str_replace(merged_taxid, "\\{", "{'"),
         merged_taxid = str_replace_all(merged_taxid, ":", "':"),
         merged_taxid = str_replace_all(merged_taxid, ", ", ", '"))

ecopcr_combined2 <- combine_ecpocr_with_merged(ecopcr_data2, merged_data2) 

df_soortenlijst <-  read_sheet(metadata_gdrive_key, "Soortenlijst") %>% 
  mutate(taxid = Taxid, priority = Priority, rank = Rank) %>% 
  filter(priority != 9)

#! multihitsoorten
df_multihit <-  read_sheet(metadata_gdrive_key, "Multihitlist_Riaz", range = "A:H", 
                           col_types = 'ccccccnn') %>% 
  rename_all(., .funs = tolower)

#!toegelaten merges op hoger niveau
df_allowed_merges <- read_sheet(metadata_gdrive_key, "Toegelaten_merges_Riaz") %>% 
  rename_all(., .funs = tolower)


###

df_conflicts2 <- genereer_conflicten(ecopcr_combined2, df_soortenlijst, df_allowed_merges)
write_excel_csv2(df_conflicts2, 
                 file = file.path('database', db_name_riaz,"niet_op_soort_gebracht_obi2.csv"))

df_soortenevaluatie2 <- genereer_soortenevaluatie(ecopcr_combined2, 
                                                 merged_data2,
                                                 input_data_before_multihit,
                                                 df_soortenlijst, 
                                                 df_multihit, 
                                                 df_allowed_merges, 
                                                 df_conflicts2)
write_excel_csv2(df_soortenevaluatie2, 
                 file = file.path('database', db_name_riaz, "soortenevaluatie_riaz_obi2.csv"))

#--------------------------------------------------------------------------







# 
# 
# 
# 
# 
# 
# ecopcr_file <- "amplified_clean2.fasta" #zelfde filename ook indien niet gecleand
# merged_file <- "amplified_clean_uniq2.fasta"
# input_fasta_file <- "input.fasta"
# 
# df_soortenlijst_all <- read_sheet(metadata_gdrive_key, sheet = "Soortenlijst") %>% 
#   filter(Taxid > 0) %>% 
#   rename(taxid = Taxid, priority = Priority)
# df_soortenlijst <- df_soortenlijst_all %>% 
#   filter(priority %in% c(1,2,3,4))
# df_ok_merges <- read_sheet(metadata_gdrive_key, sheet = "Toegelaten_merges_Riaz") %>% 
#   rename(taxid = TAXID, rank = RANK)
# df_multihits <- read_sheet(metadata_gdrive_key, sheet = "Multihitlist_Riaz") %>% 
#   rename(taxid = TAXID, pref_taxid = PREF_TAXID)
# 
# ## ECOPCR INHOUD
# 
# input_data2 <- 
#   parse_refdb_fasta(file.path("database", db_name, "kept_input.fasta"), 
#                     is_merged_file = FALSE)
# 
# input_data_before_multihit <- readRDS(file.path("database", db_name, 'inputs_before_multihit.RDS'))
# 
# ecopcr_data2 <- 
#   parse_refdb_fasta(file.path("database", db_name, ecopcr_file), 
#                     is_merged_file = FALSE) %>% 
#   rename(taxid = taxid,
#          dna_hash = DNA_HASH)
# 
# merged_data2 <- merged_data2_backup <- 
#   parse_refdb_fasta(file.path("database", db_name, merged_file), 
#                     is_merged_file = TRUE) %>% 
#   rename(merged_count = merged_count, 
#          dna_hash = DNA_HASH,
#          merged_taxid = merged_taxid,
#          taxid = taxid)
# merged_data2 <- merged_data2 %>% 
#   mutate(merged_taxid = str_replace(merged_taxid, "\\{", "{'"),
#          merged_taxid = str_replace_all(merged_taxid, ":", "':"),
#          merged_taxid = str_replace_all(merged_taxid, ", ", ", '"))
# 
# ecopcr_combined2 <- combine_ecpocr_with_merged(ecopcr_data2, merged_data2) 
# 
# table(ecopcr_combined2$is_merged)
# table(ecopcr_combined2$obi_rank)
# table(ecopcr_combined2$obi_rank, ecopcr_combined2$is_merged)
# 
# 
# df_conflicts2 <- genereer_conflicten(ecopcr_combined2, df_soortenlijst, df_ok_merges)
# write_excel_csv2(df_conflicts2, file = paste0(output_path,"/", "niet_op_soort_gebracht_obi2.csv"))
# 
# df_soortenevaluatie2 <- genereer_soortenevaluatie(ecopcr_combined2, 
#                                                   merged_data2,
#                                                   input_data_before_multihit,
#                                                   df_soortenlijst, 
#                                                   df_multihits, 
#                                                   df_ok_merges, 
#                                                   df_conflicts2)
# write_excel_csv2(df_soortenevaluatie2, file = paste0(output_path,"/", "soortenevaluatie_riaz_obi2.csv"))
# 

