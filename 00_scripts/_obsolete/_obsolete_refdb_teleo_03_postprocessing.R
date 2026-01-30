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

df_ok_merges <- read_sheet(metadata_gdrive_key, sheet = "Toegelaten_merges_Teleo") %>% 
  rename(taxid = TAXID, rank = RANK)

df_multihits <- read_sheet(metadata_gdrive_key, sheet = "Multihitlist_Teleo") %>% 
  rename(taxid = TAXID, pref_taxid = PREF_TAXID)

## ECOPCR INHOUD

input_data <- 
  parse_refdb_fasta(file.path("database", db_name, "kept_input.fasta"), 
                    is_merged_file = FALSE)

input_data_before_multihit <- readRDS(file.path("database", db_name, 'inputs_before_multihit.RDS'))

#INDIEN GECLEAND
# ecopcr_data <- 
#   parse_refdb_fasta(file.path("database", db_name, ecopcr_file), 
#                     is_merged_file = FALSE) %>% 
#   rename(taxid = TAXID,
#          dna_hash = DNA_HASH,
#          count = COUNT) %>% 
#   mutate(taxid = str_replace(taxid, ";", ""))

#INDIEN NIET GECLEAND
ecopcr_data <-
  parse_refdb_fasta(file.path("database", db_name, before_ecopcr_file),
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
              select(taxid, priority, NameScientific, NameEnglish))
mtch <- NULL
for (i in seq_len(nrow(niet_geamplificeerd)))
  mtch[i] <- find_matching_seq(niet_geamplificeerd$dna_sequence[i])[8]
table(mtch)
geweigerd_ondanks_match <- niet_geamplificeerd %>% 
  mutate(gematched = mtch) %>% 
  filter(mtch == 11)
view(geweigerd_ondanks_match)
write_excel_csv2(geweigerd_ondanks_match, 
                 file = file.path(output_path, "geweigerd_ondanks_perfecte_match.csv"))


inputs_used <- read_rds(file.path("database", db_name, 'inputs_sent_to_obitools.RDS'))
test <- inputs_used %>% 
  mutate(taxid = as.numeric(taxid)) %>% 
  left_join(ecopcr_combined %>% mutate(is_amplified = TRUE)) %>% 
  mutate(is_amplified = ifelse(is.na(is_amplified), FALSE, is_amplified))
inputampli <- test %>% 
  group_by(source) %>% 
  summarise(geamplificeerd = sum(is_amplified), niet_geamplificeerd = sum(!is_amplified))
write_excel_csv2(inputampli, file.path("database", db_name, "output", "amplified_inputs_overview.csv"))


###############################################################################
### CONFLICTEN OPSPOREN EN SOORTENTABEL AANMAKEN
###############################################################################

#importeren data

from_rds <- FALSE
if (from_rds) {
  library(tidyverse)
  source("scripts/_functions_fasta.R")
  source("scripts/_functions_postprocessing.R")
  source("scripts/refdb_teleo_00_initialisation.R")
  df_soortenlijst <- read_rds(file.path("database", db_name, "df_soortenlijst.RDS"))
  df_ok_merges <- read_rds(file.path("database", db_name, "df_ok_merges.RDS"))
  ecopcr_data <- read_rds(file.path("database", db_name, "ecopcr_data.RDS"))
  merged_data <- read_rds(file.path("database", db_name, "merged_data.RDS"))
  ecopcr_combined <- read_rds(file.path("database", db_name, "ecopcr_combined.RDS")) 
  input_data_before_multihit <- read_rds(file.path("database", db_name,  'inputs_before_multihit.RDS'))
  df_multihits <- read_rds(file.path("database", db_name, "df_multihtitlist.RDS"))
  #inputs_used <- read_rds(file.path("database", db_name, 'inputs_sent_to_obitools.RDS'))
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
write_excel_csv2(df_soortenevaluatie, file = paste0(output_path,"/", "soortenevaluatie_teleo.csv"))
