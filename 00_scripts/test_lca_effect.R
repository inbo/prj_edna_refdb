
source("scripts/_functions_fasta.R")
library(tidyverse)
library(googlesheets4)

splist <- "1NidznDq9EVHN4_wfrhv0W5gZ0Q-S8wEn5jvXbb9rf8k"

file99 <- "database/refdb_teleo_2023-07-31/final_db_0.99.fasta"
file97 <- "database/refdb_teleo_2023-07-31/final_db_0.97.fasta"
file96 <- "database/refdb_teleo_2023-07-31/final_db_0.96.fasta"
file10 <- "database/refdb_teleo_2023-07-31/final_db_1.fasta"

#LCA score cutoff op 0.99 of 1.0 is identiek, omdat er niets tussen 0.99 en 1 ligt omdat een teleo fragment maximaal 100 baseparen telt (bij Riaz is dat wel mogelijk). 2 LCA_TAXID waren niet gelijk, maar het was wel dezelfde soort
#Bij lagere LCA scores zijn meerdere matches mogelijk: dan krijg je meerdere LCA_scores en meerdere LCA_Taxid voor 1 record
#Aangezien bij Annelies er geen meerdere matches waren, werkte obitools2 waarschijnlijk niet met deze scores en enkel met exacte matches, wat we hier gaan aanhouden. Alhoewel ze werkte wel met een cutoff van 0.96, maar dat is mogelijks een heel ander getal.
#Waarom wijzigt LCA score, een score vqn 0.98,0.97 in de 0.96 fasta is 1.0 in de 0.99 fasta


test99  <- parse_refdb_fasta(file99)
test97 <- parse_refdb_fasta(file97)
test10  <- parse_refdb_fasta(file10)
df_species_all <- read_sheet(splist, "Soortenlijst") %>% filter(Taxid > 0)


test97 <- test97 %>% 
  mutate(LCA_N_TAXID = stringr::str_count(LCA_TAXID, ",") + 1)

test97b <- test97 %>% 
  filter(LCA_N_TAXID > 1) %>% 
  transmute(TAXID = as.numeric(TAXID), LCA_TAXID, LCA_SCORE, species, genus, family, MERGED_TAXID, LCA_N_TAXID) %>% 
  arrange(LCA_TAXID) %>% 
  distinct() %>% 
  inner_join(df_species_all %>% select(Taxid, NameScientific, NameEnglish, Priority ), by = c("TAXID"="Taxid"))




test <- bind_rows(test96,test99)
table(duplicated(test %>% select(genlab_id)))

table(test99$LCA_SCORE == test96$LCA_SCORE)
table(test99$LCA_SCORE == test10$LCA_SCORE)
table(test99$TAXID == test96$TAXID)
table(test99$LCA_TAXID == test96$LCA_TAXID)
table(test99$LCA_TAXID == test10$LCA_TAXID)
table(test99$genus == test10$genus)
table(test99$family == test10$family)

hist(test99$LCA_SCORE)
