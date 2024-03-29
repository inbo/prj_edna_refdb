
library(tidyverse)
library(googlesheets4)
gs4_auth("pieter.verschelde@inbo.be")

###

inputfile <- "1NidznDq9EVHN4_wfrhv0W5gZ0Q-S8wEn5jvXbb9rf8k"


### GBIF TAXIDS
scinames <- read_sheet(inputfile, sheet = "Soortenlijst", range = "A:J") %>% 
  transmute(RN = 1 + (1:nrow(.)), Taxid, 
            scientificName = NameScientific, kingdom = NA, GBIF_Taxid,
            sciNameUpp = toupper((scientificName))) %>% 
  filter(!is.na(scientificName) & Taxid > 0 & is.na(GBIF_Taxid))

write_csv(scinames %>% select(scientificName, kingdom), 
          file = "output/gbif_input.csv", na = "")
#uploaden naar https://www.gbif.org/tools/species-lookup
#csv exporteren en hernoemen tot gbif_result_normalized
gbif_names <- read_csv("output/gbif_result_normalized.csv")
gbif_names <- gbif_names %>% 
  mutate(sciNameUpp = toupper(verbatimScientificName))

#result
linked <- scinames %>% select(-kingdom, scientificNameINBO = scientificName) %>% 
  left_join(gbif_names, by = "sciNameUpp") %>% 
  transmute(RN, Taxid, GBIF_Taxid = key, Rank = tolower(rank), Name = canonicalName,
            Kingdom = kingdom, Phylum = phylum, Class = class, 
            Order = order, Family = family, Genus = genus, Species = species)
write_excel_csv2(linked, file = "output/gbif_to_specieslist.csv")

############################
#### MANUEEL ##############
############################

indata <- readLines(n = 11)
Acipenser oxyrinchus oxyrinchus
Acipenser oxyrinchus desotoi
Carassius auratus auratus
Carassius auratus grandoculis
Carassius auratus ssp. 'Pingxiang'
Carassius auratus 'high back crucian carp'
Barbatula aff. barbatula lineage I
Barbatula aff. barbatula lineage II
Salmo caspius
Alburnus chalcoides aralensis
Paramisgurnus dabryanus ssp. DLY-2014

write_csv(data.frame(scientificName = indata, kingdom = NA), 
          file = "output/sp_20230829b.csv")
# TAAK: export naar gbif en resultaat opslaan als csv
gbif_names <- read_csv("output/gbif_result_20230829b.csv")
