
library(tidyverse)
library(googlesheets4)
gs4_auth("pieter.verschelde@inbo.be")

###

inputfile <- "1NidznDq9EVHN4_wfrhv0W5gZ0Q-S8wEn5jvXbb9rf8k"


### GBIF TAXIDS
scinames <- read_sheet(inputfile, sheet = "Soortenlijst", range = "D:H") %>% 
  transmute(RN = 1 + (1:nrow(.)), Taxid, 
            scientificName = NameScientific, kingdom = NA,
            sciNameUpp = toupper((scientificName))) %>% 
  filter(!is.na(scientificName))

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
  