library(tidyverse)
library(googlesheets4)
gs4_auth("pieter.verschelde@inbo.be")

inputfile     <- "1NidznDq9EVHN4_wfrhv0W5gZ0Q-S8wEn5jvXbb9rf8k"

soortenlijst <- read_sheet(inputfile, "Soortenlijst") %>% 
  mutate(RN = 1:n())

multihitlist <- read_sheet(inputfile, "Multihitlist_riaz") %>% 
  select(TAXID, PREF_TAXID) %>% 
  mutate(RN = 1:n()) 

multihitcheck <- multihitlist %>% 
  left_join(soortenlijst %>% 
              select(Taxid, SCI_NAME = NameScientific, 
                     DUTCH_NAME=NameDutch, 
                     ENG_NAME = NameEnglish),
            by = c("TAXID" = "Taxid"))
  left_join(soortenlijst %>% 
              select(Taxid, SCI_NAME_PREF=NameScientific, 
                     DUTCH_NAME_PREF=NameDutch,
                     )
              )




