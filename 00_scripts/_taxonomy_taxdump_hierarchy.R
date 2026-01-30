# SETUP

library(tidyverse)
curr_wd <- getwd()
#periodiek: downloaden van de volledige taxdump
source <- "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.zip"
source_tux <- "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz"
date <- Sys.Date()
filename <- paste0("new_taxdump_", date, ".zip")
filename_tux <- paste0("new_taxdump_", date, ".tar.gz")
destdir <- file.path("taxonomy", paste0("taxdump_", date))
if (!dir.exists("taxonomy")) 
  dir.create('taxonomy')
if (!dir.exists(destdir)) {
  dir.create(destdir) 
  download_new <- TRUE
} else {
  download_new <- FALSE  
}
  

# DOWNLOAD TAXONOMY


if (download_new) {
  #don't run if not necessary
  setwd(destdir)
  cat("download dir: ", getwd(), "\n")
  cat("downloading taxonomy.zip ...\n")
  download.file(source, destfile = filename)
  cat("downloading taxonomy.tar.gz ...\n")
  download.file(source_tux, destfile = filename_tux)
  cat("extracting ...\n")
  system2('c:/Program Files/7-Zip/7z.exe', 
          args = c("e", filename))
  
  # kopieer naar DOCKER (er vanuit gaande dat docker container obitools runt)
  cat("copying to docker ...\n")
  try({
    taxdir <- paste0("taxdump_", date)
    #system2("docker", args = c("exec", "obitools", "mkdir", taxdir))
    system2("docker", args = c("cp", file.path("../", taxdir), "obitools:/"))    
  })

  setwd(curr_wd);getwd()  
}

# NAMES

df_names <- 
  read_delim(file.path(destdir, "names.dmp"), delim = "|",
             trim_ws = TRUE, col_names = FALSE, col_select = c(1,2,4)) %>% 
  filter(X4 == "scientific name") %>% 
  select(taxid = X1, scientific_name = X2)

df_nodes <- 
  read_delim(file.path(destdir, "nodes.dmp"), delim = "|",
             trim_ws = TRUE, col_names = FALSE, col_select = c(1,2,3)) %>% 
  select(taxid = X1, rank = X3, , parent_taxid = X2) %>% 
  filter(rank != "no rank" ) %>% 
  left_join(df_names) 

#HIERARCHY

df_hierarchy_rd <-read_delim(file.path(destdir, "taxidlineage.dmp"), 
                             delim = "|", trim_ws = TRUE, col_names = FALSE)

## lees de hierarchie in en splits in kolommen (duurt een minuutje)
df_hierarchy <- df_hierarchy_rd %>% 
  select(-X3) %>% 
  separate(X2, into = paste0("L", sprintf("%02d", 1:36)))

## converteer naar lang formaat (duurt een 15-tal seconden)
df_hierarchy_long <- df_hierarchy %>%  
  pivot_longer(names_to = "order", cols = !X1, values_to = "taxidh") %>% 
  filter(!is.na(taxidh)) %>% 
  transmute(taxid = as.numeric(X1), taxidh = as.numeric(taxidh))

## link de hikerarchy met de wetenschappelijke namen (duurt een 15-tal seconden)
df_hierarchy_full <- df_hierarchy_long %>% 
  inner_join(df_nodes %>% select(taxid, rank, sci_name_taxid = scientific_name)) %>% 
  inner_join(df_nodes %>% select(taxidh = taxid, rankh = rank, scientific_name), 
            by = c("taxidh")) %>% 
  filter(rankh %in% c("kingdom", "phylum", "class", "order", 
                     "family", "genus", "species", "subspecies")) %>% 
  select(taxid, rank, sci_name_taxid, taxidh, rankh, scientific_name)

# SAVE IN TAXONDB
con <- DBI::dbConnect(RSQLite::SQLite(), 
                 file.path(destdir, "taxondb.sqlite")) 
DBI::dbWriteTable(con, "tblNames", df_nodes, overwrite = TRUE)
DBI::dbWriteTable(con, "tblHierarchy", df_hierarchy_full, overwrite = TRUE)

