library(tidyverse)
library(googlesheets4)

root_gdrive <- file.path("G:", 
                         ".shortcut-targets-by-id",
                         "0B1XJuciaZSENZG55ZnlDQ0FvT0E", 
                         "PRJ_eDNA",
                         "PRJ_eDNA_Refdb_2023" )
fasta_inputs <- file.path(root_gdrive, 'input_seqs', 'import')

teleo_path <- "C:\\_GIT_PROJECTS\\prj_edna_refdb\\database\\refdb_teleo_2023-09-01"
riaz_path  <- "C:\\_GIT_PROJECTS\\prj_edna_refdb\\database\\refdb_riaz_2023-09-01"
inputchecks_name <- "inputchecks_2023_09_01.Rdata"
source("scripts/_functions_fasta.R")

files <- sort(list.files(fasta_inputs, pattern = ".fasta"))
 all_inputs <- NULL
for (file in files) {
  cat("\n\nINLEZEN VAN ", file, "\n--------------------------------\n")
  parsed <- parse_refdb_fasta(file.path(fasta_inputs, file))
  all_inputs <- all_inputs %>% 
    bind_rows(parsed)
}
 
find_fragments <- function(x) {
  
}

amplified_teleo <- parse_refdb_fasta(file.path(teleo_path, 'amplified.fasta'))
amplified_riaz  <- parse_refdb_fasta(file.path(riaz_path,  'amplified.fasta'))
amplified_clean_teleo <- parse_refdb_fasta(file.path(teleo_path, 'amplified_clean.fasta'))
amplified_clean_riaz  <- parse_refdb_fasta(file.path(riaz_path,  'amplified_clean.fasta'))
#amplified_clean_uniq_teleo <- parse_refdb_fasta(file.path(teleo_path, 'amplified_clean_uniq.fasta'))
#amplified_clean_uniq_riaz  <- parse_refdb_fasta(file.path(riaz_path,  'amplified_clean_uniq.fasta'))
refdb_teleo <- parse_refdb_fasta(file.path(teleo_path, 'final_db_0.99.fasta'))
refdb_riaz  <- parse_refdb_fasta(file.path(riaz_path,  'final_db_0.99.fasta'))

save(all_inputs, amplified_teleo, amplified_riaz,
     amplified_clean_teleo, amplified_clean_riaz, 
#    amplified_clean_uniq_teleo, amplified_clean_uniq_riaz, 
     refdb_teleo, refdb_riaz, file = inputchecks_name)
#load(inputchecks_name)

whidup <- which(duplicated(all_inputs$genbank_id))
all_inputs %>% 
  filter(genbank_id %in% all_inputs$genbank_id[whidup]) %>% 
  mutate(seq_len = nchar(dna_sequence)) %>% 
  view()

overview <- all_inputs %>% 
  left_join(amplified_teleo %>% transmute(genbank_id, taxid = TAXID, teleo_amplified = TRUE)) %>% 
  left_join(amplified_riaz %>% transmute(genbank_id, taxid = TAXID, riaz_amplified = TRUE)) %>% 
  mutate(teleo_amplified = ifelse(is.na(teleo_amplified), FALSE, TRUE),
         riaz_amplified = ifelse(is.na(riaz_amplified), FALSE, TRUE)) %>% 
  mutate(amplified_status = ifelse(teleo_amplified & riaz_amplified, 
                                   'both',
                                   ifelse(teleo_amplified, 
                                          "teleo",
                                          ifelse (riaz_amplified,
                                                  'riaz',
                                                  'none'))))

not_amplified <- overview %>% filter(amplified_status == "none")
write_excel_csv2(not_amplified, "output/niet_geamplificeerd_voor_teleo_en_riaz.csv")
table(not_amplified$source)
  


