source("scripts/update_refdb_00_init.R")

ecopcr_riaz_file <- "amplified_clean.fasta"
merged_riaz_file <- "amplified_clean_uniq.fasta"
input_fasta_file <- "input_obi3.fasta"

outpath <- file.path("output", db_name)

con <- DBI::dbConnect(RSQLite::SQLite(), 
                      file.path("taxonomy", 
                                taxdump_name, 
                                "taxondb.sqlite")) 

species_names <- tbl(con, "tblNames")

soortenlijst <- read_sheet(input_source, sheet = "Soortenlijst") %>% 
  filter(Priority %in% c(1,2), Taxid > 0) %>% 
  rename(taxid = Taxid, priority = Priority)

allowed_merges <- read_sheet(input_source, sheet = "Toegelaten_merges") 

input_data <- 
  parse_refdb_fasta(file.path("database", db_name, input_fasta_file), 
                    is_merged_file = FALSE) %>% 
  mutate(taxid = as.numeric(taxid))

ecopcr_data <- 
  parse_refdb_fasta(file.path("database", db_name, ecopcr_riaz_file), 
                    is_merged_file = FALSE)
colnames(ecopcr_data) <- tolower(colnames(ecopcr_data))

merged_data <- 
  parse_refdb_fasta(file.path("database", db_name, merged_riaz_file), 
                    is_merged_file = TRUE) %>% 
  mutate(merged_count = as.numeric(COUNT), 
         COUNT = NULL) 
colnames(merged_data) <- tolower(colnames(merged_data))

save(file = "ecopcr_base_data.Rdata", soortenlijst, allowed_merges, 
     input_data, ecopcr_data, merged_data)
#load("ecopcr_base_data.Rdata")

###

ecopcr_combined <- combine_ecpocr_with_merged(ecopcr_data, merged_data)

#behou enkel soorten in de soortenlijst
#enkel records waar merged_rank slechter dan species is 
#combineer identieke amplicons uit inputs 
ecopcr_combined_filtered_sp <- ecopcr_combined %>% 
  inner_join(soortenlijst %>% select(taxid, priority)) %>% 
  select(genlab_id, taxid, rank, priority, species_name,
         amplicon_hash, obi_rank, obi_taxid, genus_name, family_name,
         merged_overview, obi_count) %>% 
  filter(!(obi_rank %in% c("subspecies", "species"))) %>% 
  group_by(across(-genlab_id)) %>% 
  summarise(genlab_id = paste(genlab_id, collapse = ";")) %>% 
  arrange(obi_taxid, taxid, amplicon_hash)


table(ecopcr_combined$is_merged)
sum(table(ecopcr_combined$is_merged)) == nrow(ecopcr_data)
table(ecopcr_combined$obi_rank)

### teloor gegane inputsequenties
lost_inputs <- input_data %>% anti_join(ecopcr_combined, by = "genlab_id") %>% 
  left_join(species_names, by = "taxid", copy = TRUE) %>% 
  select(genlab_id, taxid, scientific_name, dna_sequence)
write_excel_csv2(lost_inputs, file.path(outpath, "lost_inputs.csv"))

### overzicht soorten in de referentiedb
ecopcr_combined %>% 
  group_by(family_name) %>% 
  summarise(verschillende_genera = n_distinct(genus), 
            verschillende_taxa = n_distinct(taxid)) %>% 
  arrange(desc(verschillende_taxa)) %>% 
  write_excel_csv2(file = file.path(outpath, "aantal_taxa_per_familie.csv"))

ecopcr_combined %>% 
  group_by(taxid) %>% 
  summarise(species = paste(unique(species_name), sep = ";"),
            genus = paste(unique(genus_name), sep = ";"),
            family = paste(unique(family_name), sep = ";")) %>% 
  arrange(family, genus, species) %>% 
  write_excel_csv2(file = file.path(outpath, "aanwezige_taxa.csv"))




### evaluatie gemergde sequenties

#selecteer enkel de gemergde sequenties (obi_count > 1)
#die niet minstens op soortniveau gemerged zijn
#en die niet op een toegelaten merge op genus samengevoegd zijn
ecopcr_merges <- ecopcr_combined %>% 
  filter(obi_count > 1, !(obi_rank %in% c("subspecies", "species"))) %>% 
  filter(!(obi_rank %in% c("subgenus", "genus", "species group") &
           genus %in% (allowed_merges %>% 
                         filter(RANK == "genus") %>% 
                         pull(TAXID)))) %>% 
  left_join(soortenlijst %>% select(taxid, priority)) %>% 
  select(genlab_id, taxid, rank, priority, species_name,
         amplicon_hash, obi_rank, obi_taxid, genus_name, family_name,
         merged_overview, obi_count) %>% 
  arrange(obi_taxid, taxid, amplicon_hash)

#hieruit alle merges halen waarvoor er meer dan 1 soort in de soortenlijst zit
#behou enkel een conflict indien met meer dan 1 relevante soort
relevant_amplicons <- ecopcr_merges %>% 
  group_by(amplicon_hash) %>% 
  summarize(relevant_species = sum(!is.na(priority))) %>% 
  filter(relevant_species > 1) 


#vector met alle amplicons die een merge hebben met minstens 1 soortenlijst soort
#alle records met deze probleemamplicons
amplicons_with_merges <- ecopcr_merges %>% 
  filter(amplicon_hash %in% (relevant_amplicons %>% pull(amplicon_hash))) %>% 
  arrange(amplicon_hash)

#haal  de probleemsoorten op
taxa_in_amplicons_with_merges <- amplicons_with_merges %>% 
  select(taxid, species_name, priority) %>% 
  filter(!is.na(priority)) %>% 
  distinct()

#haal alle amplicons op die 1 van deze probleemsoorten heeft
#filter de records weg zonder merged_overview, 
#alsook die die op soortniveau gemerged zijn in een ander amplicon
#want die zitten per definitie ook in een record met een aanwezige merged_overview
amplicons_multihit <- ecopcr_combined_filtered_sp %>% 
  filter(taxid %in% taxa_in_amplicons_with_merges$taxid)

amplicons_multihit_eval <- amplicons_multihit %>% 
  group_by(obi_taxid) %>% 
  do(check_multihits(., specieslist = soortenlijst))

multihit_list_input <- amplicons_multihit_eval %>% 
  rowwise() %>% 
  do(multihit_list_update(.)) %>% 
  filter(pref_taxid != -99) %>% 
  distinct() %>% 
  left_join(input_data %>% select(genlab_id, taxid)) %>% 
  group_by(across(-genlab_id)) %>% 
  summarise(genlab_id = paste(genlab_id, collapse = ";"))

multihit_tabel <- 
  link_species_multihit_check(multihit_list_input, db_name, soortenlijst) %>% 
  arrange(pref_taxid, taxid)
  
write_excel_csv2(multihit_tabel, 
                 file = file.path("output", db_name, "new_multihits.csv"))

