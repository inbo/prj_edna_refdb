#ZORG DAT DOCKER DESKTOP RUNT en OBITOOLS geladen is
#run de container obitools_pv 

#################
## INIT
#################

source(file.path("scripts", "refdb_teleo_00_initialisation.R"))

docker_container <- "obitools_pv"
fasta_name <- "input.fasta"
db_namepath <- paste0(db_name, "/", db_name) 

e <- 5
l <- 20
L <- 100
p1 <- select_primer_tags("TELEO")["oligo1"]
p2 <- select_primer_tags("TELEO")['oligo2']

#######################
### DOCKER PROCESSING
#######################

## Maak structuur aan en kopieer input.fasta
system2("docker", args = c("exec", docker_container, "mkdir", db_name))
system2("docker", args = c("cp", file.path("database", db_name, fasta_name), 
                           paste0(docker_container,":/",db_name)))

#KOPIEREN TAXONOMY indien nog niet gebeurd
system2("docker", args = c("exec", docker_container, "rm", "-r", "taxdump"))
system2("docker", args = c("exec", docker_container, "mkdir", "taxdump"))
system2("docker", args = c("cp", file.path(taxdump_name), 
                           paste0(docker_container,":/","taxdump/",  "taxdump.tar.gz")))
system2("docker", args = c("exec", docker_container, "tar", "-xzvf", "taxdump/taxdump.tar.gz", "-C", "taxdump"))

## converteer input+taxonomie naar voor obitools bruikbaar formaat

obiconvertstring <- paste("obiconvert", "--fasta", paste0("--ecopcrdb-output=",db_namepath),
                          "-t", "taxdump", 
                          paste0(db_name, "/", fasta_name))
system2('docker', args = c('exec', docker_container, 'sh', '-c', paste0('"', obiconvertstring, '"')))

## voer de in silico PCR uit
#als ecoPCR niet gekend is ./home/ecopcr/src/ecoPCR
ecopcrstring <- paste("/home/ecopcr/src/ecoPCR", "-d", db_namepath, 
                      "-e", e, "-l", l, "-L", L, p1, p2, 
                      ">", paste0(db_name, "/", "amplified2.ecopcr"))
system2('docker', args = c('exec', docker_container, 'sh', '-c', paste0('"', ecopcrstring, '"')))

## clean:(hou enkel die over die een parent rank hebben die soort en familie zijn)

grepstring <- paste("obigrep", "-d", db_namepath, 
                    "--require-rank=species",
                    "--require-rank=family",
                    paste0(db_name, "/", "amplified2.ecopcr"),
                    ">", 
                    paste0(db_name, "/", "amplified_clean2.fasta"))
system2('docker', args = c('exec', docker_container, 'sh', '-c', paste0('"', grepstring, '"')))

## Keep unique sequences (merge identical amplicon sequences)

uniqstring <- paste("obiuniq", "-d", db_namepath, 
                    paste0(db_name, "/", "amplified_clean2.fasta"),
                    ">", 
                    paste0(db_name, "/", "amplified_clean_uniq2.fasta"))
system2('docker', args = c('exec', docker_container, 'sh', '-c', paste0('"', uniqstring, '"')))

## Copy back to windows
system2("docker", args = c("cp", paste0(docker_container, ':', db_name),file.path("database")))


#############################
### POSTPROCESSING
#############################

ecopcr_file <- "amplified_clean2.fasta"
merged_file <- "amplified_clean_uniq2.fasta"
input_fasta_file <- "input.fasta"

df_soortenlijst_all <- read_sheet(metadata_gdrive_key, sheet = "Soortenlijst") %>% 
  filter(Taxid > 0) %>% 
  rename(taxid = Taxid, priority = Priority)
df_soortenlijst <- df_soortenlijst_all %>% 
  filter(priority %in% c(1,2,3,4))
df_ok_merges <- read_sheet(metadata_gdrive_key, sheet = "Toegelaten_merges_Teleo") 

## ECOPCR INHOUD

input_data2 <- 
  parse_refdb_fasta(file.path("database", db_name, "kept_input2.fasta"), 
                    is_merged_file = FALSE)

ecopcr_data2 <- 
  parse_refdb_fasta(file.path("database", db_name, ecopcr_file), 
                    is_merged_file = FALSE) %>% 
  rename(taxid = taxid,
         amplicon_hash = AMPLICON_HASH)

merged_data2 <- 
  parse_refdb_fasta(file.path("database", db_name, merged_file), 
                    is_merged_file = TRUE) %>% 
  rename(merged_count = merged_count, 
         amplicon_hash = AMPLICON_HASH,
         merged_taxid = merged_taxid,
         taxid = taxid
  )
ecopcr_combined2 <- combine_ecpocr_with_merged(ecopcr_data2, merged_data2) 

table(ecopcr_combined2$is_merged)
table(ecopcr_combined2$obi_rank)
table(ecopcr_combined2$obi_rank, ecopcr_combined2$is_merged)

ecopcr_filtered_sp2 <- ecopcr_combined2 %>% 
  inner_join(df_soortenlijst %>% select(taxid, priority)) %>% 
  select(genbank_id, taxid, rank, priority, species_name,
         amplicon_hash, obi_rank, obi_taxid, genus_name, family_name,
         merged_overview, obi_count) %>% 
  filter(!(obi_rank %in% c("subspecies", "species"))) %>% 
  group_by(across(-genbank_id)) %>% 
  summarise(genbank_id = paste(genbank_id, collapse = ";")) %>% 
  arrange(obi_taxid, taxid, amplicon_hash) %>% 
  dplyr::filter(!(obi_rank == "genus" & 
                    obi_taxid %in% (df_ok_merges %>% filter(RANK == "genus") %>% pull(TAXID))), 
                !(obi_rank == "family" & 
                    obi_taxid %in% (df_ok_merges %>% filter(RANK == "family") %>% pull(TAXID))  ))

df_spm2 <- ecopcr_filtered_sp2 %>%
  arrange(obi_taxid) %>%
  group_by(obi_taxid) %>%
  do ({
    get_taxa_from_merged(., df_soortenlijst)
  }
  )

df_beoordeeld2 <- df_spm2 %>% group_by(obi_taxid) %>%
  do ({
    judge_species(.)
  })

beoordeling2 <- df_beoordeeld2 %>% group_by(taxid, pref_taxid) %>%
  summarise(oordeel = paste(unique(oordeel), collapse = "|"))

df_soortentabel2 <- df_soortenlijst %>% 
  group_by(Group, NameScientific, NameEnglish, NameDutch, taxid, priority) %>% 
  summarise(n_input_orig = sum(input_data_before_multihit$taxid == .data$taxid[1]),
            #n_input_seqs = sum(as.numeric(input_data$TAXID) == .data$taxid[1]),
            n_input_seqs = sum(as.numeric(input_data2$TAXID == .data$taxid[1]) > 0, na.rm = TRUE),
            n_amplicons = sum(ecopcr_combined2$taxid == .data$taxid[1]),
            n_determined = sum(ecopcr_combined2$obi_taxid == .data$taxid[1]))

ecopcr_soorten2 <- ecopcr_combined2 %>% 
  right_join(df_soortenlijst %>% select(taxid)) %>% 
  group_by(taxid) %>% 
  summarise(merged = sum(is_merged)>0, 
            obi_taxid = paste(unique(obi_taxid), collapse = '|'),
            obi_rank = paste(unique(obi_rank), collapse = "|")) %>% 
  mutate(merged = merged & taxid != obi_taxid)

df_multihits2 <-
  df_beoordeeld2 %>% group_by(taxid) %>%
  do({
    taxid <- .$taxid[1]
    obitaxids <- (unique(.$obi_taxid))
    taxids <- df_beoordeeld2$taxid[df_beoordeeld2$obi_taxid %in% obitaxids]
    data.frame(multihit_taxa = paste("|", sort(unique(taxids)), collapse = '', sep = ''))
  })

df_soortenevaluatie2 <- df_soortentabel2 %>% 
  left_join(ecopcr_soorten2) %>% 
  left_join(beoordeling2) %>% 
  left_join(df_multihits2) %>% 
  mutate(pref_taxid = ifelse(merged == FALSE, taxid, pref_taxid)) %>% #klopt niet, merged is fout
  left_join(df_soortenlijst %>% select(taxid, NameEnglish, NameScientific), by = c("pref_taxid" = "taxid"))

write_excel_csv2(df_soortenevaluatie2, file = paste0(output_path,"/", "soortenevaluatie_obi2_ruw.csv"))


