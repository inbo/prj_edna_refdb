source("scripts/refdb_00_initialisation.R")
library(tidyverse)

## SPECIES

df_species_all <- read_sheet(metadata_gdrive_key, "Soortenlijst") %>% 
  mutate(ROWNR = 1:n())

df_species <- df_species_all %>% 
  filter(Taxid > 0, !is.na(Priority))

df_priority <- df_species %>% 
  select(Taxid, Priority)

df_passlist <- read_sheet(metadata_gdrive_key, "Pass_list_Riaz")

df_seq_errors <- read_sheet(metadata_gdrive_key, "Sequentie_fouten")

df_multihit <-  read_sheet(metadata_gdrive_key, "Multihitlist_riaz")

df_allowed_merges <- read_sheet(metadata_gdrive_key, "Toegelaten_merges")


### INPUTS

#get the source files
files <- sort(list.files(fasta_inputs, pattern = ".fasta"))

#read the source files
df_inputs_all <- NULL
for (file in files) {
  cat("\n\nINLEZEN VAN ", file, "\n--------------------------------\n")
  parsed <- parse_refdb_fasta(file.path(fasta_inputs, file))
  df_inputs_all <- df_inputs_all %>% 
    bind_rows(parsed)
}
whidup <- which(duplicated(df_inputs_all %>%  select(genlab_id, taxid, dna_sequence)))
df_inputs_all <- df_inputs_all[-whidup, ]

#find duplicate ids with different species
dup_ids <- df_inputs_all %>% group_by(genlab_id) %>% 
  summarise(aantal_seqs = n(), 
            sources = paste(source, collapse = ","),
            aantal_taxa = n_distinct(taxid), 
            taxa = paste(unique(taxid), collapse = ",")) %>% 
  filter(aantal_seqs > 1 | aantal_taxa > 1)

#Er zijn nog 7 duplicaten met andere sequenties, waarvan 1 met een andere soort (zelfde genus) nl 14-015707

## CURATE INPUTS

df_inputs <- df_inputs_all %>% 
  filter(!(genlab_id %in% (df_seq_errors %>% pull(ENTRY_ID)))) %>% 
  filter(!taxid %in% (df_seq_errors %>% 
                        filter(TYPE == "taxid") %>% pull(TAXID)))

### remove alle sequenties die in multihit taxa staan en waar taxid != pref_taxid
df_inputs <- df_inputs %>% 
  filter(!taxid %in% (df_multihit %>% 
                        filter(TAXID != PREF_TAXID) %>% 
                        pull(TAXID)))

### add passlist sequences
df_inputs <- df_inputs %>% bind_rows(df_passlist) 

### remove full duplicates
df_inputs <- df_inputs %>% distinct(genlab_id, taxid, AMPLICON_HASH, 
                                    .keep_all = TRUE)

## create fasta input for ecopcr (obitools2)

#create_input_fasta(file = file.path("database", db_name, fasta_name), data = df_inputs)

#voor obitools3
create_input_fasta(file = file.path("database", db_name, fasta_name), 
                   lowercase = TRUE, 
                   data = df_inputs)

################################################

make_shellscript_refdb(script = "generate_refdb.sh",
                       db_location = file.path("database",db_name), 
                       db_name = db_name,
                       input_file = "input.fasta",
                       taxonomy_location = file.path("taxonomy"),
                       taxonomy_file = "taxdump.tar.gz", 
                       environment_call = "/app/obi3-env/bin/activate")


#docker files
image_name <- "obitools3pv"
container_name <- "obitools3testcontainer"
obi_script <- "generate_refdb.sh"
refdb_location <- "database/refdb_2023-06-01/"
input_fasta <- paste0(refdb_location, "input.fasta")
taxdump_name <- "taxonomy/new_taxdump_2022-05-02.tar.gz"
docker_path <- 'obitools3testcontainer:/app/'
docker_script <- paste0('/app/', obi_script)
docker_taxonomy <- paste0(docker_path, "taxdump.tar.gz")
docker_input_fasta <- paste0(docker_path, "input.fasta")
pscommand_build <- paste("docker", "build", "-t", image_name, ".")
pscommand_run <- paste("docker run --rm -it --name", container_name, image_name)
pscommand_run <- paste("docker run -it --name", container_name, image_name)
#prep commands

system2("powershell", args = pscommand_build)
system2("powershell", args = pscommand_run) #mogelijks werkt dit niet --> via docker desktop
system2('docker',  c('cp', obi_script, paste0(container_name, ':', docker_script)))
system2('docker',  c('cp', taxdump_name, docker_taxonomy))
system2('docker',  c('cp', input_fasta, docker_input_fasta))

#execution
system2('docker', c('exec', container_name, 'bash', '-c', docker_script))

#return results
system2('docker', c("cp", paste0(docker_path, "logfile.txt"), refdb_location))
system2('docker', c("cp", paste0(docker_path, obi_script), refdb_location))
system2('docker', c("exec", container_name, "chmod", "777", obi_script))
system2('docker', c("cp", paste0(docker_path, "kept_input.fasta"), refdb_location))
system2('docker', c("cp", paste0(docker_path, "amplified_clean.fasta"), refdb_location))
system2('docker', c("cp", paste0(docker_path, "amplified_clean_uniq.fasta"), refdb_location))
system2('docker', c("cp", paste0(docker_path, "final_db_0.99.fasta"), refdb_location))
system2('docker', c("cp", paste0(docker_path, "refdb.obidms"), refdb_location))        