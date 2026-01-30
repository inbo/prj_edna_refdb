#let op voor _mod of modified
source("scripts/refdb_teleo_00_initialisation.R")

# 

path <- "G:\\.shortcut-targets-by-id\\0B1XJuciaZSENZG55ZnlDQ0FvT0E\\PRJ_eDNA\\PRJ_eDNA_Refdb_2023\\input_seqs\\import"

files <- list.files(path, pattern = "*.fasta")
df_inputs_all <- NULL
for (file in files) {
  cat("\n\nINLEZEN VAN ", file, "\n--------------------------------\n")
  parsed <- parse_refdb_fasta(file.path(path, file))
  df_inputs_all <- df_inputs_all %>% 
    bind_rows(parsed)
}
df_inputs_all <- df_inputs_all %>% 
  mutate(len = nchar(dna_sequence),
         id8 = substring(genlab_id,1,8),
         id9 = substring(genlab_id, 1, 9),
         id10 = substring(genlab_id, 1, 10))


dupids <- df_inputs_all$genlab_id[which(duplicated(df_inputs_all$id9))]
df_inputs_all %>% filter(genlab_id %in% dupids ) %>% view()

df_inputs_all %>% 
  group_by(id8) %>% 
  summarise(n_taxid = n_distinct(taxid), 
            taxids = paste(taxid, collapse = "|"), 
            sources = paste(unique(source), collapse = ","),
            genlab_ids = paste(unique(genlab_id), collapse = ",")) %>% 
  filter(n_taxid > 1) %>% 
  arrange(id8) %>% 
  view()



fulldups <- unique(df_inputs_all$genlab_id[which(duplicated(df_inputs_all[c('id10', 'taxid', 'DNA_HASH')]))])
partdups <- unique(df_inputs_all$genlab_id[which(duplicated(df_inputs_all[c('genlab_id', 'taxid')]))])


df_inputs_all %>% filter(genlab_id %in% fulldups) %>% arrange(genlab_id, source) %>% view()


dups <- df_inputs_all %>% 
  filter(genlab_id %in% fulldups, 
         source %in% c("2019_08_30.2_Teleo2019", "2020_03_11.2 Riaz_March2020"),
         !(genlab_id %in% c("NC_000860", "EU880339"))) %>% 
  arrange(genlab_id, source) %>% view()


shareds <- df_inputs_all %>% filter(source == "2020_03_12 Shared refseqs_Riaz_Teleo")

origs <- parse_refdb_fasta(file.path(path, "2019_08_30.2a_Teleo2019.fasta"))

repair <- shareds %>% left_join(origs %>% select(idcor= genlab_id, taxid, DNA_HASH))






# s2017 <- parse_refdb_fasta(file.path(path1, "2017_12_31 Riaz_2017.fasta"))
# s2020 <- parse_refdb_fasta(file.path(path1, "2020_03_11.2 Riaz_March2020.fasta"))
# table(substring(s2017$genlab_id, 1, 8) %in% substring(s2020$genlab_id, 1, 8))
# 
# s2017 %>%
#   transmute(genlab_id, taxid, DNA_HASH, id = substring(genlab_id, 1, 8)) %>%
#   inner_join(s2020 %>%
#                transmute(genlab_id2=genlab_id, taxid2=taxid, DNA_HASH2=DNA_HASH, id = substring(genlab_id, 1, 8))) %>%
#   transmute(genlab_id, genlab_id2, taxid, taxid2, taxid==taxid2, DNA_HASH, DNA_HASH2, DNA_HASH==DNA_HASH2) %>%
#   view()

teleoseqs <- parse_refdb_fasta(file.path(path2, "2019_08_30.2_Teleo2019.fasta"))
riazseqs <- parse_refdb_fasta(file.path(path2, "2020_03_11.2 Riaz_March2020.fasta"))


teleoseqs %>%
  transmute(genlab_id, taxid, DNA_HASH, id = substring(genlab_id, 1, 10)) %>%
  inner_join(riazseqs %>%
               transmute(genlab_id2=genlab_id, taxid2=taxid, DNA_HASH2=DNA_HASH, id = substring(genlab_id, 1, 10))) %>%
  transmute(genlab_id, genlab_id2, taxid, taxid2, taxid==taxid2, DNA_HASH, DNA_HASH2, DNA_HASH==DNA_HASH2) %>%
  view()







# 
# mgseqs <- parse_refdb_fasta(file.path(path, "2019_03_11 Matthias Geiger 1087 NCBI with taxids.fasta"))
# teleoseqs <- parse_refdb_fasta(file.path(path2, "2019_08_30.2_Teleo2019.fasta"))
# 
# dups <- mgseqs %>% filter(substring(genlab_id, 1, 8) %in% substring(teleoseqs$genlab_id, 1, 8))
# 
# mgseqs %>% 
#   transmute(genlab_id, taxid, DNA_HASH, id = substring(genlab_id, 1, 8)) %>% 
#   inner_join(teleoseqs %>% 
#                transmute(genlab_id2=genlab_id, taxid2=taxid, DNA_HASH2=DNA_HASH, id = substring(genlab_id, 1, 8))) %>% 
#   transmute(genlab_id, genlab_id2, taxid, taxid2, taxid==taxid2, DNA_HASH, DNA_HASH2, DNA_HASH==DNA_HASH2) %>% 
#   view()
# 
# 
# ampseqs <- parse_refdb_fasta(file.path(path, "20220104 Amphibian new refseqs.fasta")) %>% 
#   mutate(source = "20220104 Amph")
# fishseqs <- parse_refdb_fasta(file.path(path, "20220104 Fish new refseqs.fasta")) %>% 
#   mutate(source = "20220104 Fish")
# 
# gobiseqs <- parse_refdb_fasta(file.path(path, "Copy of 2022_02_10 5 modified Gobiosoma bosc.fasta")) %>% 
#   mutate(source = "Gobi")
# salvseqs <- parse_refdb_fasta(file.path(path, "Copy of 2022_02_11 Salvelinus fontinalis.fasta")) %>% 
#   mutate(source = "Salvelinus")
# ponderseqs <- parse_refdb_fasta(file.path(path, "Copy of 2022_05_17 Ponderful update.fasta")) %>% 
#   mutate(source = "Ponderful")
# newriazseqs <- parse_refdb_fasta(file.path(path, "Copy of 2022_10_18 Riaz update.fasta")) %>% 
#   mutate(source = "NewRiaz")
# 
# smallseqs <- bind_rows(ampseqs, fishseqs, gobiseqs, salvseqs, ponderseqs, newriazseqs) %>% 
#   arrange(genlab_id)
# 
# mgseqs %>% 
#   bind_rows(teleoseqs) %>% 
#   transmute(genlab_id, taxid, DNA_HASH, id = substring(genlab_id, 1, 8)) %>% 
#   inner_join(smallseqs %>% 
#                transmute(genlab_id2=genlab_id, taxid2=taxid, DNA_HASH2=DNA_HASH, id = substring(genlab_id, 1, 8))) %>% 
#   transmute(genlab_id, genlab_id2, taxid, taxid2, taxid==taxid2, DNA_HASH, DNA_HASH2, DNA_HASH==DNA_HASH2) %>% 
#   view()


inputseqs_orig <- parse_refdb_fasta(file.path(path1, "2020_03_11 Riaz_March2020.fasta"))

inputseqs <- inputseqs_orig %>%
  mutate(seq_len = nchar(dna_sequence),
         new_id = NA)

for (i in seq_len(nrow(inputseqs))) {
  name <- nameT <- len <- NA
  genlab_id <- inputseqs$genlab_id[i]
  len <- inputseqs$seq_len[i]
  dots <- gregexpr("\\.", genlab_id)
  first_dot <- dots[[1]][1]
  second_dot <- dots[[1]][2]
  if (!is.na(second_dot)) {
    name <- paste(substring(genlab_id, 1, second_dot - 1))
  } else {
    name <- genlab_id
  }
  name <- gsub("\\(modified\\)", ".mod", name)
  name <- gsub("modRiaz", ".mod", name)
  if (len < 200) {
    nameT <- paste(name, "R", sep = ".")
  } else {
    nameT <- name
  }
  inputseqs$new_id[i] <- substring(nameT, 1, 20)
}

view(inputseqs %>% select(genlab_id, new_id, seq_len))

tmp1 <- regexpr("MOD", inputseqs$genlab_id)>0
tmp2 <- regexpr("mod", inputseqs$genlab_id)>0
tmp3 <- regexpr("Mod", inputseqs$genlab_id)>0
table(tmp1)
table(tmp2)
table(tmp3)
inputseqs[tmp2,]

### Find duplicate seqs
table(duplicated(inputseqs %>% select(DNA_HASH)))
table(duplicated(inputseqs %>% select(DNA_HASH, taxid)))
table(duplicated(inputseqs %>% select(genlab_id)))
table(duplicated(inputseqs %>% select(new_id)))

which(duplicated(inputseqs %>% select(genlab_id)))
which(duplicated(inputseqs %>% select(new_id)))

inputseqs %>% 
  filter(genlab_id %in% inputseqs$genlab_id[which(duplicated(inputseqs %>% select(new_id)))]) %>% 
  arrange(genlab_id) %>% 
  view()

teleoseqs <- parse_refdb_fasta(file.path(path2, "2019_08_30.2_Teleo2019.fasta"))


inputseqs %>%
  transmute(genlab_id, taxid, DNA_HASH, id = substring(genlab_id, 1, 10)) %>%
  inner_join(teleoseqs %>%
               transmute(genlab_id2=genlab_id, taxid2=taxid, DNA_HASH2=DNA_HASH, id = substring(genlab_id, 1, 10))) %>%
  transmute(genlab_id, genlab_id2, taxid, taxid2, taxid==taxid2, DNA_HASH, DNA_HASH2, DNA_HASH==DNA_HASH2) %>%
  view()
inputseqs_mod <- inputseqs

# 
# inputseqs %>%  filter(new_id == "AB777216.1.T") #HOE MET DIT GEVAL OMGAAN? VERSCHILLENDE LOCATIE, 1 BP VERSCHIL
# inputseqs <- inputseqs %>%
#   mutate(new_id = ifelse(genlab_id == "AB777216.1.9418.10356", "AB777216.1.B.T", genlab_id),
#          new_id = ifelse(genlab_id == "AB777216.1.6671.7609",  "AB777216.1.A.T", genlab_id))
# 
# ### Foutieve sequenties corrigeren
# 
# correcties <- read_csv2("database/refdb_teleo_2023-08-10/te_controleren_seqs.csv") %>% select(1:4)
# 
# inputseqs_mod <- inputseqs %>% filter(!genlab_id %in% c("JN627425", "EU880339"))
# inputseqs_mod <- inputseqs_mod %>% left_join(correcties %>% select(new_id, gevonden_taxon))
# inputseqs_mod <- inputseqs_mod %>%
#   mutate(new_id = ifelse(is.na(inputseqs_mod$gevonden_taxon), new_id, str_replace(new_id,'.T', ".mod.T" )),
#          taxid =  ifelse(is.na(inputseqs_mod$gevonden_taxon), taxid, gevonden_taxon))
create_input_fasta("newriaz.fasta",
                   data = inputseqs_mod,
                   seq = "dna_sequence",
                   taxid = "taxid",
                   entry_id = "new_id"
                   )


