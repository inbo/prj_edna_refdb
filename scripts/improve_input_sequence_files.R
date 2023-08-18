#let op voor _mod of modified

path <- "G:\\.shortcut-targets-by-id\\0B1XJuciaZSENZG55ZnlDQ0FvT0E\\PRJ_eDNA\\PRJ_eDNA_Refdb_2023\\input_seqs\\Teleo\\import"

inputseqs_orig <- parse_refdb_fasta(file.path(path, "Teleo_RefDB_aug2019.fasta"))

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
    nameT <- paste(name, "T", sep = ".")
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
which(duplicated(inputseqs %>% select(new_id)))
inputseqs %>%  filter(new_id == "AB777216.1.T") #HOE MET DIT GEVAL OMGAAN? VERSCHILLENDE LOCATIE, 1 BP VERSCHIL 


inputseqs %>% filter(DNA_HASH == "a85ce5775b2a8dadf8f3f13463ecc33b5c3de584")



