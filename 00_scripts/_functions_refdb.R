#' Maak ecopcr commando's
#'
#' @param script 
#' @param db_name 
#' @param db_location 
#' @param input_file 
#' @param taxonomy_location 
#' @param taxonomy_file 
#' @param max_errors 
#' @param min_length 
#' @param max_length 
#' @param primer 
#' @param ecotag_min_similarity 
#'
#' @return
#' @export
#'
#' @examples
make_shellscript_refdb <- 
  function(script = "", 
          db_name,
          db_location,
          input_file,
          taxonomy_location,
          taxonomy_file,
          max_errors = 4, 
          min_length = 50, 
          max_length = 160,
          primer = "RIAZ",
          ecotag_min_similarity = 0.99,
          environment_call =  "/app/obi3-env/bin/activate") {
  primers <- select_primer_tags(which = "RIAZ")
  p1 <- primers[1]
  p2 <- primers[2]
  wdwin <- file.path(getwd(), db_location)
  wdtux <- gsub(substring(wdwin, 1,2), 
                paste0("/mnt/", tolower(substring(wdwin, 1, 1))),
                wdwin)
  
  tdtux <- file.path(getwd(), taxonomy_location)
  tdtux <- gsub(substring(tdtux,1,2), 
                paste0("/mnt/", tolower(substring(tdtux, 1, 1))),
                tdtux)
  tdtux <- file.path(tdtux, taxonomy_file)
  commands <- 
    paste(sep = "\n",
          paste0("#!/bin/bash"),
          paste0("cd /app"),
          paste0(paste0("source ","/app/obi3-env/bin/activate")),
          paste0("echo 'aantal inputs: ' > logfile.txt"),
          paste0("grep -E  -i '>' ", 
                 input_file, 
                 " | wc -l >> logfile.txt"), #telt de lijnen die gegrept zijn
          paste0("obi import ", input_file, " refdb/input"),
          paste0("obi export", " refdb/input", " -o ", "kept_input.fasta"),
          paste0("echo 'aantal behouden inputs: ' >> logfile.txt"),
          paste0("grep -E  -i '>' ",
                 " kept_input.fasta", 
                 " | wc -l >> logfile.txt"),
          paste0("obi import --taxdump ", taxonomy_file,  " refdb/taxonomy/dump"),
          paste0("obi ecopcr -e ", max_errors, 
                 " -l ", min_length, 
                 " -L ", max_length, 
                 " -F ", p1, 
                 " -R ", p2, 
                 " --taxonomy ", "refdb/taxonomy/dump",
                 " refdb/input ", "refdb/ecopcr"),
          
          paste0("obi grep", 
                 " --require-rank=species",
                 " --require-rank=genus",
                 " --require-rank=family",
                 " --taxonomy refdb/taxonomy/dump", 
                 " refdb/ecopcr", 
                 " refdb/ecopcr_clean"),
          paste0("obi export", 
                 " refdb/ecopcr_clean",
                 " -o ", "amplified_clean.fasta"),
          paste0("echo 'aantal geamplificeerde inputs: ' >> logfile.txt"),
          paste0("grep -E  -i '>' ",
                 " amplified_clean.fasta", 
                 " | wc -l >> logfile.txt"), #telt de lijnen die gegrept zijn
          
          paste0("obi uniq", 
                 " --taxonomy refdb/taxonomy/dump", 
                 " refdb/ecopcr_clean",
                 " refdb/ecopcr_uniq"),
          paste0("obi grep",
                 " --require-rank=family",
                 " --taxonomy refdb/taxonomy/dump", 
                 " refdb/ecopcr_uniq ", 
                 " refdb/ecopcr_uniq_clean"),
          paste0("obi export", 
                 " refdb/ecopcr_uniq_clean",
                 " -o ", "amplified_clean_uniq.fasta"),          
          paste0("obi build_ref_db", 
                 " -t ", ecotag_min_similarity, 
                 " --taxonomy refdb/taxonomy/dump",
                 " refdb/ecopcr_uniq_clean", 
                 " refdb/ecopcr_final_", ecotag_min_similarity),
          paste0("obi export ",
                 " --fasta-output refdb/ecopcr_final_", ecotag_min_similarity,
                 " -o ", "final_db_", ecotag_min_similarity, ".fasta"),
          paste0("echo 'aantal unieke inputs na merging: ' >> logfile.txt"),
          paste0("grep -E  -i '>' ",
                 " final_db_", ecotag_min_similarity, ".fasta", 
                 " | wc -l >> logfile.txt"), #telt de lijnen die gegrept zijn
          paste0("\n")
    )
  if (script == "console") {
    script <- ""
    cat(commands, file = script)
  }
  #convert to unix enconding
  f <- file(script, open = "wb")
  cat(commands, file = f)
  close(f)
}


################################################################################

#' Title
#'
#' @param ecopcr 
#' @param merged 
#'
#' @return
#' @export
#'
#' @examples
combine_ecpocr_with_merged <- function(ecopcr, merged) {
  merged <- merged %>%   
    mutate(merged_count = count) %>% 
    filter(merged_count > 1)  
  
  combined <- ecopcr %>% 
    left_join(merged %>% 
                transmute(dna_hash, merged_count, 
                          merged_taxa = merged_taxid,
                          merged_rank = rank, obi_taxid = taxid,
                          merged_taxid, IS_MERGED = 1),
              by = "dna_hash")
  
  rv <- combined %>% 
    transmute(label = "INIT", genbank_id, strand,
              seq_len_input = seq_length_ori, 
              seq_len_amplicon = nchar(dna_sequence), 
              amplicon = dna_sequence, dna_hash,
              rank, taxid = as.numeric(taxid),
              species, species_name, 
              genus = as.numeric(genus), genus_name, 
              family = as.numeric(family), family_name, 
              forward_match, forward_error, forward_tm,
              reverse_match, reverse_error, reverse_tm, 
              is_merged = IS_MERGED, obi_count = as.numeric(merged_count), 
              obi_rank = merged_rank, obi_taxid = as.numeric(obi_taxid),
              merged_overview = merged_taxa)
  
  rv <- rv %>% 
    mutate(obi_count = ifelse(is.na(obi_count), 1, obi_count), 
           obi_rank = ifelse(is.na(obi_rank), rank, obi_rank), 
           obi_taxid = ifelse(is.na(obi_taxid), taxid, obi_taxid), 
           is_merged = ifelse(is.na(is_merged), 0, is_merged))
  
  rv
}

######################################################

check_multihits <- function(df, specieslist) {
  merges <- paste(df$merged_overview, collapse = ";")
  locs <- gregexpr("'", merges)[[1]]
  starts <- locs[seq(1, length(locs), by = 2)] + 1
  ends <- locs[seq(2, length(locs), by = 2)] - 1
  taxa <- NULL
  for (k in 1:length(starts)) {
    taxa <- rbind(taxa, 
                  data.frame(taxid = as.numeric(substring(merges, 
                                                          starts[k], ends[k]))))
  }
  taxa <- taxa %>% 
    distinct() %>% 
    left_join(specieslist %>% select(taxid, priority))
  
  if (all(is.na(taxa$priority))) {
    print(taxa)
    df$eval <- "OK"
    df$remove <- NA
    df$remark <- "No multihit between species from species list"
    return(df)
  } 
  
  taxa_ok <- taxa %>% 
    filter(!is.na(priority))
  print(taxa_ok)
  
  if (nrow(taxa_ok) < 2) {
    txid <- taxa_ok$taxid[1]
    df$eval <- NA
    df$eval[df$taxid == txid] <- "OK"
    df$eval[df$taxid != txid] <- "REMOVE"
    df$remove<- paste(taxa$taxid[taxa$taxid != txid], collapse = ";")
    df$remark <- "Multihit with a species not on specieslist "
    return(df)
  }
  
  n_p1 <- taxa %>% filter(priority == 1) %>% nrow()
  n_p2 <- taxa %>% filter(priority == 2) %>% nrow()
  
  if (n_p1 == 0) { #conflict tss prioriteit-2 soorten
    df$eval <- "P2_CLASH"
    df$remove<- paste(taxa$taxid, collapse = ";")
    df$remark <- "Multihit between priority-2 species"
    return(df)
  }
  
  if (n_p1 == 1) { #exact 1 prioriteit-1 soort --> kies deze
    txid <- taxa_ok$taxid[taxa_ok$priority == 1]
    df$eval <- NA
    df$eval[df$taxid == txid] <- "OK"
    df$eval[df$taxid != txid] <- "REMOVE"
    df$remove<- paste(taxa$taxid, collapse = ";")
    df$remark <- "Multihit between priority-1 and priority-2 species"
    return(df)
  }
  if (n_p1 > 1 & n_p2 > 0) { #conflict tss prioriteit-1 soorten en prioriteit-2
    txid <- taxa_ok$taxid[taxa_ok$priority == 1]
    df$eval <- NA
    df$eval[df$taxid %in% txid] <- "P1-CLASH"
    df$remove<- paste(taxa$taxid, collapse = ";")
    df$eval[!(df$taxid %in% txid)] <- "REMOVE"
    df$remark <- "Multihit between multiple p-1 species (with additional p-2 species)"
    return(df)
    
  }
  if (n_p1 > 1 & n_p2 == 0) {
    df$eval <- "P1_CLASH"
    df$remove<- paste(taxa$taxid, collapse = ";")
    df$remark <- "Multihit between p-1 species"
    return(df)
    
  }
  
  df$eval <- "not evaluated"
  df$remove <- NA
  df$remark <- "not caught in eval routine"
  return(df)
}

#####################################################################################
check_multihits2 <- function(data) {
  np1 <- sum(data$priority == 1)
  np2 <- sum(data$priority == 2)
  np3 <- sum(data$priority == 3)
  nrows <- nrow(data)
  data$EVAL <- NA
  data$NEWID <- NA
  data$NEWRANK <- NA
  if (np1 + np2 + np3 == 1) {
    data$EVAL <- "OK, unique"
    data$NEWID <- data$taxid
    data$NEWRANK <- data$rank   
  }
  if (np1 == 1 & np2 + np3 > 0) {
    whi1 <- which(data$priority == 1)
    ntx <- data$taxid[whi1]
    data$EVAL[whi1] <- "OK, Pref species"
    data$EVAL[-whi1] <- "multihit_with_p1"
    data$NEWID <- ntx
    data$NEWRANK[whi1] <- data$rank[whi1]
  }
  if (np1 > 1) {
    data$EVAL <- "CLASH P1"
    data$NEWID <- NA
    data$NEWRANK <- NA
  }
  if (np1 == 0 & np2 + np3 == 1 ){
    data$EVAL <- "OK, in this db"
    data$NEWID <- data$taxid
    data$NEWRANK <- data$rank
  }
  if (np1 == 0 & np2 + np3 > 1) {
    #to be implemented
  }
  data
}

########################################################

multihit_list_update <- function(df){
  if (df$eval == "OK") {
    pref_taxid <- df$taxid
    tmp <- data.frame(taxid = pref_taxid, pref_taxid = pref_taxid)
    rmtaxid <- as.numeric(unlist(str_split(df$remove, pattern = ";")))
    tmp2 <- data.frame(taxid = rmtaxid, pref_taxid = pref_taxid)
    rv <- bind_rows(tmp, tmp2)
    return(rv)
  }
  if (df$eval == "REMOVE") {
    pref_taxid <- -99
    rv <- data.frame(taxid = df$taxid, pref_taxid = pref_taxid)
    return(rv)
  } 
  
  pref_taxid <- -1
  tmp <- data.frame(taxid = df$taxid, pref_taxid = pref_taxid)
  rmtaxid <- unlist(str_split(df$remove, pattern = ";"))
  tmp$clash <- paste0(substring(df$eval,1,2), ": ", 
                      paste(rmtaxid, collapse = ";"))
  rv <- tmp
  return(rv)
}

######################################################

link_species_multihit_check <- function(df, db_name, specieslist){
  rv <- df %>% 
    left_join(specieslist %>% 
                select(taxid, taxid_name = NameScientific)) %>% 
    left_join(specieslist %>% 
                select(taxid, pref_taxid_name = NameScientific),
              by = c("pref_taxid" = "taxid")) %>% 
    transmute(label = db_name,  remark = clash, 
              taxid_name, taxid, 
              pref_taxid, pref_taxid_name, genbank_id)
  
  rv
}

#LABEL	REMARK	BLOCKED_ON	SCI_NAME	DUTCH_NAME	ENG_NAME	TAXID	PREF_TAXID	SCI_NAME_PREF

##########################################################



