#' Title
#'
#' @param x 
#'
#' @return
#' @export
#'
#' @examples
sha1hash <- function(x) {
  rv <- character(length(x))
  for(i in 1:length(x)) {
    rv[i] <- digest::sha1(x[i])
  }
  rv
}

#################################################################


#' Title
#'
#' @param record 
#' @param is_merged 
#'
#' @return
#' @export
#'
#' @examples
#' fr <- NULL
#' fr[1] <- ">AF038468.1 count=6; merged_taxid={58317: 6}; species_name=Blicca bjoerkna; family=2743726; family_name=Leuciscidae; scientific_name=Blicca bjoerkna; reverse_match=CTTCCGGTACACTTACCATG; forward_error=0; rank=species; taxid=58317; forward_tm=60.26; genus_name=Blicca; forward_match=ACACCGCCCGTCACTCT; reverse_tm=nan; genus=58316; reverse_error=0; species=58317; strand=D;"
#' fr[2] <- "cccctgtcaaaatgcaataaagttacttaacaccaaagcgctgacaaggggaggcaagtc" 
#' fr[3] <- "gtaa"
#' testrecord <- parse_fasta_record(fr, is_merged = TRUE, dna_type = "uppercase")
#' 
#' rec2 <- NULL
#' rec2[1] <- ">MF326939 taxid=214916;"
#' rec2[2] <- "CAAAGGCTTGGTCCTGACTTTATTATCAACTTTAGCCAAACTT"
#' testrecord2 <- parse_fasta_record(rec2, is_merged = FALSE)

parse_fasta_record <- function(record, is_merged, dna_type = "uppercase") {
  out <- NULL
  
  #>>> verkrijg de sequentie
  sequentie <- paste(record[2:length(record)], collapse = "")
  if (dna_type == "uppercase") {
    sequentie <- toupper(sequentie)
  } else if (dna_type == "lowercase") {
    sequentie <- tolower(sequentie)
  } else {
    #don't change anything
  }
  
  #>>> verkrijg de metadatarij
  metarow <- record[1]
  spatie <- regexpr("\\s", metarow) #vind de eerste spatie
  merged_count <- -1
  recordname <- substring(metarow, 2, spatie - 1) 
  tmpmeta <-  substring(metarow, spatie) #data met recordnaam eruit gestript, spatie moet behouden blijven
  
  #>>> Vind Taxid
  taxidpos <- regexpr(" taxid=", tmpmeta)
  seppos_mult <- gregexpr(";", tmpmeta)[[1]]
  seppos <- min(seppos_mult[seppos_mult > taxidpos+7]) #eerstvolgend na taxid
  taxid <- NULL #initial value
  if (taxidpos > 0) {
    taxid <- substring(tmpmeta, taxidpos+7, seppos-1)
  }
  
  #>>> Nieuws: vind count
  
  #>>> bepaal de genlab recordnaam en merged count
  if (is_merged) {
    #indien gemerged staat count= bij de recordnaam, dus hier ophalen
    cntpos <- regexpr("count=", tmpmeta)
    if (cntpos == -1) cntpos <- regexpr("COUNT=", tmpmeta)
    seppos_mult <- gregexpr(";", tmpmeta)[[1]]
    seppos <- min(seppos_mult[seppos_mult > cntpos+6])  
    merged_count <- substring(tmpmeta, cntpos+6, seppos-1)
    if (cntpos > 0) {
      out <- c(out,  'merged_count' = merged_count)     
    }
  }
  out <- c(out,  'genbank_id' = recordname)
  out <- c(out, "taxid" = taxid)
  out <- c(out, 'dna_sequence' = sequentie)

  
  #>>> splits de metadata op in verschillende velden
  splits <- str_split(tmpmeta, "; ", simplify = FALSE)[[1]]  
  for (j in 1:length(splits)) {
    piece <- str_trim(splits[j])
    sep <- regexpr("=", piece) #de recordnaam en waarde zijn gescheiden door =
    if (sep > 0) {
      keyname <- substring(piece, 1, sep-1)
      value <- substring(piece, sep + 1)
      if (!(keyname %in% c("taxid", "count"))) { #taxid, count wordt al eerder toegekend
        out <- c(out, value )
        names(out)[length(out)] <- keyname        
      }
    }
  }
  out
}

###########################################


#' Title
#'
#' @param file 
#' @param is_merged_file 
#'
#' @return
#' @export
#'
#' @examples
parse_refdb_fasta <- 
  function(file, is_merged_file = FALSE, dna_type = "uppercase") {
  ### >>>lees de file lijn per lijn in
  content <- readLines(file)
  lastslash <- max(gregexpr("/", file)[[1]])
  filestem <- substring(file, lastslash+1)
  inputfilestem <- gsub(".fasta", "", filestem)
  
  #content: vind waar een record start
  newrecords <- which(startsWith(content, ">"))
  endrecords <- c(newrecords[-1] - 1, length(content))
  recordlist <- data.frame(seq_nr = 1:length(newrecords), 
                           start = newrecords, 
                           end = endrecords)
  
  #>>> loop door alle records en haal de informatie op
  refdb <- NULL
  for (i in 1:nrow(recordlist)) {
    divs <- nrow(recordlist) %/% 50
    if (divs == 0) cat("=") else if (i %% divs == 0) cat("=")
    record <- content[recordlist[i,"start"]:recordlist[i,"end"]]
    parsed_record <- parse_fasta_record(record, 
                                        is_merged = is_merged_file,
                                        dna_type =  dna_type)
    parsed_record <- as.data.frame(t(parsed_record))
    parsed_record$source <- inputfilestem
    whi <- which(colnames(parsed_record) == "merged_count")
    if (length(whi)) {
      parsed_record$merged_count <- as.character(parsed_record$merged_count)
    }
    if (i == 1) { 
      test <<- parsed_record
      cat("\n\nfirst record parsed:\n-------------------------\n")
      cat("GENBANK_ID:", parsed_record$genbank_id, "\n")
      cat("TAXID:    ", parsed_record$taxid, "\n")
      cat("LEN:      ", nchar(parsed_record$dna_sequence), "\n")
      cat("SEQ:      ", substring(parsed_record$dna_sequence, 1, 30), "...\n")
    }
    if (i == 1) cat("\nimporting", nrow(recordlist), "records:\n")
    refdb <- bind_rows(refdb, parsed_record)
  }
  
  ### >>> VOEG DE FILES SAMEN
  for ( i in 1:nrow(refdb)) {
    refdb$DNA_HASH[i] <- sha1hash(refdb$dna_sequence[i])
  }
  cat("\nDONE READING", nrow(refdb), "OF", nrow(recordlist),  "RECORDS\n\n")
  refdb
}

#####################################################################

#' Title
#'
#' @param file 
#' @param data 
#'
#' @return
#' @export
#'
#' @examples
create_input_fasta <- function(file, data, lowercase = FALSE, 
                               seq = "dna_sequence",
                               taxid = "taxid",
                               entry_id = "genbank_id") {
  if(lowercase) {
    data$SEQ_12S <-  tolower(data[[seq]])
  } else {
    data$SEQ_12S <-  toupper(data[[seq]])
  }
  data$txt <- paste0(">", data[[entry_id]], 
                     " taxid=", data[[taxid]], ";\n", 
                     data$SEQ_12S,
                     "\n")
  alltxt <- paste(data$txt, collapse = "")
  f <- file(file, open = "wb") #open as binary (geen CRLF maar LF)
  cat(alltxt, file = f, append = FALSE)
  close(f)
}

#########################################################################

#' Title
#'
#' @param which 
#'
#' @return
#' @export
#'
#' @examples
select_primer_tags <- function(which = "RIAZ") {
  if (which == "RIAZ") {
    p1 <- "ACTGGGATTAGATACCCC"
    p2 <- "TAGAACAGGCTCCTCTAG"
    p2rc <- reverse_complement(p2)
    p1rc <- reverse_complement(p1)   
  } else if (which == "TELEO") {
    p1 <- "ACACCGCCCGTCACTCT"
    p2 <- "CTTCCGGTACACTTACCRTG"
    p2rc <- reverse_complement(p2)
    p1rc <- reverse_complement(p1)     
  } else {
    stop("not implemented")
  }

  c(oligo1 = p1, oligo2 = p2, oligo2rc = p2rc, oligo1rc = p1rc)
}

###########################################################################

#' Title
#'
#' @param x 
#'
#' @return
#' @export
#'
#' @examples
reverse_complement <- function(x) {
  rv <- NULL
  for (i in 1:length(x)){
    dna <- Biostrings::DNAString(x[i])
    rv[i] <- as.character(Biostrings::reverseComplement(dna))
  } 
  rv
}

#############################################################################


#' Title
#'
#' @param size 
#' @param prob 
#' @param preamble 
#' @param postamble 
#' @param h1 
#' @param h2 
#' @param p1 
#' @param p2 
#'
#' @return
#' @export
#'
#' @examples
create_random_seq <- function(size = 100, prob = c(A=0.25, C=0.25, G=0.25, T=0.25),
                              preamble = 10,  postamble = 10, 
                              h1 = NULL, h2 = NULL,
                              p1 = select_primer_tags("TELEO")[1],
                              p2 = select_primer_tags("TELEO")[2]) {
  #R = GA
  #Y = TC
  #N = GATC
  pre = paste(sample(c("A", "C", "G", "T"), size = preamble, replace = TRUE, prob = prob), collapse = "")
  post = paste(sample(c("A", "C", "G", "T"), size = postamble, replace = TRUE, prob = prob), collapse = "")
  dna =  paste(sample(c("A", "C", "G", "T"), size = size, replace = TRUE, prob = prob), collapse = "")
  seq <- paste0(pre, h1, p1, dna, p2, h2, post)
  len <- nchar(seq)
  seq <- str_split_1(seq, "")
  R <- sample(c("G", "A"), size = len, replace = TRUE)
  Y <- sample(c("T", "C"), size = len, replace = TRUE)
  N <- sample(c("A", "C", "G", "T"), size = size, replace = TRUE)
  seq_fin <- character(len)
  seq_fin[seq == "A"] <- "A"
  seq_fin[seq == "C"] <- "C"
  seq_fin[seq == "G"] <- "G"
  seq_fin[seq == "T"] <- "T"
  seq_fin[seq == "R"] <- R[seq == "R"]
  seq_fin[seq == "Y"] <- R[seq == "Y"]
  seq_fin[seq == "N"] <- R[seq == "N"]
  seq_fin <- paste(seq_fin, collapse = "")
  seq_fin
}

############################################################################

#' Title
#'
#' @param dna 
#' @param which 
#'
#' @return
#' @export
#'
#' @examples
find_matching_seq <- function(dna, which = "TELEO") {
  dna <- toupper(dna)
  if (toupper(which) == 'RIAZ') {
    R1_match <- regexpr(pattern = "ACTGGGATTAGATACCCC", text = dna)
    R2_match <- regexpr(pattern = "CTAGAGGAGCCTGTTCTA", text = dna)
    rv <- paste0(as.numeric(R1_match>0), "_", as.numeric(R2_match>0), 
                 "_len:",R2_match - R1_match )
  } else if (toupper(which) == "TELEO") {
    T1_match <- regexpr(pattern = "ACACCGCCCGTCACTCT",  text = dna)
    T2a_match <- regexpr(pattern = "CATGGTAAGTGTACCGGAAG", text = dna)
    T2b_match <- regexpr(pattern = "CACGGTAAGTGTACCGGAAG", text = dna)
    rv <- paste0(as.numeric(T1_match>0), "_", 
                 max(T2a_match>0, T2b_match>0), 
                 "_len:", max(T2a_match, T2b_match) - T1_match + 1)
  } else {
    stop('which is not TELEO or RIAZ')
  }
  rv
}



