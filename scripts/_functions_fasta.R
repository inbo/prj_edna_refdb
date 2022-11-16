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
  tmpmeta <-  substring(metarow, spatie+1) #data met recordnaam eruit gestript
  
  #>>> Vind Taxid
  taxidpos <- regexpr("taxid=", tmpmeta)
  seppos_mult <- gregexpr(";", tmpmeta)[[1]]
  seppos <- min(seppos_mult[seppos_mult > taxidpos+6]) #eerstvolgend na taxid
  taxid <- NULL #initial value
  if (taxidpos > 0) {
    taxid <- substring(tmpmeta, taxidpos+6, seppos-1)
  }
  
  #>>> Nieuws: vind count
  
  #>>> bepaal de genlab recordnaam en merged count
  if (is_merged) {
    #indien gemerged staat count= bij de recordnaam, dus hier ophalen
    cntpos <- regexpr("count=", tmpmeta)
    seppos_mult <- gregexpr(";", tmpmeta)[[1]]
    seppos <- min(seppos_mult[seppos_mult > cntpos+6])  
    merged_count <- substring(tmpmeta, cntpos+6, seppos-1)
    if (cntpos > 0) {
      out <- c(out,  'merged_count' = merged_count)     
    }
  }
  out <- c(out,  'genlab_id' = recordname)
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
      cat("GENLAB_ID:", parsed_record$genlab_id, "\n")
      cat("TAXID:    ", parsed_record$taxid, "\n")
      cat("LEN:      ", nchar(parsed_record$dna_sequence), "\n")
      cat("SEQ:      ", substring(parsed_record$dna_sequence, 1, 30), "...\n")
    }
    if (i == 1) cat("\nimporting", nrow(recordlist), "records:\n")
    refdb <- bind_rows(refdb, parsed_record)
  }
  
  ### >>> VOEG DE FILES SAMEN
  for ( i in 1:nrow(refdb)) {
    refdb$AMPLICON_HASH[i] <- sha1hash(refdb$dna_sequence[i])
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
                               entry_id = "genlab_id") {
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
