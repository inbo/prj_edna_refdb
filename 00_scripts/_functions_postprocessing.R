
get_taxa_from_merged <- function(x, soortenlijst) {
  #NA probleem met 1414750, omdat de rank 'no rank' is
  x <- paste(x$merged_overview, collapse = ",")
  y <- gsub("[{}]", "",x)
  y <- gsub("'", "", y)
  y <- unlist(strsplit(y, split = ","))
  y <- strsplit(y, ": ")
  df <- as.data.frame(do.call("rbind", y))
  if(ncol(df) == 1) df$V2 = 0 #om errors te voorkomen indien taxid NA
  colnames(df) <- c('taxid', 'aantal')
  df$taxid <- as.numeric(df$taxid)
  df$aantal <- as.numeric(df$aantal)
  df <- left_join(df, soortenlijst %>% select(taxid, priority, NameScientific, NameEnglish), by = "taxid")
  df <- df %>% 
    group_by(taxid, priority, NameScientific, NameEnglish) %>% 
    summarise(sum_aantal = sum(aantal), .groups = "drop")
}
#get_taxa_from_merged("{'9823': 1, '273789': 1, '41807': 1, '159856': 1}", df_soortenlijst)


###############################################################################################

judge_species <- function(df) {
  aantalp1 <- sum(df$priority == 1, na.rm = TRUE)
  aantalp2 <- sum(df$priority == 2, na.rm = TRUE)
  aantalp3 <- sum(is.na(df$priority) | df$priority %in% c(3:8))
  df$priority[is.na(df$priority)] <- 3
  df$oordeel <- NA
  if (aantalp1 > 0) {
    if (aantalp1 == 1) {
      pref_taxid <- df$taxid[df$priority == 1]
      df$oordeel[df$taxid == pref_taxid] <- "multihit, CHOSEN"
      df$oordeel[df$taxid != pref_taxid] <- "multihit, NOT chosen"
      df$pref_taxid <- pref_taxid
    } else if (aantalp1 > 1) {
      pref_taxid <- NA
      df$oordeel[df$priority == 1] <- "UNRESOLVED"
      df$oordeel[df$priority != 1] <- "multihit, NOT chosen"
      df$pref_taxid <- pref_taxid
    } 
    return(df)
  }
  if (aantalp2 == 1) {
    pref_taxid <- df$taxid[df$priority == 2]
    df$oordeel[df$taxid == pref_taxid] <- "multihit p2, CHOSEN"
    df$oordeel[df$taxid != pref_taxid] <- "multihit p2, NOT chosen"
    df$pref_taxid <- pref_taxid
  } else if (aantalp2 > 1) {
    pref_taxid <- NA
    df$oordeel[df$priority == 2] <- "UNRESOLVED p2"
    df$oordeel[df$priority != 2] <- "multihit p2, NOT chosen"
    df$pref_taxid <- pref_taxid
  } else {
    df$oordeel = "IGNORED, no p1 or p2"
    df$pref_taxid <- NA
  }
  return(df)
}

################################################################################

genereer_conflicten <- function(combined, soortenlijst, okmerges) {
  ecopcr_filtered_sp <- combined %>% 
    inner_join(soortenlijst %>% select(taxid, priority)) %>% 
    select(genbank_id, taxid, rank, priority, species_name,
           dna_hash, obi_rank, obi_taxid, genus_name, family_name,
           merged_overview, obi_count) %>% 
    filter(!(obi_rank %in% c("subspecies", "species", "species;"))) %>% 
    group_by(across(-genbank_id)) %>% 
    summarise(genbank_id = paste(genbank_id, collapse = ";"), .groups = "drop") %>% 
    arrange(obi_taxid, taxid, dna_hash)
  
  allowed_genera <- okmerges %>% dplyr::filter(.data$rank == "genus") %>% pull(.data$taxid)
  allowed_families <- okmerges %>% dplyr::filter(.data$rank == "family") %>% pull(.data$taxid)
  ecopcr_filtered_sp <- ecopcr_filtered_sp %>%  
    dplyr::filter(!(substring(obi_rank,1,5) == "genus" & obi_taxid %in% allowed_genera), 
                  !(substring(obi_rank,1,6) == "family" & obi_taxid %in% allowed_families))
  
  beoordeeld <- 
    ecopcr_filtered_sp %>%
    arrange(obi_taxid) %>%
    group_by(obi_taxid) %>%
    do ({
      df1 <- get_taxa_from_merged(., soortenlijst)
      #genbank_ids <- paste(.$genbank_id, collapse = ";")
      #bind_cols(df1, genbank_id = genbank_ids)
      df1
    }) %>% 
    group_by(obi_taxid) %>% 
    do({
      judge_species(.)    
    })
  
  conflicts <- beoordeeld %>% 
    group_by(taxid) %>%
    do({
      taxid <- .$taxid[1]
      obitaxids <- (unique(.$obi_taxid))
      taxids <- unique(beoordeeld$taxid[beoordeeld$obi_taxid %in% obitaxids])
      taxid_order <- order(taxids) 
      scinams <- unique(beoordeeld$NameScientific[beoordeeld$obi_taxid %in% obitaxids])
      data.frame(has_conflicts = TRUE, 
                 conflict_taxa = paste(taxids[taxid_order], collapse = '|'),
                 conflict_names = paste(scinams[taxid_order], collapse = '|'),
                 conflict_oordeel = paste(unique(.$oordeel), collapse = "|"))
    }) %>% 
    arrange(conflict_taxa)
  
  conflicts
}

#################################################################################


genereer_soortenevaluatie <- function(combined, merged, inputs, soortenlijst, multihits, okmerges, conflicts) {
  
  ### join de soortenlijst met de ecopcr-output en bereken de maximale merge rank
  ecopcr_soorten <- combined %>% 
    right_join(soortenlijst %>% select(taxid), by = "taxid") 
  
  ecopcr_soorten <- ecopcr_soorten %>% 
    group_by(taxid) %>% 
    do({
      merged = sum(.$is_merged)>0
      obi_ranks = factor(.$obi_rank, 
                         levels = c("subspecies", "species", "subgenus", "genus", "subfamily", "family",
                                    "no rank"))
      whimaxrank = which.max(obi_ranks)
      if (length(whimaxrank)) {
        corresptaxid = as.numeric(.$obi_taxid[whimaxrank])
        corresprank = as.character(.$obi_rank[whimaxrank])      
      } else {
        corresptaxid = as.numeric(.$obi_taxid[1])
        corresprank = as.character(.$obi_rank[1])       
      }
      merged = .$taxid[1] != corresptaxid
      rv <- data.frame(merged = merged, obi_taxid = corresptaxid, obi_rank = corresprank)
      rv
    })
  
  ### join met de multihit lijst zodat het pref_taxid ingevuld is, en de vlag is_chosen gezet wordt
  multihit_species <- multihits %>% 
    select(taxid = taxid, pref_taxid = pref_taxid) %>% 
    inner_join(soortenlijst %>% select(taxid)) %>% 
    mutate(is_multihit = TRUE,
           is_chosen = taxid == pref_taxid) 
  multihit_species <- multihit_species %>% 
    group_by(taxid) %>% 
    do({
      pref = .$pref_taxid[1]
      taxids <- sort(multihit_species$taxid[pref == multihit_species$pref_taxid])
      whi <- which(taxids == pref )
      if(length(whi)) {
        taxids = c(taxids[whi], taxids[-whi])
      }
      cbind(., hitlist = paste(taxids, collapse = " | "))
    })
  
  ### ok_merges (namen komen dan hieruit ipv multihitlijst)
  merges <- okmerges %>%
    left_join(combined %>%  select(obi_taxid, merged_overview),
              by = c("taxid" = "obi_taxid")) %>%
    distinct()
  
  merges <- merges %>% 
    group_by(taxid) %>%
    do({
      taxids <- paste(.$merged_overview, collapse = ",")
      taxids <- str_replace_all(taxids, "\\{", "")
      taxids <- str_replace_all(taxids, "\\}", "")
      taxids <- str_replace_all(taxids, "\\'", "")
      separated <- unlist(str_split(taxids, ", "))
      separated <- substring(separated, 1, regexpr(":", separated)-1)
      combined <- paste(sort(as.numeric(unique(separated))), collapse = " | ")
      data.frame(taxid_okm = .$taxid[1], rank = .$rank[1], name_okm = .$sci_name[1], hitlist_okm = combined)
    })
  
  ### voorkomen van soorten in de refdb als merged species zoek naar '#####' met ##### = taxid in de merged fasta file
  
  n_in_amplified <- soortenlijst %>%
    group_by(taxid) %>%
    do({
      aantal <- sum(str_detect(merged$merged_taxid, paste0("'", .$taxid[1], "'"))) #werkt niet voor obi2
      data.frame(n_merged_taxid = aantal)
    })

  ### Het genereren van de effectieve tabel
  soortenevaluatie <- soortenlijst %>% 
    group_by(Group, NameScientific, NameEnglish, NameDutch, taxid, rank, priority) %>% 
    summarise(n_input_orig = sum(inputs$taxid == .data$taxid[1]),
              n_obi_taxid = sum(merged$taxid == .data$taxid[1]),
              .groups = "drop")
    
  soortenevaluatie <- soortenevaluatie %>%   
    left_join(n_in_amplified, by = "taxid") %>% 
    
    left_join(ecopcr_soorten, by = "taxid") %>% 
    
    left_join(multihit_species, by = "taxid") %>% 
    mutate(pref_taxid = ifelse(is.na(pref_taxid) & (obi_taxid == taxid | obi_rank == "species"), 
                               taxid, 
                               pref_taxid)) %>%
    mutate(is_multihit = ifelse(!(obi_rank %in% c("species", "subspecies")), TRUE, is_multihit)) %>% 
    
    left_join(soortenlijst %>% select(taxid, 
                                      Pref_NameEnglish = NameEnglish, 
                                      Pref_NameScientific = NameScientific), 
              by = c("pref_taxid" = "taxid"))
    
    soortenevaluatie <- soortenevaluatie %>% 
        left_join(merges %>% transmute(taxid, taxid_okm, name_okm, hitlist_okm), 
              by = c("obi_taxid" = "taxid")) %>% 
    mutate(Pref_NameScientific = ifelse(is.na(Pref_NameScientific), name_okm, Pref_NameScientific),
           pref_taxid = ifelse(is.na(hitlist), taxid_okm, pref_taxid),
           hitlist = ifelse(is.na(hitlist), hitlist_okm, hitlist)
    )
    
  soortenevaluatie <- soortenevaluatie %>% 
    left_join(conflicts, by = "taxid") %>% 
    
    select(Group, priority, n_input_orig, n_merged_taxid, n_obi_taxid,
           NameScientific, NameEnglish, NameDutch,
           is_multihit, is_chosen, Pref_NameScientific, 
           rank_obi_taxid = obi_rank, obi_taxid, rank, taxid, pref_taxid,
           hitlist, Pref_NameEnglish, has_conflicts, conflict_taxa, conflict_oordeel, conflict_names) %>% 
    arrange(Group, pref_taxid, desc(is_chosen))
  
  soortenevaluatie %>% distinct()
}

####################################################################################


