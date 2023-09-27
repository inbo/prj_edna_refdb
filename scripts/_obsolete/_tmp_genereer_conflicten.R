
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
  
  ecopcr_filtered_sp <- ecopcr_filtered_sp %>%  
    dplyr::filter(!(substring(obi_rank,1,5) == "genus" & 
                      obi_taxid %in% (okmerges %>% dplyr::filter(rank == "genus") %>% pull(taxid))), 
                  !(substring(obi_rank,1,6) == "family" & 
                      obi_taxid %in% (okmerges %>% dplyr::filter(rank == "family") %>% pull(taxid))  ))
  
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
