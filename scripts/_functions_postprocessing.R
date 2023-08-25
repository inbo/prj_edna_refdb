
get_taxa_from_merged <- function(x, soortenlijst) {
  x <- paste(x$merged_overview, collapse = ",")
  y <- gsub("[{}]", "",x)
  y <- gsub("'", "", y)
  y <- unlist(strsplit(y, split = ","))
  y <- strsplit(y, ": ")
  df <- as.data.frame(do.call("rbind", y))
  colnames(df) <- c('taxid', 'aantal')
  df$taxid <- as.numeric(df$taxid)
  df$aantal <- as.numeric(df$aantal)
  df <- left_join(df, soortenlijst %>% select(taxid, priority, NameScientific, NameEnglish), by = "taxid")
  df <- df %>% 
    group_by(taxid, priority, NameScientific, NameEnglish) %>% 
    summarise(sum_aantal = sum(aantal), .groups = "drop")
}
#get_taxa_from_merged("{'9823': 1, '273789': 1, '41807': 1, '159856': 1}", df_soortenlijst)


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
