#######################################
### HELPER FUNCTIONS FOR EVALUATION ###
#######################################

# Conflicts
# Show conflicting amplicons from all taxids
# Need to have the "ecopcr_combined" loaded in the env
conflict_seqs = function(conflict_string, ecopcr_combined){
  mycols=c(2,4,5,6,7,9,10,11, 22:26)
  conf_taxids=unlist(strsplit(conflict_string,split = "\\|"))
  xx = ecopcr_combined[ecopcr_combined$taxid %in% conf_taxids, mycols]
  return(xx)
}

# Examples
View(conflict_seqs("150711|318395", ecopcr_combined))
View(conflict_seqs("446431|446432", ecopcr_combined))
View(conflict_seqs("356229|635147", ecopcr_combined))
