# Extention of metabolites when an extention has already been performed
further_extend_metabolites <- function(model, shortmodel, mylist, result_mets, max_rxns, previous_index_removed, step){
  index <- mylist
  index_removed <- mylist
  index_removed2 <- mylist
  for(i in 1:length(max_rxns)){
    if (nrow(result_mets[[i]]) !=0){
      index[[i]] <- which(as.vector(unlist(model$mets)) %in% as.vector(result_mets[[i]][,"met_short"]))
      index_removed[[i]] = removeMetsFromSet(index[[i]],model)
    }
  }
  
  index_removed2 <- lapply(1:length(max_rxns), function(x) {
    if(nrow(result_mets[[x]]) != 0){
      min_rxn <- ifelse(x==1, 0, max_rxns[x-1])
      max_rxn <- max_rxns[x]
      tmpindex <- filter_indices(model, shortmodel, index_removed[[x]], max_rxns = max_rxn, min_rxns = min_rxn)
      # Also skip any indices that were already extended in previous step
      tmpindex <- tmpindex[!tmpindex %in% previous_index_removed[[x]]]
      if(length(tmpindex) != 0 ) tmpindex
    } 
  })
  
  result_mets <- extend_metabolites(model, max_rxns, index = index_removed2, mylist = mylist, step = step)
  
  return(list("result_mets"=result_mets, "index_removed" = index_removed2))
  
}