filter_indices <- function(model, shortmodel, met_index, max_rxns){
  # take only the metabolites NOT exceeding the number of reactions and take the index of those in the original model (recon3D) datamatrix
  # This script takes into account that multiple rows may contain the same metabolite
  metNames <- unique(as.vector(unlist(model$metNames[met_index])))
  metNames_u_thres <- metNames[rowSums(shortmodel[metNames,,drop = FALSE] != 0) <= max_rxns]
  indices_filt <- unlist(lapply(metNames_u_thres, function(x) which(x == as.vector(unlist(model$metNames)))))
  return(indices_filt)
}



