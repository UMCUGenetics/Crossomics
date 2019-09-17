getMets2ReconID <- function(mets, model){
  # mets = consumed
  
  chebi = lapply(model$metCHEBIID, function(x) {
    
    if (length(unlist(x))==0) {
      return(NA)
    } else {
      return(as.vector(unlist(x)))
    }
    
  })
  
  kegg = lapply(model$metKeggID, function(x) {
    
    if (length(unlist(x))==0) {
      return(NA)
    } else {
      return(as.vector(unlist(x)))
    }
    
  })

  pubchem = lapply(model$metPubChemID, function(x) {
    
    if (length(unlist(x))==0) {
      return(NA)
    } else {
      return(as.vector(unlist(x)))
    }
    
  })
  
  chebi_recon = as.vector(unlist(chebi))
  
  # chebi_kegg = unlist(lapply(strsplit(mets[,"CheBI_id"], " ", fixed = T), function(x) x[2]))
  chebi_kegg = mets[,"CheBI_id"]

  index.tmp = which(is.na(chebi_kegg))
  if (length(index.tmp)>0){
    chebi_kegg = chebi_kegg[-index.tmp]  
  }

  index=which(chebi_recon %in% chebi_kegg) 

  #####################################
  
  kegg_recon = as.vector(unlist(kegg))
  
  kegg_kegg = mets[,"KEGG_id"]
  
  index.tmp = which(is.na(kegg_kegg))
  if (length(index.tmp)>0){
    kegg_kegg = kegg_kegg[-index.tmp]  
  }

  index2=which(kegg_recon %in% kegg_kegg)

  
    
  ####################################

  pubchem_recon = as.vector(unlist(pubchem))
  
  pubchem_kegg = mets[,"PubChem_id"]
  
  index.tmp = which(is.na(pubchem_kegg))
  if (length(index.tmp)>0){
    pubchem_kegg = pubchem_kegg[-index.tmp]  
  }
  
  index3=which(pubchem_recon %in% pubchem_kegg)

  return(union(union(index, index2), index3))
}    
