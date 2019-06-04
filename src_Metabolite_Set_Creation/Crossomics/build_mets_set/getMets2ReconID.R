getMets2ReconID <- function(mets, model){
# mets = result_mets
  
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
  
  chebi_set = mets[,"chebi"]
  
  index.tmp = which(is.na(chebi_set))
  if (length(index.tmp)>0){
    chebi_set = chebi_set[-index.tmp]  
  }

  index=which(chebi_recon %in% chebi_set) 

  #####################################
  
  kegg_recon = as.vector(unlist(kegg))
  
  kegg_set = mets[,"kegg"]
  
  index.tmp = which(is.na(kegg_set))
  if (length(index.tmp)>0){
    kegg_set = kegg_set[-index.tmp]  
  }

  index2=which(kegg_recon %in% kegg_set)

    
  ####################################

  pubchem_recon = as.vector(unlist(pubchem))
  
  pubchem_set = mets[,"pubchem"]
  
  index.tmp = which(is.na(pubchem_set))
  if (length(index.tmp)>0){
    pubchem_set = pubchem_set[-index.tmp]  
  }
  
  index3=which(pubchem_recon %in% pubchem_set)

  return(union(union(index, index2), index3))
}    
