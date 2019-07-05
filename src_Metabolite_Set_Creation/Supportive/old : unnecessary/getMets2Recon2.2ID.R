getMets2Recon2.2ID <- function(mets, recon2chebi){
#mets = consumed
  
  chebi_recon = recon2chebi[,"CHEBI"]
  chebi_kegg = mets[,"CheBI_id"]
  
  index.tmp = which(is.na(chebi_kegg))
  if (length(index.tmp)>0){
    chebi_kegg = chebi_kegg[-index.tmp]  
  }
  
  index=which(chebi_recon %in% chebi_kegg) 
  
  recon2chebi[index,]
  
  return(index)
}  