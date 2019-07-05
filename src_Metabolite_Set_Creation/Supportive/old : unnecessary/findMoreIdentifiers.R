findMoreIdentifiers <- function(metaboliteSet){
  library("BridgeDbR")
  
  # A lot of missing HMDB identifiers!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (!is.null(metaboliteSet)){
    if (is.null(dim(metaboliteSet))){
      metaboliteSet=data.frame(t(metaboliteSet))
    }
  }  
  
  ##################################################
  # temporarily work around to be fixed in findMetabolicEnvironment
  index = which(is.na(metaboliteSet[,"hmdb"]))
  if (length(index)>0) metaboliteSet[index,"hmdb"] = "character(0)"
  ##################################################
  
  ########### Recon 2.0 ############################
  index = which(metaboliteSet[,"hmdb"] == "character(0)") 
  
  # pressent in BridgeDB?
  index.sub = which(metaboliteSet[index,"kegg"] != "character(0)")
  kegg_id = metaboliteSet[index[index.sub],"kegg"]
  
  mapper = loadDatabase("./db/BridgeDB/metabolites_20150717.bridge")
  hmdb = getSystemCode("HMDB")
  kegg = getSystemCode("KEGG Compound")
  
  if (length(kegg_id)>0){
    for (k in 1:length(kegg_id)){
      if (!is.null(unlist(map(mapper, kegg, kegg_id[k], hmdb)[1]))) metaboliteSet[index[index.sub[k]],"hmdb"] = unlist(map(mapper, kegg, kegg_id[k], hmdb)[1]) 
    }
  }  
  
  index = which(metaboliteSet[,"hmdb"] == "character(0)")
  index.sub = which(metaboliteSet[index,"chebi"] != "character(0)")
  chebi_id = metaboliteSet[index[index.sub],"chebi"]
  
  chebi = getSystemCode("ChEBI")
  
  if (length(chebi_id)>0){
    for (k in 1:length(chebi_id)){
      if (!is.null(unlist(map(mapper, chebi, chebi_id[k], hmdb)[1]))) metaboliteSet[index[index.sub[k]],"hmdb"] = unlist(map(mapper, chebi, chebi_id[k], hmdb)[1]) 
    }
  }
  
  index = which(metaboliteSet[,"hmdb"] == "character(0)") 
  if (length(index)>0) metaboliteSet = metaboliteSet[-index,,drop=FALSE]
  #######################################################

  return(metaboliteSet)
  
}