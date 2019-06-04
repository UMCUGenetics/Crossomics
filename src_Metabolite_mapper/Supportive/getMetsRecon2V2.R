# From gene to metabolite
getMetsRecon2V2 <- function(rxn_index, level, met, model){
# rxn_index=rxn_index_consumed[k]

  S.sub = model$S[,rxn_index]
  if (length(S.sub)>0) {
    met_index_left=which(S.sub<0)
    met_index_right=which(S.sub>0)
  }
  
  rxn_long = as.vector(unlist(model$rxnNames[rxn_index], use.names = FALSE))
  rxn_ss = as.vector(unlist(model$subSystems[rxn_index], use.names = FALSE))
  
  met_long_left = as.vector(unlist(model$metNames[met_index_left], use.names = FALSE)) # 5062 ????????????????????
  met_short_left = as.vector(unlist(model$mets[met_index_left], use.names = FALSE))
  met_long_right = as.vector(unlist(model$metNames[met_index_right], use.names = FALSE)) # 5062 ????????????????????
  met_short_right = as.vector(unlist(model$mets[met_index_right], use.names = FALSE))
  
  getMets <- function(index, label, level, rxn_index, met, rxn_long, rxn_ss, met_short, met_long) {
  
#     index =met_index_right
#     label="right"
#     met_short=met_short_right
#     met_long=met_long_right
    
    if (is.null(met_short)) met_short=NA
    if (is.null(met_long)) met_long=NA
    
    hmdb = model$metHMDB[index]
    tmp = rep(NA, length(hmdb))
    if (length(hmdb)>0){
      for (i in 1:length(hmdb)){
        if (identical(as.character(hmdb[[i]]), character(0))) tmp[i] = NA else tmp[i] = as.character(hmdb[[i]])   
      }
      hmdb=tmp
    } else {
      hmdb=NA
    }

    kegg = model$metKeggID[index]
    tmp = rep(NA, length(kegg))
    if (length(kegg)>0){
      for (i in 1:length(kegg)){
        if (identical(as.character(kegg[[i]]), character(0))) tmp[i] = NA else tmp[i] = as.character(kegg[[i]])   
      }
      kegg=tmp
    } else {
      kegg=NA
    }
    
    chebi = model$metCHEBIID[index]
    tmp = rep(NA, length(chebi))
    if (length(chebi)>0){
      for (i in 1:length(chebi)){
        if (identical(as.character(chebi[[i]]), character(0))) tmp[i] = NA else tmp[i] = as.character(chebi[[i]])   
      }
      chebi=tmp
    } else {
      chebi=NA
    }
    
    pubchem = model$metPubChemID[index]
    tmp = rep(NA, length(pubchem))
    if (length(pubchem)>0){
      for (i in 1:length(pubchem)){
        if (identical(as.character(pubchem[[i]]), character(0))) tmp[i] = NA else tmp[i] = as.character(pubchem[[i]])   
      }
      pubchem=tmp
    } else {
      pubchem=NA
    }
    
    rxn_id = rep(rxn_index, length(hmdb))
    step = rep(level, length(hmdb))
    met_in = rep(met, length(hmdb))
    left_right = rep(label, length(hmdb))
    rxn_name = rep(rxn_long, length(hmdb))
    rxn = rep(rxn_ss, length(hmdb))
    
    return(list("metabolites"=cbind(rxn_id, step, met_in, left_right,  met_short, met_long, hmdb, kegg, chebi, pubchem, rxn_name, rxn), "rxns"=rxn_index, "mets"=index))

  }
  
  left = getMets(met_index_left,"left", level, rxn_index, met, rxn_long, rxn_ss, met_short_left, met_long_left)
  right = getMets(met_index_right,"right",level, rxn_index, met, rxn_long, rxn_ss, met_short_right, met_long_right)
  
  return(list("metabolites"=rbind(left$metabolites, right$metabolites), "matrix"=model$S[c(met_index_left, met_index_right), rxn_index], "rxns"=rxn_index, "mets"=c(met_index_left, met_index_right)))
  
}
