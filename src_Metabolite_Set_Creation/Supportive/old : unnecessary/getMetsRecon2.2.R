# From gene to metabolite
getMetsRecon2.2 <- function(rxn_index, level, met, model, recon2chebi){
# rxn_index=rxn_index[k]
# level=1

  S.sub = model@S[,rxn_index]
  if (length(S.sub)>0) {
    met_index_left=which(S.sub<0)
    met_index_right=which(S.sub>0)
  }
  
  # which(model@subSys[rxn_index,]==TRUE)
  
  rxn_long = as.vector(unlist(model@react_name[rxn_index], use.names = FALSE))
  rxn_ss = names(model@subSys[rxn_index,])[which(model@subSys[rxn_index,]==TRUE)]
  
  met_long_left = as.vector(unlist(model@met_name[met_index_left], use.names = FALSE)) # 5062 ????????????????????
  met_short_left = as.vector(unlist(model@met_id[met_index_left], use.names = FALSE))
  met_long_right = as.vector(unlist(model@met_name[met_index_right], use.names = FALSE)) # 5062 ????????????????????
  met_short_right = as.vector(unlist(model@met_id[met_index_right], use.names = FALSE))
  
  getMets <- function(index, label, level, rxn_index, met, rxn_long, rxn_ss, met_short, met_long) {
#     index =met_index_left
#     label="left"
#     met_short=met_short_left
#     met_long=met_long_left
    
    if (is.null(met_short)) met_short=NA
    if (is.null(met_long)) met_long=NA
    
    chebi = recon2chebi[index,2]
    tmp = rep(NA, length(chebi))
    if (length(chebi)>0){
      for (i in 1:length(chebi)){
        if (identical(as.character(chebi[[i]]), character(0))) tmp[i] = NA else tmp[i] = as.character(chebi[[i]])   
      }
      chebi=tmp
    } else {
      chebi=NA
    }
    
    InChI = recon2chebi[index,3]
    tmp = rep(NA, length(InChI))
    if (length(InChI)>0){
      for (i in 1:length(chebi)){
        if (identical(as.character(InChI[[i]]), character(0))) tmp[i] = NA else tmp[i] = as.character(InChI[[i]])   
      }
      InChI=tmp
    } else {
      InChI=NA
    }

    rxn_id = rep(rxn_index, length(index))
    step = rep(level, length(index))
    met_in = rep(met, length(index))
    left_right = rep(label, length(index))
    rxn_name = rep(rxn_long, length(index))
    rxn = rep(rxn_ss, length(index))
    
    return(list("metabolites"=cbind(rxn_id, step, met_in, left_right,  met_short, met_long, chebi, InChI, rxn_name, rxn), "rxns"=rxn_index, "mets"=index))

  }
  
  left = getMets(met_index_left,"left", level, rxn_index, met, rxn_long, rxn_ss, met_short_left, met_long_left)
  right = getMets(met_index_right,"right",level, rxn_index, met, rxn_long, rxn_ss, met_short_right, met_long_right)
  
  if (length(met_index_left)==0 & length(met_index_right)==0){
    return(NULL)
  } else if (length(met_index_left)==0) {
    return(list("metabolites"=right$metabolites, "matrix"=model@S[met_index_right, rxn_index], "rxns"=rxn_index, "mets"=met_index_right))
  } else if (length(met_index_right)==0) {
    return(list("metabolites"=left$metabolites, "matrix"=model@S[met_index_left, rxn_index], "rxns"=rxn_index, "mets"=met_index_left))
  } else {  
    return(list("metabolites"=rbind(left$metabolites, right$metabolites), "matrix"=model@S[c(met_index_left, met_index_right), rxn_index], "rxns"=rxn_index, "mets"=c(met_index_left, met_index_right)))
  }
}
