getPreviousNext <- function(model, met_index, max_rxns, step, path){

  mets=NULL
  # Get a fraction of the number of reactions
  max_rxns <- fr_max_rxns*ncol(model$S) 
  test <- data.table("name" = as.vector(unlist(model$metNames[met_index])), 
                     "index"=met_index, 
                     "rxns" = rowSums(model$S[met_index,] != 0)
                     )
  setkey(test,name)
  test2 <- test[, .(sum = sum(rxns)), by = name]
  if(!all(test2$sum > max_rxns)){
    return(NULL)
  }
  
  
  for(j in 1:length(met_index)){

    
    met=as.vector(unlist(model$mets[met_index[j]], use.names = FALSE)) 
    
    # find in which reaction used
    rxn_index = try({which(model$S[as.numeric(met_index[j]),]!=0)},silent = TRUE)
    if (class(rxn_index) == "try-error") {
      next
    }  
    
    label = "used in: "
    
    # if (length(rxn_index)>0){  
    # Metabolites present in more than max_rxns reactions are not executed here
    if (max_rxns>=length(rxn_index)&length(rxn_index)>0){ 
      
      # message(paste("Yes",j, met))
      
      for(k in 1:length(rxn_index)){
        
        # message(paste(label, as.vector(unlist(model$rxnNames[rxn_index[k]], use.names = FALSE)), sep=""))
        
        level=step
        
        tmp = getMetsRecon2V2(rxn_index[k], level, met, model)
        mets_sub = tmp$metabolites
        

        if (is.na(path[j])){
          begin = as.vector(unlist(model$metNames[met_index[j]], use.names = FALSE))
        } else {
          begin = path[j]
        }
        
        rxns=getReactionsRecon(index = tmp$mets, model)
        mets_sub = cbind(mets_sub, "resource"=rep("Recon",dim(mets_sub)[1]), "rxn_formula"=rxns[,"rxns"], "path"=paste(begin,
                                                                                                                       " =>(", as.vector(unlist(model$rxnNames[rxn_index[k]], use.names = FALSE)),")=>",mets_sub[,"met_long"]))
        mets = rbind(mets, mets_sub)
      }
    } 
  }
  
  if (is.null(mets)){
    return(NULL)
  } else {
    return(list("mets"=mets))
  }
}
