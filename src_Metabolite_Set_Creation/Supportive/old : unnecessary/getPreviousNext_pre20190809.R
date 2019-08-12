getPreviousNext <- function(model, met_index, step, path){

  mets=NULL
  
  for(j in 1:length(met_index)){
    if(j == 1) cat("index out of", length(met_index), ": 1")
    if(j%%10 == 0 ) cat(", ", j)

    met=as.vector(unlist(model$mets[met_index[j]], use.names = FALSE)) 
    
    # find in which reaction used
    rxn_index = try({which(model$S[as.numeric(met_index[j]),]!=0)},silent = TRUE)
    if (class(rxn_index) == "try-error") {
      next
    }  
    
    label = "used in: "
    
    if (length(rxn_index)>0){ 
      
      for(k in 1:length(rxn_index)){
        
        level=step
        
        tmp = getMetsRecon2V2(rxn_index[k], level, met, model)
        mets_sub = tmp$metabolites
        

        if (is.na(path[j])){
          begin = as.vector(unlist(model$metNames[met_index[j]], use.names = FALSE))
        } else {
          begin = path[j]
        }
        
        rxns=getReactionsRecon(index = tmp$mets, model)
        mets_sub = cbind(mets_sub, 
                         "resource"=rep("Recon",dim(mets_sub)[1]), 
                         "rxn_formula"=rxns[,"rxns"], 
                         "path"=paste(begin," =>(", as.vector(unlist(model$rxnNames[rxn_index[k]], use.names = FALSE)),")=>",mets_sub[,"met_long"])
                         )
        mets = rbind(mets, mets_sub)
      }
    } 
  }
  cat("\n")
  if (is.null(mets)){
    return(NULL)
  } else {
    return(list("mets"=mets))
  }
}
