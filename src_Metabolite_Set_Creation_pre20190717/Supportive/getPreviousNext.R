getPreviousNext <- function(model, met_index, max_rxns, step, path){
# met_index=index_primary_removed
# step=1
# path=rep(NA,length(index_primary_removed))
  
  # leave this message it avoids a bug!!!!!!!!!!!!!!!!
  # message(dim(model$S))
  
  # as.vector(unlist(model$metNames[met_index], use.names = FALSE))
  
  mets=NULL
  
  for(j in 1:length(met_index)){
    # message("====================================")
    # message(paste("===>", as.vector(unlist(model$metNames[met_index[j]], use.names = FALSE))))
    # message("====================================")
    #     message(met_index)
    #     message(produced)

    met=as.vector(unlist(model$mets[met_index[j]], use.names = FALSE)) 

    # find in which reaction used
    rxn_index = try({which(model$S[as.numeric(met_index[j]),]!=0)},silent = TRUE)
    if (class(rxn_index) == "try-error") {
      # message(rxn_index)
      # message(as.numeric(met_index[j]))
      # message(j)
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

        # message("################################")
        # message(paste(mets_sub[,"met_long"],collapse = "\n"))
        # message("################################")

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
    
    #else {
    #   message(paste("NO", met))
    # } 
  }
  
  if (is.null(mets)){
    return(NULL)
  } else {
    return(list("mets"=mets))
  }
}
