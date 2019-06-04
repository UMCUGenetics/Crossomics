getPreviousNextRecon2.2 <- function(model, met_index, produced, max_rxns, recon2chebi){
#     met_index=index_consumed
#     produced=FALSE

  mets=NULL
  
  for(j in 1:length(met_index)){
    
    message(as.vector(unlist(model@met_name[met_index[j]], use.names = FALSE)))
    
    #     message(met_index)
    #     message(produced)
    
    met=as.vector(unlist(model@met_id[met_index[j]], use.names = FALSE)) 
    
    # leave this message it avoids a bug!!!!!!!!!!!!!!!!
    message(dim(model@S))
    
    # find in which reaction consumed
    if (produced){
      
      rxn_index = try({which(model@S[as.numeric(met_index[j]),]<0)},silent = TRUE)
      if (class(rxn_index) == "try-error") {
        message(rxn_index)
        message(as.numeric(met_index[j]))
        message(j)
        next
      }  
      
      label = "consumed in: "
    } else {
      # find in which reaction produced  
      
      rxn_index = try({which(model@S[as.numeric(met_index[j]),]>0)},silent = TRUE)
      if (class(rxn_index) == "try-error") {
        message(rxn_index)
        message(as.numeric(met_index[j]))
        message(j)
        next
      }  
      
      label = "produced in: "
    }  
    
    
    # if (length(rxn_index)>0){  
    if (max_rxns>=length(rxn_index)&length(rxn_index)>0){  
      for(k in 1:length(rxn_index)){
        
        #message(paste("k", k, sep=" "))
        
        message(paste(label, as.vector(unlist(model@react_name[rxn_index[k]], use.names = FALSE)), sep=""))
        
        # Get mets of reaction via S not through genes!!!!
        if (produced) {
          mets = rbind(mets, getMetsRecon2.2(rxn_index[k], level=1, met, model, recon2chebi)$metabolites)
          # print(result2$metabolites)
        } else {
          mets = rbind(mets, getMetsRecon2.2(rxn_index[k], level=-1, met, model, recon2chebi)$metabolites)
          # print(result2$metabolites)
        }
      }
    }  
  }
  
  if (is.null(mets)){
    return(NULL)
  } else {
    return(list("mets"=mets))
  }
}
