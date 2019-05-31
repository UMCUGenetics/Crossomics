getPreviousNext <- function(model, met_index, max_rxns, step, path){
# met_index=index_primary_removed
# step=1
# path=rep(NA,length(index_primary_removed))
  
  # leave this message it avoids a bug!!!!!!!!!!!!!!!!
  # message(dim(model$S))
  
  # as.vector(unlist(model$metNames[met_index], use.names = FALSE))
  
  usedIn <- function(model, met_index){
    # find in which reaction used
    rxn_index = try({which(model$S[as.numeric(met_index),]!=0)},silent = TRUE)
    if (class(rxn_index) == "try-error") {
      rxn_index = -1  
    }
    return(rxn_index)
  }
  
  mets=NULL
  
  for(j in 1:length(met_index)){
    # message("====================================")
    # message(paste("===>", as.vector(unlist(model$metNames[met_index[j]], use.names = FALSE))))
    # message("====================================")
    #     message(met_index)
    #     message(produced)

    met=as.vector(unlist(model$mets[met_index[j]], use.names = FALSE)) 

    # # find in which reaction used
    # rxn_index = try({which(model$S[as.numeric(met_index[j]),]!=0)},silent = TRUE)
    # if (class(rxn_index) == "try-error") {
    #   # message(rxn_index)
    #   # message(as.numeric(met_index[j]))
    #   # message(j)
    #   next
    # }  
    
    rxn_index = usedIn(model, met_index[j])
    if (rxn_index[1]==-1) next
    
    label = "used in: "


    for(k in 1:length(rxn_index)){

      # message(paste(label, as.vector(unlist(model$rxnNames[rxn_index[k]], use.names = FALSE)), sep=""))

      level=step

      tmp = getMetsRecon2V2(rxn_index[k], level, met, model)
      mets_sub = tmp$metabolites

      # get Recon Metabolite ID
      recon_met_ids = which(model$S[,as.numeric(rxn_index[k])]!=0)
      # model$metNames[recon_met_ids[l]]
      rownames(mets_sub)=recon_met_ids 
      
      for (l in 1:length(recon_met_ids)) {
        
        # only add if met is used in less then max_rxns
        if (max_rxns>=length(usedIn(model,recon_met_ids[l]))){
          
           if (is.na(path[j])){
             begin = as.vector(unlist(model$metNames[met_index[j]], use.names = FALSE))
           } else {
             begin = path[j]
           }
          
           rxns=getReactionsRecon(recon_met_ids[l], model)
           row = mets_sub[toString(recon_met_ids[l]),,drop=FALSE]
           
           row = cbind(row, "resource"=rep("Recon",1), "rxn_formula"=rxns[,"rxns"], "path"=paste(begin,
                                                                    "=>(", as.vector(unlist(model$rxnNames[rxn_index[k]], use.names = FALSE)),")=>",row[,"met_long"]))
           mets = rbind(mets, row)
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
