# Function to extend metabolite sets, supplying the step size, and max reactions
extend_metabolites <- function(model, max_rxns, index, mylist, step){
  
  extended_mets <- NULL
  
  for(index_rxn in 1:length(max_rxns)){
    if(is.null(extended_mets)){
      allindices <- index[[index_rxn]]
      if(length(allindices) != 0) {
        result_tmp <- getPreviousNext(model,
                                      met_index = allindices,
                                      step = step,
                                      path = rep(NA,length(allindices))
        )
        # Record which metabolites could already be extended when alowing less reactions (and skip them later)
        extended_mets <- result_tmp$extended_mets
        result_tmp <- result_tmp$mets
      } else {
        result_tmp <- data.frame(
          "rxn_id" = NULL, "step" = NULL,"met_in" = NULL,"left_right" = NULL,"met_short" = NULL,"met_long" = NULL,"hmdb" = NULL,
          "kegg" = NULL,"chebi"  = NULL,"pubchem"  = NULL,"rxn_name" = NULL,"rxn" = NULL,"resource" = NULL,"rxn_formula" = NULL,"path" = NULL)
      }
        
    } else {
      someindices <- index[[index_rxn]][!index[[index_rxn]] %in% extended_mets]
      previous_result <- result_tmp
      if(length(someindices) != 0) {
        result_tmp <- getPreviousNext(model,
                                      met_index = someindices,
                                      step = step,
                                      path = rep(NA,length(someindices))
        )
        # Record which metabolites could already be extended when alowing less reactions (and skip them later)
        extended_mets <- c(extended_mets, result_tmp$extended_mets)
        result_tmp <- rbind(previous_result, result_tmp$mets)
      } else {
        result_tmp <- previous_result
      }
        
    }

  #   if(length(index[[index_rxn]]) != 0) {
  #     result_tmp <- getPreviousNext(model,
  #                                   met_index = index[[index_rxn]],
  #                                   step = step,
  #                                   path = rep(NA,length(index[[index_rxn]]))
  #     )
  #     # Record which metabolites could already be extended when alowing less reactions (and skip them later)
  #     extended_mets <- c(extended_mets, result_tmp$extended_mets)
  #     result_tmp <- result_tmp$mets
  #   } else {
  #     result_tmp <- data.frame(
  #       "rxn_id" = NULL, "step" = NULL,"met_in" = NULL,"left_right" = NULL,"met_short" = NULL,"met_long" = NULL,"hmdb" = NULL,
  #       "kegg" = NULL,"chebi"  = NULL,"pubchem"  = NULL,"rxn_name" = NULL,"rxn" = NULL,"resource" = NULL,"rxn_formula" = NULL,"path" = NULL)
  #   }
  # }
    
    if (index_rxn == 1) {
      mylist[[index_rxn]] <- result_tmp[,, drop = FALSE]
    } else {
      mylist[[index_rxn]] <- result_tmp
    }
  }
  return(mylist)
}

  