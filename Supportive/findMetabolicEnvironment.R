findMetabolicEnvironment <- function(gene_in, model, HGNC, onlyPrimary=FALSE){
#     gene_in=entrezgene
#     HGNC=hgnc
#     onlyPrimary=FALSE
  
#     gene_in=79048
#     HGNC="SECISBP2"
  
  message(gene_in)
  message(HGNC)
  
  #S=model@S
  S=model$S
  max_rxns=4
  
  convert <- function(a) {as.vector(unlist(lapply(a, function(x) {
    if (identical(unlist(x),character(0))){
      return(NA)
    } else {
      return(x)
    }
  })))}
  
  rvalMets=NULL
  rvalRxns=NULL
  index_consumed=NULL
  index_produced=NULL
  
  # result = getMetsRecon2(59272, -1, "ser_L[c]")
  #   consumed = getMetsKEGG(gene_in,consumedMets=TRUE)
  #   produced = getMetsKEGG(gene_in,consumedMets=FALSE)
  
  # Get primary reaction
  rval = getMetsPathwayCommons(gene_in)
  if (is.null(rval)) {
    return(NULL)  
  } else {
    
    consumed = rval$left
    produced = rval$right
    
    if (length(consumed)>0) {
      # Do not convert to ReconID in primary reaction!
      amount = dim(consumed)[1]
      result_mets_consumed = cbind("rxn_id"=rep(toString(gene_in),amount),
                                   "step"=rep(0,amount),
                                   "met_in"=rep(NA,amount),
                                   "left_right"=rep("left",amount),
                                   "met_short"=consumed[,"Compound"],
                                   "met_long"=consumed[,"Compound"],
                                   "hmdb"=rep(NA,amount),
                                   "kegg"=consumed[,"KEGG_id"],
                                   "chebi"=consumed[,"CheBI_id"],
                                   "pubchem"=consumed[,"PubChem_id"],
                                   "rxn_name"=rep(NA,amount),
                                   "rxn"=rep(NA,amount))
      
      if (!onlyPrimary) {
        
        index_consumed = getMets2ReconID(consumed, model)
        
        # replace empty entries by NA
        for (j in 1:length(index_consumed)){
          if (length(convert(model$metHMDB[index_consumed[j]]))==0) model$metHMDB[index_consumed[j]][[1]]=list(matrix(NA))
        }
        
        # check mets in
        amount = length(convert(model$mets[index_consumed]))
        result_mets_consumed = rbind(result_mets_consumed,
                                     cbind("rxn_id"=rep(toString(gene_in),amount),
                                           "step"=rep(0,amount),
                                           "met_in"=rep(NA,amount),
                                           "left_right"=rep("left",amount),
                                           "met_short"=convert(model$mets[index_consumed]),
                                           "met_long"=convert(model$metNames[index_consumed]),
                                           "hmdb"=convert(model$metHMDB[index_consumed]),
                                           "kegg"=convert(model$metKeggID[index_consumed]),
                                           "chebi"=convert(model$metCHEBIID[index_consumed]),
                                           "pubchem"=convert(model$metPubChemID[index_consumed]),
                                           "rxn_name"=rep(NA,amount),
                                           "rxn"=rep(NA,amount)))
      }  
    }
    
    if (length(produced)>0) {
      # Do not convert to ReconID if only primary reaction used!
      amount = dim(produced)[1]
      result_mets_produced = cbind("rxn_id"=rep(toString(gene_in),amount),
                                   "step"=rep(0,amount),
                                   "met_in"=rep(NA,amount),
                                   "left_right"=rep("right",amount),
                                   "met_short"=produced[,"Compound"],
                                   "met_long"=produced[,"Compound"],
                                   "hmdb"=rep(NA,amount),
                                   "kegg"=produced[,"KEGG_id"],
                                   "chebi"=produced[,"CheBI_id"],
                                   "pubchem"=produced[,"PubChem_id"],
                                   "rxn_name"=rep(NA,amount),
                                   "rxn"=rep(NA,amount))
      
      if (!onlyPrimary) {
        
        index_produced = getMets2ReconID(produced, model)
        
        # replace empty entries by NA
        for (j in 1:length(index_produced)){
          if (length(convert(model$metHMDB[index_produced[j]]))==0) model$metHMDB[index_produced[j]][[1]]=list(matrix(NA))
        }
        
        amount = length(convert(model$mets[index_produced]))
        result_mets_produced = rbind(result_mets_consumed,  
                                     cbind("rxn_id"=rep(toString(gene_in),amount),
                                           "step"=rep(0,amount),
                                           "met_in"=rep(NA,amount),
                                           "left_right"=rep("right",amount),
                                           "met_short"=convert(model$mets[index_produced]),
                                           "met_long"=convert(model$metNames[index_produced]),
                                           "hmdb"= convert(model$metHMDB[index_produced]),
                                           "kegg"= convert(model$metKeggID[index_produced]),
                                           "chebi"= convert(model$metCHEBIID[index_produced]),
                                           "pubchem"= convert(model$metPubChemID[index_produced]),
                                           "rxn_name"=rep(NA,amount),
                                           "rxn"=rep(NA,amount)))
      } 
    }
    
    if (length(consumed)>0 && length(produced)>0) {
      result_mets = rbind(result_mets_consumed, result_mets_produced)
    } else if (length(consumed)==0 && length(produced)>0) {
      result_mets = result_mets_produced
    } else if (length(consumed)>0 && length(produced)==0) {
      result_mets = result_mets_consumed
    } else {
      result_mets = NULL
    }
    
    if (!onlyPrimary & ((length(index_consumed)==0) | (length(index_produced)==0))){
      result_mets = NULL
    }
    
    # rvalMets = result_mets 
    
    result_mets_consumed=NULL
    result_mets_produced=NULL
    
    if (!is.null(result_mets)) {
      
      #       rval = removeMetsFromSet(result_mets, model, index_consumed, index_produced)
      #       
      #       result_mets = rval$result_mets
      #       index_produced = rval$index_produced
      #       index_consumed = rval$index_consumed
      
      if (!onlyPrimary) {
        
        if (length(index_consumed)>0) {
          # step -1 ###################################################################################
          A = getPreviousNext(model, index_consumed, produced=FALSE, max_rxns)
          B = getPreviousNext(model, index_consumed, produced=TRUE, max_rxns)
          result_mets_consumed = rbind(A$mets,B$mets)
          
          # result_mets_consumed = removeMetsFromSet(result_mets_consumed, model)$result_mets
        }
        
        if (length(index_produced)>0) {
          # step 1 ####################################################################################
          A = getPreviousNext(model, index_produced, produced=TRUE, max_rxns)
          B = getPreviousNext(model, index_produced, produced=FALSE, max_rxns)
          result_mets_produced = rbind(A$mets,B$mets)
          
          # result_mets_produced = removeMetsFromSet(result_mets_produced, model)$result_mets
        }
        
        result_from_consumed_plus_produced = rbind(result_mets,result_mets_consumed,result_mets_produced)
        
        if (dim(result_from_consumed_plus_produced)[1]>0){  
          # dir.create("./results", showWarnings = FALSE)
          # library(XLConnect)
          # genExcelFileShort(result_from_consumed_plus_produced, paste("./results/mets",HGNC,"xls", sep="."))
          
          
          # remove: "Exchange/demand reaction", "Transport, extracellular",
          index=grep("Exchange/demand reaction", result_from_consumed_plus_produced[,"rxn"], fixed=TRUE)
          if (length(index)>0) result_from_consumed_plus_produced = result_from_consumed_plus_produced[-index,]
          index=grep("Transport, extracellular", result_from_consumed_plus_produced[,"rxn"], fixed=TRUE)
          if (length(index)>0) result_from_consumed_plus_produced = result_from_consumed_plus_produced[-index,]
          
          # genExcelFileShort(as.data.frame(result_from_consumed_plus_produced), paste("./results/mets.min.", HGNC,".xls",sep=""))
          
          #           A = result_from_consumed_plus_produced[,"met_short"]
          A = result_from_consumed_plus_produced[,"met_long"]
          unic = !duplicated(A)  ## logical vector of unique values
          index = seq_along(A)[unic]  ## indices
          result_from_consumed_plus_produced = result_from_consumed_plus_produced[index,]
          
          # genExcelFileShort(as.data.frame(result_from_consumed_plus_produced), paste("./results/mets.min.unic.", HGNC,".xls",sep=""))
          
          rvalMets=result_from_consumed_plus_produced
          
        }
      } else {
        
        # A = result_mets[,"met_short"]
        A = result_mets[,"met_long"]
        unic = !duplicated(A)  ## logical vector of unique values
        index = seq_along(A)[unic]  ## indices
        result_mets = result_mets[index,]
        
        rvalMets=result_mets
      }
    }
    
    return(list("mets"=rvalMets))
  }
}
