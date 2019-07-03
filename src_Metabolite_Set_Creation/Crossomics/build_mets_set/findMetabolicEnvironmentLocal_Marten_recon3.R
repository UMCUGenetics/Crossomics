findMetabolicEnvironmentLocal <- function(gene_in, model, recon2chebi, src, rval){
# gene_in=hgnc

  removeExchange <- function(mets){
    #mets=result_mets_0
    
    # remove: "Exchange/demand reaction", "Transport, extracellular",
    index=grep("Exchange/demand reaction", mets[,"rxn"], fixed=TRUE)
    if (length(index)>0) mets = mets[-index,]

    index=grep("Transport, extracellular", mets[,"rxn"], fixed=TRUE)
    if (length(index)>0) mets = mets[-index,]

    # genExcelFileShort(as.data.frame(mets), paste("./results/mets.min.", HGNC,".xls",sep=""))

    # A = mets[,"met_short"]
    A = mets[,"met_long"]
    unic = !duplicated(A)  ## logical vector of unique values
    index = seq_along(A)[unic]  ## indices
    mets = mets[index, , drop = FALSE]

    # genExcelFileShort(as.data.frame(mets), paste("./results/mets.min.unic.", HGNC,".xls",sep=""))

    return(mets)
  }
  
  # message(gene_in)

  # if (reconVersion==2.2){
  #   S=model@S
  # } else if (reconVersion==2.0){
  #   S=model$S
  # }

  S=model$S
  max_rxns=10

  rvalMets=NULL
  rvalRxns=NULL
  index_consumed=NULL
  index_produced=NULL

  if (!is.null(rval)) {
    consumed = rval$left
    produced = rval$right
    rxns = rval$rxns
    
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
                                   # "InChI_key"=consumed[,"InchI_key"],
                                   "rxn_name"=rep(NA,amount),
                                   "rxn"=rep(NA,amount),
                                   "resource"=consumed[,"Data Source"])
      tmp=NULL
      for (i in 1:dim(result_mets_consumed)[1]){
        source = rxns[grep(result_mets_consumed[i,"resource"], names(rxns), fixed=TRUE)]
        tmp=c(tmp,
        paste(source[grep(result_mets_consumed[i,"met_short"], source, fixed = TRUE)],collapse = ";"))
      }
      
      result_mets_consumed = cbind(result_mets_consumed,"rxn_formula"=tmp)
      result_mets_consumed = cbind(result_mets_consumed,"path"=rep(NA,amount))
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
                                   # "InChI_key"=produced[,"InchI_key"],
                                   "rxn_name"=rep(NA,amount),
                                   "rxn"=rep(NA,amount),
                                   "resource"=produced[,"Data Source"])
      tmp=NULL
      for (i in 1:dim(result_mets_produced)[1]){
        source = rxns[grep(result_mets_produced[i,"resource"], names(rxns), fixed=TRUE)]
        tmp=c(tmp,
              paste(source[grep(result_mets_produced[i,"met_short"], source, fixed = TRUE)],collapse = ";"))
      }
      
      result_mets_produced = cbind(result_mets_produced,"rxn_formula"=tmp)
      result_mets_produced = cbind(result_mets_produced,"path"=rep(NA,amount))
    }      

    # Redundant as another ordering is done after the following block of code
    # result_mets_consumed=result_mets_consumed[order(result_mets_consumed[,"rxn_formula"]),]
    # result_mets_produced=result_mets_produced[order(result_mets_produced[,"rxn_formula"]),]
    
    if (length(consumed)>0 && length(produced)>0) {
      result_mets <- rbind(result_mets_consumed, result_mets_produced)
    } else if (length(consumed)==0 && length(produced)>0) {
      result_mets <- result_mets_produced
    } else if (length(consumed)>0 && length(produced)==0) {
      result_mets <- result_mets_consumed
    } else {
      result_mets <- NULL
    }
    
    result_mets <- result_mets[order(result_mets[,"rxn_formula"]),, drop=FALSE]
    index = getMets2ReconID(mets = result_mets, model)

    if (length(index)>0){
      # replace empty entries by NA
      for (j in 1:length(index)){
        if (length(convert(model$metHMDB[index[j]]))==0) model$metHMDB[index[j]][[1]] <- list(matrix(NA))
      }
      
      amount = length(convert(model$mets[index]))
      rxn=NA
      if (length(index)>0) rxn <- getReactionsRecon(index, model)
      
      result_mets = rbind(result_mets,
                                   cbind("rxn_id"=rep(toString(gene_in),amount),
                                         "step"=rep(0,amount),
                                         "met_in"=rep(NA,amount),
                                         "left_right"=rep(NA,amount),
                                         "met_short"=convert(model$mets[index]),
                                         "met_long"=convert(model$metNames[index]),
                                         "hmdb"=convert(model$metHMDB[index]),
                                         "kegg"=convert(model$metKEGGID[index]),
                                         "chebi"=convert(model$metCHEBIID[index]),
                                         "pubchem"=convert(model$metPubChemID[index]),
                                         # "InChI_key"=rep(NA,amount),
                                         "rxn_name"=rep(NA,amount),
                                         "rxn"=rep(NA,amount),
                                         "resource"=rep("Recon",amount),
                                         "rxn_formula"=rxn[,"rxns"],
                                         "path"=rep(NA,amount)))
    }
    
    ################################################################################################
    
    # unique(as.vector(unlist(model$metNames[index_primary_removed])))

    # step -1,1 ###################################################################################
    index_primary = unique(c(index))
    index_primary_removed = removeMetsFromSet(index_primary,model)
    
# unique(unlist(model$metNames[index_primary]))
# unique(unlist(model$metNames[index_primary_removed]))
    
    result_mets_1 = NULL
    if (length(index_primary_removed)>0) {
      result_mets_1 = getPreviousNext(model, index_primary_removed, max_rxns, step=1, rep(NA,length(index_primary_removed)))$mets
      # B = getPreviousNext(model, index_primary, produced=TRUE, max_rxns, step=1)
    }
    if (is.null(result_mets_1)) result_mets_1 = result_mets  
    
# unique(result_mets_1[,"met_long"])
    
    # step -2,2 ###################################################################################
    # tmp=result_mets_1
    # index_1 = getMets2ReconID(tmp, model)
    index_1 = which(as.vector(unlist(model$mets)) %in% as.vector(result_mets_1[,"met_short"]))
    # checked => oke
    # cbind(sort(unique(as.vector(unlist(model$mets[index_1])))),
    # sort(unique(as.vector(result_mets_1[,"met_short"]))))
    index_1_removed = removeMetsFromSet(index_1,model)
    
# unique(unlist(model$metNames[index_1]))
# unique(unlist(model$metNames[index_1_removed]))

    met_short_recon = as.vector(unlist(model$mets[index_1_removed]))
    index_path = NULL
    for (i in 1:length(met_short_recon)){
      tmp = as.numeric(which(result_mets_1[,"met_short"] == met_short_recon[i]))
      if (length(tmp)>1) tmp=tmp[1] 
      index_path = c(index_path,tmp)
    }
    # checked => oke
    # cbind(as.vector(unlist(model$mets[index_1_removed])),result_mets_1[index_path,c("met_short","path")])

    result_mets_2 = NULL
    if (length(index_1_removed)>0) {
      result_mets_2 = getPreviousNext(model, index_1_removed, max_rxns, step=2, as.vector(result_mets_1[index_path,"path"]))$mets
      # A = getPreviousNext(model, index_1, produced=TRUE, max_rxns, step=2)
    }
    if (is.null(result_mets_2)) result_mets_2 = result_mets_1

# unique(result_mets_2[,"met_long"])

    # step -3,3 ###################################################################################
    # tmp=result_mets_2
    # index_2 = getMets2ReconID(tmp, model)
    # index_2_removed = removeMetsFromSet(index_2,model)
    index_2 = which(as.vector(unlist(model$mets)) %in% as.vector(result_mets_2[,"met_short"]))
    index_2_removed = removeMetsFromSet(index_2,model)
    
    met_short_recon = as.vector(unlist(model$mets[index_2_removed]))
    index_path = NULL
    for (i in 1:length(met_short_recon)){
      tmp = as.numeric(which(result_mets_2[,"met_short"] == met_short_recon[i]))
      if (length(tmp)>1) tmp=tmp[1] 
      index_path = c(index_path,tmp)
    }

    result_mets_3 = NULL
    if (length(index_2_removed)>0) {
      result_mets_3 = getPreviousNext(model, index_2_removed, max_rxns, step=3, as.vector(result_mets_2[index_path,"path"]))$mets
      # A = getPreviousNext(model, index_1, produced=TRUE, max_rxns, step=2)
    }
    if (is.null(result_mets_3)) result_mets_3 = result_mets_2

    # step -4,4 ###################################################################################
    index_3 = which(as.vector(unlist(model$mets)) %in% as.vector(result_mets_3[,"met_short"]))
    index_3_removed = removeMetsFromSet(index_3,model)
    
    met_short_recon = as.vector(unlist(model$mets[index_3_removed]))
    index_path = NULL
    for (i in 1:length(met_short_recon)){
      tmp = as.numeric(which(result_mets_3[,"met_short"] == met_short_recon[i]))
      if (length(tmp)>1) tmp=tmp[1] 
      index_path = c(index_path,tmp)
    }
    
    result_mets_4 = NULL
    if (length(index_3_removed)>0) {
      result_mets_4 = getPreviousNext(model, index_3_removed, max_rxns, step=4, as.vector(result_mets_3[index_path,"path"]))$mets
      # A = getPreviousNext(model, index_1, produced=TRUE, max_rxns, step=2)
    }
    if (is.null(result_mets_4)) result_mets_4 = result_mets_3
    # result_mets=result_mets[,-c(which(colnames(result_mets)=="InChI_key"))]

    result_mets_4 = rbind(result_mets,result_mets_1,result_mets_2,result_mets_3,result_mets_4)
    result_mets_3 = rbind(result_mets,result_mets_1,result_mets_2,result_mets_3)
    result_mets_2 = rbind(result_mets,result_mets_1,result_mets_2)
    result_mets_1 = rbind(result_mets,result_mets_1)
    result_mets_0 = result_mets
    
    if (dim(result_mets_4)[1]>0) result_mets_4 = removeExchange(result_mets_4)
    if (dim(result_mets_3)[1]>0) result_mets_3 = removeExchange(result_mets_3)
    if (dim(result_mets_2)[1]>0) result_mets_2 = removeExchange(result_mets_2)
    if (dim(result_mets_1)[1]>0) result_mets_1 = removeExchange(result_mets_1)
    if (dim(result_mets_0)[1]>0) result_mets_0 = removeExchange(result_mets_0)

    save(result_mets_4, file=paste(src, "../../Results/mss_4",paste(gene_in, "RData", sep="."), sep="/"))
    save(result_mets_3, file=paste(src, "../../Results/mss_3",paste(gene_in, "RData", sep="."), sep="/"))
    save(result_mets_2, file=paste(src, "../../Results/mss_2",paste(gene_in, "RData", sep="."), sep="/"))
    save(result_mets_1, file=paste(src, "../../Results/mss_1",paste(gene_in, "RData", sep="."), sep="/"))
    save(result_mets_0, file=paste(src, "../../Results/mss_0",paste(gene_in, "RData", sep="."), sep="/"))
    
    # load(paste(getwd(),src, "../../Results/metabolite_sets_step_0,1,2_filter_1.0/mss_0",paste(hgnc, "RData", sep="."), sep="/"))
    # library(XLConnect)
    # genExcelFileShort(result_mets_3, "RNF123_metabolite_set_step_3.xls")
  }  
}