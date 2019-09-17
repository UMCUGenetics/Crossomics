findMetabolicEnvironmentLocal <- function(gene_in, model, shortmodel, recon2chebi, src, outdir, rval, max_rxns, steps){

###########################################################################
# Functions ---------------------------------------------------------------
###########################################################################
  
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
  
  # save_files <- function(r0, r1, r2, r3, r4, r5){
  save_files <- function(r0, r1 = NULL){
    cat("Finishing up\n\n")
    
    # result_mets_0 <- result_mets
    # result_mets_1 <- rbind(result_mets_0,result_mets_1)
    # result_mets_2 <- rbind(result_mets_1,result_mets_2)
    # result_mets_3 <- rbind(result_mets_2,result_mets_3)
    # result_mets_4 <- rbind(result_mets_3,result_mets_4)
    # result_mets_5 <- rbind(result_mets_4,result_mets_5)
    result_mets_final <- vector("list", length(1:steps))
    result_mets_0 <- r0
    # For when all metabolites have been filtered out before the first extention step
    if(is.null(r1)){
      for(i in 1:length(max_rxns)){
        for (step in 1:steps){
          result_mets_final[[step]][[i]] <- result_mets_0
        }
        for(step in steps:1){
          if (nrow(result_mets_final[[step]][[i]])>0) result_mets_final[[step]][[i]] <- removeExchange(result_mets_final[[step]][[i]])
        }
      }
    } else {
      for(i in 1:length(max_rxns)){
        for (step in 1:steps){
          if(step == 1) {
            result_mets_final[[step]][[i]] <- rbind(result_mets_0, r1[[step]][[i]])
          } else {
            result_mets_final[[step]][[i]] <- rbind(result_mets_final[[step-1]][[i]], r1[[step]][[i]])
          }
        }
        for(step in steps:1){
          if (nrow(result_mets_final[[step]][[i]])>0) result_mets_final[[step]][[i]] <- removeExchange(result_mets_final[[step]][[i]])
        }
      }
    }
    
    if (nrow(result_mets_0)>0) result_mets_0 <- removeExchange(result_mets_0)
    # for(i in 1:length(max_rxns)){
    #   result_mets_1[[i]] <- rbind(result_mets_0,result_mets_1[[i]])
    #   result_mets_2[[i]] <- rbind(result_mets_1[[i]],result_mets_2[[i]])
    #   result_mets_3[[i]] <- rbind(result_mets_2[[i]],result_mets_3[[i]])
    #   result_mets_4[[i]] <- rbind(result_mets_3[[i]],result_mets_4[[i]])
    #   result_mets_5[[i]] <- rbind(result_mets_4[[i]],result_mets_5[[i]])
    #   if (nrow(result_mets_5[[i]])>0) result_mets_5[[i]] <- removeExchange(result_mets_5[[i]])
    #   if (nrow(result_mets_4[[i]])>0) result_mets_4[[i]] <- removeExchange(result_mets_4[[i]])
    #   if (nrow(result_mets_3[[i]])>0) result_mets_3[[i]] <- removeExchange(result_mets_3[[i]])
    #   if (nrow(result_mets_2[[i]])>0) result_mets_2[[i]] <- removeExchange(result_mets_2[[i]])
    #   if (nrow(result_mets_1[[i]])>0) result_mets_1[[i]] <- removeExchange(result_mets_1[[i]])
    # }
    # if (nrow(result_mets_0)>0) result_mets_0 <- removeExchange(result_mets_0)
    
    # if (nrow(result_mets_5)>0) result_mets_5 <- removeExchange(result_mets_5)
    # if (nrow(result_mets_4)>0) result_mets_4 <- removeExchange(result_mets_4)
    # if (nrow(result_mets_3)>0) result_mets_3 <- removeExchange(result_mets_3)
    # if (nrow(result_mets_2)>0) result_mets_2 <- removeExchange(result_mets_2)
    # if (nrow(result_mets_1)>0) result_mets_1 <- removeExchange(result_mets_1)
    # if (nrow(result_mets_0)>0) result_mets_0 <- removeExchange(result_mets_0)
    


    for(i in 1:length(max_rxns)){
      for(step in 0:steps){
        if(step == 0 ) {
          # save(result_mets_0, file=paste0(outdir,"/maxrxn",max_rxns[i],"/mss_",step,"/",gene_in,".RData"))
          saveRDS(result_mets_0, file=paste0(outdir,"/maxrxn",max_rxns[i],"/mss_",step,"/",gene_in,".RDS"))
        } else {
          save_this_file <- result_mets_final[[step]][[i]]
          # save(save_this_file, file=paste0(outdir,"/maxrxn",max_rxns[i],"/mss_",step,"/",gene_in,".RData"))
          saveRDS(save_this_file, file=paste0(outdir,"/maxrxn",max_rxns[i],"/mss_",step,"/",gene_in,".RDS"))
        }
      }
      # save(result_mets_5[[i]], file=paste0(outdir,"/maxrxn",max_rxns[i],"/mss_5/",gene_in,".RData"))
      # save(result_mets_4[[i]], file=paste0(outdir,"/maxrxn",max_rxns[i],"/mss_4/",gene_in,".RData"))
      # save(result_mets_3[[i]], file=paste0(outdir,"/maxrxn",max_rxns[i],"/mss_3/",gene_in,".RData"))
      # save(result_mets_2[[i]], file=paste0(outdir,"/maxrxn",max_rxns[i],"/mss_2/",gene_in,".RData"))
      # save(result_mets_1[[i]], file=paste0(outdir,"/maxrxn",max_rxns[i],"/mss_1/",gene_in,".RData"))
      # save(result_mets_0[[i]], file=paste0(outdir,"/maxrxn",max_rxns[i],"/mss_0/",gene_in,".RData"))
      
    }
    # save(result_mets_5, file=paste(outdir,"mss_5",paste0(gene_in,".RData"), sep="/"))
    # save(result_mets_4, file=paste(outdir,"mss_4",paste0(gene_in,".RData"), sep="/"))
    # save(result_mets_3, file=paste(outdir,"mss_3",paste0(gene_in,".RData"), sep="/"))
    # save(result_mets_2, file=paste(outdir,"mss_2",paste0(gene_in,".RData"), sep="/"))
    # save(result_mets_1, file=paste(outdir,"mss_1",paste0(gene_in,".RData"), sep="/"))
    # save(result_mets_0, file=paste(outdir,"mss_0",paste0(gene_in,".RData"), sep="/"))
    
  }
  
  
  
  
###########################################################################
# Code --------------------------------------------------------------------
###########################################################################  
  
  rvalMets=NULL
  rvalRxns=NULL
  index_consumed=NULL
  index_produced=NULL
  
  consumed = rval$left
  produced = rval$right
  rxns = rval$rxns
  
  if (length(consumed)>0) {
    
    # Do not convert to ReconID in primary reaction!
    amount <- nrow(consumed)
    result_mets_consumed <- cbind("rxn_id"=rep(toString(gene_in),amount),
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
    tmp <- NULL
    for (i in 1:nrow(result_mets_consumed)){
      source <- rxns[grep(result_mets_consumed[i,"resource"], names(rxns), fixed=TRUE)]
      tmp <- c(tmp, paste(source[grep(result_mets_consumed[i,"met_short"], source, fixed = TRUE)],collapse = ";"))
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
      tmp=c(tmp, paste(source[grep(result_mets_produced[i,"met_short"], source, fixed = TRUE)],collapse = ";"))
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
  index <- getMets2ReconID(mets = result_mets, model)
  
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
  

  ###########################################################################
  # Prepare result_mets lists -----------------------------------------------
  ###########################################################################
  # Trying to make a list of result_mets, with as many entries as I want max_rxns
  list_names <- unlist(lapply(max_rxns, function(x) paste0("maxrxn",x)))
  mylist <- vector("list", length(list_names))
  names(mylist) <- list_names
  
  result_mets_all <- vector("list", length(1:steps))

  # result_mets_1 <- mylist
  # result_mets_2 <- mylist
  # result_mets_3 <- mylist
  # result_mets_4 <- mylist
  # result_mets_5 <- mylist
  
  
  
  ###########################################################################
  # Step = 1 ----------------------------------------------------------------
  ###########################################################################
  cat("step 1.. \n")
  
  index_primary <- unique(c(index))
  index_primary_removed <- removeMetsFromSet(index = index_primary, model)
  
  if(length(index_primary_removed) == 0){
    save_files(result_mets, NULL)
    return()
  }
  
  # if(length(index_primary_removed) == 0) {
  #   save_files(result_mets, result_mets_1, result_mets_2, result_mets_3, result_mets_4)
  #   return()
  # }
  
  index_primary_removed2 <- lapply(1:length(max_rxns), function(x) {
    min_rxn <- ifelse(x==1, 0, max_rxns[x-1])
    max_rxn <- max_rxns[x]
    filter_indices(model, shortmodel, index_primary_removed, max_rxns = max_rxn, min_rxns = min_rxn)
  })
    
  # result_mets_1 <- extend_metabolites(model, max_rxns, index = index_primary_removed2, mylist = mylist, step = 1)
  result_mets_all[[1]] <- extend_metabolites(model, max_rxns, index = index_primary_removed2, mylist = mylist, step = 1)

  # if(length(index_primary_removed) == 0) {
  #   save_files(result_mets, result_mets_1, result_mets_2, result_mets_3, result_mets_4)
  #   return()
  # }
    
  
  
  
  ###########################################################################
  # Step >= 2 ---------------------------------------------------------------
  ###########################################################################
  cat("step 2.. \n")
  
  # results <- further_extend_metabolites(model, shortmodel, mylist, result_mets_1, max_rxns, index_primary_removed2, step = 2)
  # result_mets_2 <- results[["result_mets"]]
  results <- further_extend_metabolites(model, shortmodel, mylist, result_mets_all[[1]], max_rxns, index_primary_removed2, step = 2)
  result_mets_all[[2]] <- results[["result_mets"]]
  index_1_removed2 <- mapply(c, index_primary_removed2, results[["index_removed"]], SIMPLIFY=FALSE)
  rm(results)
    
  #step3
  cat("step 3.. \n")
  # results <- further_extend_metabolites(model, shortmodel, mylist, result_mets_2, max_rxns, index_1_removed2, step = 3)
  # result_mets_3 <- results[["result_mets"]]
  results <- further_extend_metabolites(model, shortmodel, mylist, result_mets_all[[2]], max_rxns, index_1_removed2, step = 3)
  result_mets_all[[3]] <- results[["result_mets"]]
  index_2_removed2 <- mapply(c, index_1_removed2, results[["index_removed"]], SIMPLIFY=FALSE)
  rm(results)
  
  #step4
  cat("step 4.. \n")
  # results <- further_extend_metabolites(model, shortmodel, mylist, result_mets_3, max_rxns, index_2_removed2, step = 4)
  # result_mets_4 <- results[["result_mets"]]
  results <- further_extend_metabolites(model, shortmodel, mylist, result_mets_all[[3]], max_rxns, index_2_removed2, step = 4)
  result_mets_all[[4]] <- results[["result_mets"]]
  index_3_removed2 <- mapply(c, index_2_removed2, results[["index_removed"]], SIMPLIFY=FALSE)
  rm(results)
  
  #step5
  cat("step 5.. \n")
  # results <- further_extend_metabolites(model, shortmodel, mylist, result_mets_4, max_rxns, index_3_removed2, step = 5)
  # result_mets_5 <- results[["result_mets"]]
  results <- further_extend_metabolites(model, shortmodel, mylist, result_mets_all[[4]], max_rxns, index_3_removed2, step = 5)
  result_mets_all[[5]] <- results[["result_mets"]]
  index_4_removed2 <- mapply(c, index_3_removed2, results[["index_removed"]], SIMPLIFY=FALSE)
  rm(results)
  
  # index_1 <- mylist
  # index_1_removed <- mylist
  # index_1_removed2 <- mylist
  # for(i in 1:length(max_rxns)){
  #   if (nrow(result_mets_1[[i]]) !=0){
  #     index_1[[i]] <- which(as.vector(unlist(model$mets)) %in% as.vector(result_mets_1[[i]][,"met_short"]))
  #     index_1_removed[[i]] = removeMetsFromSet(index_1[[i]],model)
  #   }
  # }
  # 
  # index_1_removed2 <- lapply(1:length(max_rxns), function(x) {
  #   if(nrow(result_mets_1[[x]]) != 0){
  #     min_rxn <- ifelse(x==1, 0, max_rxns[x-1])
  #     max_rxn <- max_rxns[x]
  #     tmpindex <- filter_indices(model, shortmodel, index_1_removed[[x]], max_rxns = max_rxn, min_rxns = min_rxn)
  #     # Also skip any indices that were already extended in previous step
  #     tmpindex[!tmpindex %in% index_primary_removed2[[x]]]
  #   } 
  # })
  # 
  # result_mets_2 <- extend_metabolites(model, max_rxns, index = index_1_removed2, mylist = mylist, step = 2)

  # if(length(index_1_removed) == 0) {
  #   save_files(result_mets, result_mets_1, result_mets_2, result_mets_3, result_mets_4)
  #   return()
  # }
  
  # index_1_removed <- filter_indices(model, shortmodel, index_1_removed, max_rxns)
  
  # if(length(index_1_removed) == 0) {
  #   save_files(result_mets, result_mets_1, result_mets_2, result_mets_3, result_mets_4)
  #   return()
  # }

  # met_short_recon = as.vector(unlist(model$mets[index_1_removed]))
  # index_path = NULL
  # for (i in 1:length(met_short_recon)){
  #   tmp = as.numeric(which(result_mets_1[,"met_short"] == met_short_recon[i]))
  #   if (length(tmp)>1) tmp=tmp[1] 
  #   index_path = c(index_path,tmp)
  # }
  # 
  # if (length(index_1_removed)>0) {
  #   result_mets_2 = getPreviousNext(model, index_1_removed, step=2, as.vector(result_mets_1[index_path,"path"]))$mets
  # }
  # if (is.null(result_mets_2)) result_mets_2 = result_mets_1
  # 
  # 
  # 
  # 
  # ###########################################################################
  # # Step = 3 ----------------------------------------------------------------
  # ###########################################################################
  # cat("step 3.. ")
  # 
  # index_2 = which(as.vector(unlist(model$mets)) %in% as.vector(result_mets_2[,"met_short"]))
  # index_2_removed = removeMetsFromSet(index_2,model)
  # 
  # if(length(index_2_removed) == 0) {
  #   save_files(result_mets, result_mets_1, result_mets_2, result_mets_3, result_mets_4)
  #   return()
  # }
  # 
  # index_2_removed <- filter_indices(model, shortmodel, index_2_removed, max_rxns)
  # 
  # if(length(index_2_removed) == 0) {
  #   save_files(result_mets, result_mets_1, result_mets_2, result_mets_3, result_mets_4)
  #   return()
  # }
  # 
  # met_short_recon = as.vector(unlist(model$mets[index_2_removed]))
  # index_path = NULL
  # for (i in 1:length(met_short_recon)){
  #   tmp = as.numeric(which(result_mets_2[,"met_short"] == met_short_recon[i]))
  #   if (length(tmp)>1) tmp=tmp[1] 
  #   index_path = c(index_path,tmp)
  # }
  # 
  # if (length(index_2_removed)>0) {
  #   result_mets_3 = getPreviousNext(model, index_2_removed, step=3, as.vector(result_mets_2[index_path,"path"]))$mets
  # }
  # if (is.null(result_mets_3)) result_mets_3 = result_mets_2
  # 
  # 
  # 
  # 
  # ###########################################################################
  # # Step = 4 ----------------------------------------------------------------
  # ###########################################################################
  # cat("step 4.. ")
  # 
  # index_3 = which(as.vector(unlist(model$mets)) %in% as.vector(result_mets_3[,"met_short"]))
  # index_3_removed = removeMetsFromSet(index = index_3,model)
  # 
  # if(length(index_3_removed) == 0) {
  #   save_files(result_mets, result_mets_1, result_mets_2, result_mets_3, result_mets_4)
  #   return()
  # }
  # 
  # index_3_removed <- filter_indices(model, shortmodel, index_3_removed, max_rxns)
  # 
  # if(length(index_3_removed) == 0) {
  #   save_files(result_mets, result_mets_1, result_mets_2, result_mets_3, result_mets_4)
  #   return()
  # }
  # 
  # met_short_recon = as.vector(unlist(model$mets[index_3_removed]))
  # index_path = NULL
  # for (i in 1:length(met_short_recon)){
  #   tmp = as.numeric(which(result_mets_3[,"met_short"] == met_short_recon[i]))
  #   if (length(tmp)>1) tmp=tmp[1] 
  #   index_path = c(index_path,tmp)
  # }
  # 
  # if (length(index_3_removed)>0) {
  #   result_mets_4 = getPreviousNext(model, index_3_removed, step=4, as.vector(result_mets_3[index_path,"path"]))$mets
  # }
  # if (is.null(result_mets_4)) result_mets_4 = result_mets_3
  # 
  # 
  # 
  
  ###########################################################################
  # Collate results ---------------------------------------------------------
  ###########################################################################
  
  # result_mets_4 <- rbind(result_mets,result_mets_1,result_mets_2,result_mets_3,result_mets_4)
  # result_mets_3 <- rbind(result_mets,result_mets_1,result_mets_2,result_mets_3)
  # result_mets_2 <- rbind(result_mets,result_mets_1,result_mets_2)
  # result_mets_1 <- rbind(result_mets,result_mets_1)
  # result_mets_0 <- result_mets
  
  save_files(result_mets, result_mets_all)
  # save_files(result_mets, result_mets_1, result_mets_2, result_mets_3, result_mets_4, result_mets_5)
  return()
  

  
  # save(result_mets_4, file=paste(outdir,"mss_4",paste(gene_in, "RData", sep="."), sep="/"))
  # save(result_mets_3, file=paste(outdir,"mss_3",paste(gene_in, "RData", sep="."), sep="/"))
  # save(result_mets_2, file=paste(outdir,"mss_2",paste(gene_in, "RData", sep="."), sep="/"))
  # save(result_mets_1, file=paste(outdir,"mss_1",paste(gene_in, "RData", sep="."), sep="/"))
  # save(result_mets_0, file=paste(outdir,"mss_0",paste(gene_in, "RData", sep="."), sep="/"))
  
}  