# performMSEA <- function(metaboliteSet, av_int_and_z_values_matrix, patient, gene_in, thresh_F_pos, thresh_F_neg, path = NULL, top = 20, id="hmdb", patient_folder, plot = TRUE){
# 
#   
#   label <- paste0("av.z_", patient)
#   patient_z_values = av_int_and_z_values_matrix[,label]

performMSEA <- function(metaboliteSet, patient_z_values, thresh_F_pos, thresh_F_neg){
  

  # patient_z_values <- as.numeric(patient_z_values)
  
  # This is logFC_weight check, I think there shouldn't be any is.na(p)
  # patient_z_values = patient_z_values[!is.na(patient_z_values)]
  
  # Get boolean of all HMDB codes indicating exceding either one of the threshold values
  allMetsExcedingThres <- patient_z_values > thresh_F_pos | patient_z_values < thresh_F_neg
  
  # split_names <- strsplit(rownames(av_int_and_z_values_matrix), split = ",")
  
  # Get location of the metaboliteSet metabolites in the intensity/Z-scores matrix.
  indices <- NULL
  # index <- NULL
  for(row_num in 1:nrow(metaboliteSet)){
    metabolites <- metaboliteSet[row_num, "alt_hmdb"]
    # indices <- c(indices, which(metabolites == rownames(av_int_and_z_values_matrix)))
    indices <- c(indices, which(metabolites == rownames(patient_z_values)))
    # genes <- trimws(as.vector(unlist(strsplit(metaboliteSet[row_num,"hmdb"], split = ","))))
    # if(length(genes) > 1){
    #   for(g in genes){
    #     res <- lapply(split_names, function(sp_nam) grep(g, sp_nam))
    #     tmp_index <- which(sapply(res, function(x) length(x) > 0))
    #     index <- c(index, tmp_index)
    #     rm(tmp_index)
    #   }
    #   index <- unique(index)
    #   if(length(index) > 1) stop(paste("duplicate metabolites not duplicate in Z-score dataset for gene:", gene_in))
    # } else {
    #   res <- lapply(split_names, function(sp_nam) grep(genes, sp_nam))
    #   index <- which(sapply(res, function(x) length(x) > 0))
    # }
    # 
    # 
    # 
    # # which vectors contain a search term
    # if(length(index) > 0){
    #   indices <- c(indices, index)
    # } else {
    #   indices <- c(indices, NA)
    # }
    # index <- NULL
  }
  
  # Check if there are any metabolites actually present in both the sets and intensity scores
  if(all(is.na(indices))) return(list("p.value"= 1))
  
  # Make single rows of all rows that contain indistinguishable HMDBs.
  # MetSet_short <- cbind(metaboliteSet, "alt_HMDB" = rownames(av_int_and_z_values_matrix)[indices])
  # MetSet_short <- MetSet_short[!is.na(MetSet_short[,"alt_HMDB"]),, drop = FALSE]
  
  # newdf <- MetSet_short[!duplicated(MetSet_short[,"alt_HMDB"]),, drop = FALSE]
  # newdf <- newdf[order(newdf[,"alt_HMDB"]), , drop = FALSE]
  # MetSet_short[is.na(MetSet_short)] <- "NA"
  metaboliteSet[is.na(metaboliteSet)] <- "NA"
  
  # if(any(duplicated(MetSet_short[,"alt_HMDB"]))){
  #   for(colname in colnames(newdf)[1:length(colnames(newdf))-1]){
  #     newdf[, colname] <- aggregate(MetSet_short[,colname]~alt_HMDB, data=MetSet_short, toString, drop = FALSE)[,2]
  #   }
  # }

  
  # index <- names(patient_z_values) %in% metaboliteSet[,"alt_hmdb"]
  index <- rownames(patient_z_values) %in% metaboliteSet[,"alt_hmdb"]
  # index <- names(patient_z_values) %in% newdf[,"alt_HMDB"]
  # metsInset <- data.frame("InSetExceedingThresh" = allMetsExcedingThres[index])
  metsInset <- as.data.frame(allMetsExcedingThres[index, , drop = FALSE])
  colnames(metsInset) <- "InSetExceedingThresh"
  # if(nrow(metsInset) == 0)   return(list("p.value"= 1))
  
  # metNames <- newdf[match(newdf[,"alt_HMDB"], rownames(metsInset)),"met_long"]
  metNames <- metaboliteSet[match(metaboliteSet[,"alt_hmdb"], rownames(metsInset)),"met_long"]
  
  
  # TRYING TO GET THESE PATHS AND HMDB_SET FUNCTION SIMPLER AND WORKING WITH THE NEW STRUCTURE
  # This may be enough to replace the cluster#@-function of paths underneath...
  # I think this is only possible atm because I already removed duplicate hmdb names and have easier rownames (already adducts).
  # paths <- newdf[match(rownames(metsInset), newdf[,"alt_HMDB"]),"path"]
  # hmdb_set <- newdf[match(rownames(metsInset), newdf[,"alt_HMDB"]),"hmdb"]
  hmdb_set <- metaboliteSet[match(rownames(metsInset), metaboliteSet[,"alt_hmdb"]),"hmdb"]
  
  
  # Used for weighted Fishers test
  # average of the control intensity values and paste it to the end of the previous matrix
  # avg.int.c <- apply(av_int_and_z_values_matrix[,grep("C", colnames(av_int_and_z_values_matrix), fixed = TRUE)], 1, mean)
  # av_int_and_z_values_matrix <- cbind(av_int_and_z_values_matrix,"avg.int.controls"=avg.int.c)
  # metsInset=cbind(metsInset, 
  #                 "names"= metNames, 
  #                 "hmdb_set"= hmdb_set, 
  #                 "z-score"=patient_z_values[rownames(metsInset)], 
  #                 "path"=paths)
  metsInset=cbind(metsInset, 
                  "names"= metNames, 
                  "hmdb_set"= hmdb_set, 
                  "z-score"=patient_z_values[match(rownames(metsInset), rownames(patient_z_values))])
  
  
  
  # Discard double rows
  # metsInset = unique(metsInset, drop=FALSE)

  InSetExceedingThresh=sum(metsInset[,"InSetExceedingThresh"])
  inSetBelowThresh=sum(!metsInset[,"InSetExceedingThresh"])
  notInSetExceedingThresh <- sum(allMetsExcedingThres) - InSetExceedingThresh
  notInSetBelowThresh <- sum(!allMetsExcedingThres) - inSetBelowThresh
  


  # Fishers exact test 
  enrichment =
    matrix(c(InSetExceedingThresh, inSetBelowThresh, notInSetExceedingThresh, notInSetBelowThresh),
           nrow = 2,
           dimnames = list(c("exceeding_thresh", "below_thresh"),c("in_set", "not_in_set")))
  retVal = fisher.test(enrichment, alternative = "greater")
  p = retVal$p.value

  
  return(list("p.value"=p, "mets_exc_thres"=InSetExceedingThresh, "hmdb_set"=hmdb_set))
}
