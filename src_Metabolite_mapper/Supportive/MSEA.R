performMSEA <- function(metaboliteSet, av_int_and_z_values_matrix, patient, gene_in, thresh_F_pos, thresh_F_neg, path = NULL, top = 20, id="hmdb", patient_folder, plot = TRUE){

    # set the dimensions of the pictures for later
  # width=1024
  # height=768
  # Fold Change is put to the power logFC_weight (weighted_logFC = round(abs(log2(foldChange))^logFC_weight))
  # logFC_weight=1.2
  
  label <- paste0("av.z_", patient)
  patient_z_values_vector = av_int_and_z_values_matrix[,label]
  
  # This is logFC_weight check, I think there shouldn't be any is.na(p)
  patient_z_values_vector = patient_z_values_vector[!is.na(patient_z_values_vector)]
  
  # Get boolean of all HMDB codes indicating exceding either one of the threshold values
  allMetsExcedingThres <- patient_z_values_vector > thresh_F_pos | patient_z_values_vector < thresh_F_neg
  
  split_names <- strsplit(rownames(av_int_and_z_values_matrix), split = ",")
  # tmp <- strsplit(rownames(metsInset), split = ",")
  
  # Get location of the metaboliteSet metabolites in the intensity/Z-scores matrix.
  indices <- NULL
  index <- NULL
  for(row_num in 1:nrow(metaboliteSet)){
    # res <- lapply(split_names, function(sp_nam) grep(metaboliteSet[row_num,"hmdb"], sp_nam))
    genes <- trimws(as.vector(unlist(strsplit(metaboliteSet[row_num,"hmdb"], split = ","))))
    if(length(genes) > 1){
      for(g in genes){
        res <- lapply(split_names, function(sp_nam) grep(g, sp_nam))
        tmp_index <- which(sapply(res, function(x) length(x) > 0))
        index <- c(index, tmp_index)
        rm(tmp_index)
      }
      index <- unique(index)
      if(length(index) > 1) stop(paste("duplicate metabolites not duplicate in Z-score dataset for gene:", gene_in))
    } else {
      res <- lapply(split_names, function(sp_nam) grep(genes, sp_nam))
      index <- which(sapply(res, function(x) length(x) > 0))
    }

    
    
    # which vectors contain a search term
    # index <- which(sapply(res, function(x) length(x) > 0))
    if(length(index) > 0){
      indices <- c(indices, index)
    } else {
      indices <- c(indices, NA)
    }
    index <- NULL
  }
  
  # Check if there are any metabolites actually present in both the sets and intensity scores
  # if(length(indices) == 0) return(list("p.value"= 1))
  if(all(is.na(indices))) return(list("p.value"= 1))
  
  # Make single rows of all rows that contain indistinguishable HMDBs.
  # test <- testit::has_warning(cbind(metaboliteSet, "alt_HMDB" = rownames(av_int_and_z_values_matrix)[indices]))
  # if(test){cat(gene_in)}
  MetSet_short <- cbind(metaboliteSet, "alt_HMDB" = rownames(av_int_and_z_values_matrix)[indices])
  MetSet_short <- MetSet_short[!is.na(MetSet_short[,"alt_HMDB"]),, drop = FALSE]
  
  newdf <- MetSet_short[!duplicated(MetSet_short[,"alt_HMDB"]),, drop = FALSE]
  newdf <- newdf[order(newdf[,"alt_HMDB"]), , drop = FALSE]
  MetSet_short[is.na(MetSet_short)] <- "NA"
  
  if(any(duplicated(MetSet_short[,"alt_HMDB"]))){
    for(colname in colnames(newdf)[1:length(colnames(newdf))-1]){
      newdf[, colname] <- aggregate(MetSet_short[,colname]~alt_HMDB, data=MetSet_short, toString, drop = FALSE)[,2]
    }
  }

  
  
  index <- names(patient_z_values_vector) %in% newdf[,"alt_HMDB"]
  metsInset <- data.frame("InSetExceedingThresh" = allMetsExcedingThres[index])
  if(nrow(metsInset) == 0)   return(list("p.value"= 1))
  
  metNames <- newdf[match(newdf[,"alt_HMDB"], rownames(metsInset)),"met_long"]
  
  
  # TRYING TO GET THESE PATHS AND HMDB_SET FUNCTION SIMPLER AND WORKING WITH THE NEW STRUCTURE
  # This may be enough to replace the cluster#@-function of paths underneath...
  # I think this is only possible atm because I already removed duplicate hmdb names and have easier rownames (already adducts).
  paths <- newdf[match(rownames(metsInset), newdf[,"alt_HMDB"]),"path"]
  hmdb_set <- newdf[match(rownames(metsInset), newdf[,"alt_HMDB"]),"hmdb"]
  
  
  # Used for weighted Fishers test
  # average of the control intensity values and paste it to the end of the previous matrix
  avg.int.c <- apply(av_int_and_z_values_matrix[,grep("C", colnames(av_int_and_z_values_matrix), fixed = TRUE)], 1, mean)
  av_int_and_z_values_matrix <- cbind(av_int_and_z_values_matrix,"avg.int.controls"=avg.int.c)
  # foldChange = av_int_and_z_values_matrix[,2]/avg.int.c # Dangerous, logFC_weight value of 1 means no change
  # weighted_logFC = round(abs(log2(foldChange))^logFC_weight)
  metsInset=cbind(metsInset, 
                  "names"= metNames, 
                  "hmdb_set"= hmdb_set, 
                  # "weighted_logFC"=weighted_logFC[index],
                  # "FC"=foldChange[index],
                  "z-score"=patient_z_values_vector[rownames(metsInset)], 
                  "path"=paths)
  
  # Discard double rows
  metsInset = unique(metsInset, drop=FALSE)

  InSetExceedingThresh=sum(metsInset[,"InSetExceedingThresh"])
  inSetBelowThresh=sum(!metsInset[,"InSetExceedingThresh"])
  notInSetExceedingThresh <- sum(allMetsExcedingThres) - InSetExceedingThresh
  notInSetBelowThresh <- sum(!allMetsExcedingThres) - inSetBelowThresh
  


  # Fishers exact test ######################################################################
  ###########################################################################################
  enrichment =
    matrix(c(InSetExceedingThresh, inSetBelowThresh, notInSetExceedingThresh, notInSetBelowThresh),
           nrow = 2,
           dimnames = list(c("exceeding_thresh", "below_thresh"),c("in_set", "not_in_set")))
  retVal = fisher.test(enrichment, alternative = "greater")
  p = retVal$p.value
  ###########################################################################################
    
  # 
  # # Code for creating just the metsInset exceding the thresholds
  # metsInset <- metsInset[metsInset$InSetExceedingThresh,]
  # 
  # ints <- av_int_and_z_values_matrix[rownames(metsInset), -grep("av.z", colnames(av_int_and_z_values_matrix), fixed=TRUE),drop=FALSE]
  # z_values <- av_int_and_z_values_matrix[rownames(metsInset), grep("av.z", colnames(av_int_and_z_values_matrix), fixed=TRUE),drop=FALSE]
  # 
  # 
  # rownames(ints) = metsInset[,"names"]
  # rownames(z_values) = metsInset[,"names"]
  # 
  # ints <- ints[, -grep("avg.int.controls", colnames(ints), fixed=TRUE),drop=FALSE]
  # 
  # # Remove negative and NA values
  # ints[ints<0] <- NA
  # ints[is.na(ints)] <- 0
  # 
  # # Whole column zeros
  # remove <- which(colSums(ints)==0)
  # if (length(remove)>0) ints = ints[,-as.numeric(remove), drop=FALSE]
  # 
  # if (length(z_values) > 1) {
  #   
  #   if(!plot){
  #     ####
  #     # ordering of heatmap without the heatmap (as heatmap3 would have done it; look at its source code)
  #     distfun = function(x) as.dist(1 - cor(t(x),use="pa")) # use of pa is not necessary, NA values are not present
  #     # calculate means
  #     Rowv <- rowMeans(ints)
  #     Colv <- colMeans(ints)
  #     # calculate distribution and make dencrogram of rows
  #     distMatrixR <- distfun(ints)
  #     hcr <- hclust(distMatrixR,method="complete")
  #     ddr <- as.dendrogram(hcr)
  #     #reorder dendrogram
  #     reorderfun = function(d, w) reorder(d, w)
  #     ddr <- reorderfun(ddr, Rowv)
  #     
  #     # Same for columns:
  #     distMatrixC <- distfun(t(ints))
  #     hcc <- hclust(distMatrixC,method="complete")
  #     ddc <- as.dendrogram(hcc)
  #     ddc <- reorderfun(ddc, Colv)
  #     
  #     rowInd <- order.dendrogram(ddr)
  #     colInd <- order.dendrogram(ddc)
  #     
  #     # order data according to dendrogram
  #     ints <- ints[rowInd, colInd]
  #     labels = colnames(z_values)
  #     z_values = as.data.frame(z_values[rowInd,])
  #     colnames(z_values)=labels
  #     metsInset = metsInset[rowInd,]
  #     ####
  #   } else {
  #     CairoPNG(filename=paste0(path, "/", patient_folder, "/", gene_in,".png"), width, height)
  #     
  #     # hm = heatmap(as.matrix(ints),
  #     #              scale="row",
  #     #              #distfun = "euclidean", #spearman.dist,
  #     #              col=colorRampPalette(c("yellow","blue"))(100),
  #     #              margins=c(2.5,26),
  #     #              labRow = metsInset[,"hmdb_set"])
  #     #
  #     hm <- heatmap3::heatmap3(ints,
  #                              col = colorRampPalette(c("navy", "white", "firebrick3"))(1024),
  #                              balanceColor = TRUE,
  #                              labRow = metsInset[,"hmdb_set"])
  #     #
  #     dev.off()
  #     
  #     ints = ints[hm$rowInd,hm$colInd]
  #     labels = colnames(z_values)
  #     z_values = as.data.frame(z_values[hm$rowInd,])
  #     colnames(z_values)=labels
  #     metsInset = metsInset[hm$rowInd,]
  #   }
  # }
  # 
  # ints = data.frame("compound"=rownames(ints),
  #                   "HMDB_set"=metsInset[,"hmdb_set"], 
  #                   "Z score" = as.numeric(z_values[,grep(colnames(z_values),pattern=toString(patient),fixed=TRUE)]),
  #                   "path"=metsInset[,"path"],
  #                   "HMDB"=rownames(metsInset),
  #                   ints)
  # 
  # ints=ints[order(ints[,"Z.score"]),]
  # 
  # if(nrow(ints) > 0) genExcelFileShort(as.data.frame(ints), paste(path, "/", patient_folder, "/", gene_in,".xls",sep=""))
  # if(nrow(ints) > 0) save(ints, file = paste0(path, patient_folder, "/", gene_in,".RData"))
  
  
  return(list("p.value"=p,"mets_exc_thres"=InSetExceedingThresh))
}
