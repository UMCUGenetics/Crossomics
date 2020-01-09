performMSEA <- function(metaboliteSet, patient_z_values, thresh_F_pos, thresh_F_neg){
  
  # Get boolean of all HMDB codes indicating exceding either one of the threshold values
  allMetsExcedingThres <- patient_z_values > thresh_F_pos | patient_z_values < thresh_F_neg
  
  
  # Get location of the metaboliteSet metabolites in the intensity/Z-scores matrix.
  indices <- NULL
  for(row_num in 1:nrow(metaboliteSet)){
    metabolites <- metaboliteSet[row_num, "alt_hmdb"]
    indices <- c(indices, which(metabolites == rownames(patient_z_values)))
  }
  
  # Check if there are any metabolites actually present in both the sets and intensity scores
  if(all(is.na(indices))) return(list("p.value"= 1))
  
  metaboliteSet[is.na(metaboliteSet)] <- "NA"
  
  # Determine which metabolites in the metaboliteSet exceed the threshold
  index <- rownames(patient_z_values) %in% metaboliteSet[,"alt_hmdb"]
  metsInset <- as.data.frame(allMetsExcedingThres[index, , drop = FALSE])
  colnames(metsInset) <- "InSetExceedingThresh"
  
  # note metabolite names (not actually used for anything at the moment)
  # metNames <- metaboliteSet[match(metaboliteSet[,"alt_hmdb"], rownames(metsInset)),"met_long"]
  
  # note which hmdbs are present
  hmdb_set <- metaboliteSet[match(rownames(metsInset), metaboliteSet[,"alt_hmdb"]),"hmdb"]
  
  metsInset=cbind(metsInset, 
                  # "names"= metNames, 
                  "hmdb_set"= hmdb_set, 
                  "z-score"=patient_z_values[match(rownames(metsInset), rownames(patient_z_values))])
  

  # Determine Fishers exact numbers
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
