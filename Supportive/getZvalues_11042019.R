getZvalues <- function(peaklist, patients, controls){

  # Make intensity and Z-scores of controls and patients numeric
  # Common features of those colnames: start with C or P, have a number.number motif present
  peaklist[,grep("^[CP]\\d+", colnames(peaklist))] <- apply(peaklist[,grep("^[CP]\\d+", colnames(peaklist))], MARGIN = 2, as.numeric)
  # peaklist[,c(5:ncol(peaklist))] <- apply(peaklist[c(5:ncol(peaklist))], MARGIN = 2, as.numeric)

  z.values.col <- matrix(data = NA, nrow = nrow(peaklist), ncol = length(patients))
  int.values.col <- matrix(data = NA, nrow = nrow(peaklist), ncol = length(patients))
  int.values.col.ctrl <- matrix(data = NA, nrow = nrow(peaklist), ncol = length(controls))
  
  # Average Z scores and intensities for patients
  for (j in 1:length(patients)){
    if(length(grep(paste0("P",patients[j],"\\.[0-9]_Zscore"), colnames(peaklist))) == 1){
    # if(length(grep(paste0("P",patients[j],".*_Zscore"), colnames(peaklist))) == 1){
      z.values.col[,j] <- peaklist[,grep(paste0("P",patients[j],"\\.[0-9]_Zscore"), colnames(peaklist))]
      # z.values.col[,j] <- peaklist[,grep(paste0("P",patients[j],".*_Zscore"), colnames(peaklist))]
      int.values.col[,j] <- peaklist[,grep(paste0("P",patients[j],"\\.[0-9]$"), colnames(peaklist))]
    } else {
      z.values.col[,j] <- rowMeans(peaklist[,grep(paste0("P",patients[j],"\\.[0-9]_Zscore"), colnames(peaklist))])
      # z.values.col[,j] <- rowMeans(peaklist[,grep(paste0("P",patients[j],".*_Zscore"), colnames(peaklist))])
      int.values.col[,j] <- rowMeans(peaklist[,grep(paste0("P",patients[j],"\\.[0-9]$"), colnames(peaklist))])
    }
    colnames(z.values.col) <- paste0("z.average_P", patients)
    rownames(z.values.col) <- rownames(peaklist)
    colnames(int.values.col) <- paste0("av.intensity_P", patients)
    rownames(int.values.col) <- rownames(peaklist)
  }  
  
  # Remove duplicated columns, a result from having multiple DBS within peaklist
  z.values.col <- z.values.col[,!duplicated(colnames(z.values.col))]
  int.values.col <- int.values.col[,!duplicated(colnames(int.values.col))]
  
  # Average intensities for controls
  for (j in 1:length(controls)){
    if(length(grep(paste0("C",controls[j],".[0-9]$"), colnames(peaklist))) == 1){
      int.values.col.ctrl[,j] <- peaklist[,grep(paste0("C",controls[j],".[0-9]$"), colnames(peaklist))]
    } else {
      int.values.col.ctrl[,j] <- rowMeans(peaklist[,grep(paste0("C",controls[j],".[0-9]$"), colnames(peaklist))])
    }
    colnames(int.values.col.ctrl) <- paste0("av.intensity_C", controls)
    rownames(int.values.col.ctrl) <- rownames(peaklist)
  }  
  
  # Collate averaged values in one dataframe
  averaged_z_df <- cbind(z.values.col,int.values.col,int.values.col.ctrl)
  return(averaged_z_df)
}




