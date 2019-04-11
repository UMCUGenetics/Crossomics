getZvalues <- function(peaklist, patients, controls){

  peaklist[,c(5:ncol(peaklist))] <- apply(peaklist[c(5:ncol(peaklist))], MARGIN = 2, as.numeric)

  z.values.col <- matrix(data = NA, nrow = nrow(peaklist), ncol = length(patients))
  int.values.col <- matrix(data = NA, nrow = nrow(peaklist), ncol = length(patients))
  int.values.col.ctrl <- matrix(data = NA, nrow = nrow(peaklist), ncol = length(controls))
  
  for (j in 1:length(patients)){
    if(length(grep(paste0("P",patients[j],".*_Zscore"), colnames(peaklist))) == 1){
      z.values.col[,j] <- peaklist[,grep(paste0("P",patients[j],".*_Zscore"), colnames(peaklist))]
      int.values.col[,j] <- peaklist[,grep(paste0("P",patients[j],".[0-9]$"), colnames(peaklist))]
    } else {
      z.values.col[,j] <- rowMeans(peaklist[,grep(paste0("P",patients[j],".*_Zscore"), colnames(peaklist))])
      int.values.col[,j] <- rowMeans(peaklist[,grep(paste0("P",patients[j],".[0-9]$"), colnames(peaklist))])
    }
    colnames(z.values.col) <- paste0("z.average_P", patients)
    rownames(z.values.col) <- rownames(peaklist)
    colnames(int.values.col) <- paste0("av.intensity_P",patients)
    rownames(int.values.col) <- rownames(peaklist)
  }  
  for (j in 1:length(controls)){
    if(length(grep(paste0("C",controls[j],".[0-9]$"), colnames(peaklist))) == 1){
      int.values.col.ctrl[,j] <- peaklist[,grep(paste0("C",controls[j],".[0-9]$"), colnames(peaklist))]
    } else {
      int.values.col.ctrl[,j] <- rowMeans(peaklist[,grep(paste0("C",controls[j],".[0-9]$"), colnames(peaklist))])
    }
    colnames(int.values.col.ctrl) <- paste0("av.intensity_C", controls)
    rownames(int.values.col.ctrl) <- rownames(peaklist)
  }  
  averaged_z_df <- cbind(z.values.col,int.values.col,int.values.col.ctrl)
  return(averaged_z_df)
}




