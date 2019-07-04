getPvalues <- function(peaklist, n_patients, n_controls, assi.lab, patients, controls, adducts=FALSE){
  # peaklist = outlist.pos.stats.more
  # assi.lab="HMDB_code"
  # adducts=FALSE
  
  if (adducts){
    n_rows <- dim(peaklist)[1] 
    p.values = matrix(NA, nrow = n_rows, ncol = 2*n_patients+n_controls)
    rownames(p.values) = peaklist[,assi.lab]
    # int.index = c(2:(grep("HMDB_name", colnames(peaklist))-1))
    # int.index = c(2:(grep("HMDB_code", colnames(peaklist))-1)) # Changed "name" to "code"
    # int.index = c((grep("theormz_noise", colnames(peaklist),fixed = TRUE)+1):(grep("avg.ctrls", colnames(peaklist))-1))
    # get only colnames starting with C or P, without _ characters after the digits.
    int.index <- grep("^(?=[CP]\\d+.\\d)(?!.*_)", colnames(peaklist), perl = TRUE)
    # int.index <- grep("[^CP]", colnames(peaklist))
  } else {
    # assigned mets in profile
    row_index=which(peaklist[,assi.lab]!="")  # ME assi.lab= HMDB_code kolom
    n_rows <- length(row_index)
    p.values = matrix(NA, nrow = n_rows, ncol = 2*n_patients+n_controls)
    rownames(p.values) = peaklist[row_index,assi.lab]
    # int.index = c((grep("theormz_noise", colnames(peaklist),fixed = TRUE)+1):(grep("avg.ctrls", colnames(peaklist))-1))
    # get only colnames starting with C or P, without _ characters after the digits.
    int.index <- grep("^(?=[CP]\\d+.\\d)(?!.*_)", colnames(peaklist), perl = TRUE)
    # int.index <- grep("[^CP]", colnames(peaklist))
  }

  p.labels =  rep("", n_patients)
  int.labels =  rep("", n_patients)
  int.labels.c =  rep("", n_controls)
  z.index = grep("Zscore", colnames(peaklist))

  # loop over assigned mets
  for (i in 1:n_rows){
    
    if (adducts){
      z.scores = peaklist[i, z.index]
      intensities = peaklist[i, int.index]
    } else {
      # Check Z-score
      z.scores = peaklist[row_index[i], z.index]
      intensities = peaklist[row_index[i], int.index]
    }  

    p.values.row = rep(NA, n_patients)
    int.row = rep(NA, n_patients)
    int.row.c = rep(NA, n_controls)
    
    # loop over patients    
    for (j in 1:length(patients)){
      
      # average biological replicates
      # z.average = mean(as.vector(unlist(z.scores[grep(paste("P",j,"_",sep=""), names(z.scores))])))
      # z.average = mean(as.vector(unlist(z.scores[grep(paste("P",sprintf("%02d", patients[j]),".",sep=""), names(z.scores),fixed = TRUE)])))
      z.average = mean(as.vector(unlist(z.scores[grep(paste("P",patients[j],".",sep=""), names(z.scores),fixed = TRUE)])))
      
      # int.row[j] = mean(as.vector(unlist(intensities[grep(paste("P",j,"_",sep=""), names(intensities))])))
      # int.row[j] = mean(as.numeric(as.vector(unlist(intensities[grep(paste("P",sprintf("%02d", patients[j]),".",sep=""), names(intensities),fixed = TRUE)]))))
      int.row[j] = mean(as.numeric(as.vector(unlist(intensities[grep(paste("P",patients[j],".",sep=""), names(intensities),fixed = TRUE)]))))
      
      p.values.row[j] = z.average
      
      if (i==1) p.labels[j] = paste("z.average_P", patients[j], sep="")
      if (i==1) int.labels[j] = paste("P", patients[j], sep="")
    }
    
    # loop over controls    
    for (j in 1:length(controls)){
      
      # average biological replicates
      # int.row.c[j] = mean(as.vector(unlist(intensities[grep(paste("C",j,"_",sep=""), names(intensities))])))
      # int.row.c[j] = mean(as.numeric(as.vector(unlist(intensities[grep(paste("C",sprintf("%02d", controls[j]),".",sep=""), names(intensities))]))))
      int.row.c[j] = mean(as.numeric(as.vector(unlist(intensities[grep(paste("C",controls[j],".",sep=""), names(intensities),fixed = TRUE)]))))
      if (i==1) int.labels.c[j] = paste("C", controls[j], sep="")
    }
    
    p.values[i,] = c(p.values.row, int.row, int.row.c) 
  }
  
  colnames(p.values) = c(p.labels, int.labels, int.labels.c)
  
  return(p.values)
}
