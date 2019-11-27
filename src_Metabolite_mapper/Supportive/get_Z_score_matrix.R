get_Z_score_matrix <- function(xls_data, patientID, data_location, bad_mets){
  old_patient_number <- sub(xls_data[grep(patientID, xls_data$PatientID, fixed = TRUE)[1], Old.patient.number], pattern = "\\..*", replacement = "")
  
  RDS_file <- list.files(path = data_location, pattern = "*.RDS")
  Zint_values <- readRDS(file = paste(data_location, RDS_file, sep = "/"))
  Zint_values <- as.data.table(Zint_values)
  Zint_values[nchar(HMDB_code) == 9, HMDB_code := str_replace(HMDB_code, pattern = "HMDB", replacement = "HMDB00")]
  
  # Remove any metabolites that (for some reason) should not be present
  Zint_values <- Zint_values[!HMDB_code %in% bad_mets]
  
  # Collate indistinguishable metabolites
  columns_to_compare <- grep("[PC][0-9]+\\.[0-9]", colnames(Zint_values), value = TRUE)
  tmp <- as.data.frame(apply(Zint_values[,..columns_to_compare], 1, paste, collapse = ","), stringsAsFactors = FALSE)
  rownames(tmp) <- Zint_values$HMDB_code
  colnames(tmp) <- "values"
  tmp <- aggregate(rownames(tmp), by=tmp['values'], paste, collapse = ",")
  rownames(tmp) <- tmp$x
  tmp$x <- NULL
  tmp <- tidyr::separate(tmp, values, into = columns_to_compare, sep = ",", convert = TRUE)
  Zint_pruned <- as.matrix(tmp[order(rownames(tmp)),])
  rm(tmp)
  
  #  Get patient specific columns and average Z-scores if #DBS > 1
  DBS <- xls_data[tolower(PatientID) == tolower(patientID), Old.patient.number]
  
  if(length(DBS) == 1){
    Z_scores_matrix <- Zint_pruned[,grep(paste0(DBS,"_Zscore"), colnames(Zint_pruned)), drop = FALSE]
  } else {
    Z_cols <- unlist(lapply(paste0(DBS,"_Zscore"), function(x) grep(x, colnames(Zint_pruned))))
    Z_scores_matrix <- as.matrix(Zint_pruned[,Z_cols])
  }
  
  return(Z_scores_matrix)
}


