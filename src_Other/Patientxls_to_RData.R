###########################################################################
# Info --------------------------------------------------------------------
###########################################################################

# This script makes an RData file from an excel file which holds all the patient information 
# for running the cross omics pipeline

# R version:  3.6.0 (2019-04-26)
# platform:   x86_64-apple-darwin15.6.0 (64-bit)
# OS:         macOS Mojave 10.14.6
# 
# libraries:


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Libraries ---------------------------------------------------------------
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

library(rstudioapi)
library(data.table)
library(xlsx)
library(taRifx) # to remove factors


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Manual changes ----------------------------------------------------------
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
date_of_metsets <- "2019-08-12"
name_patient_data <- "Crossomics_DBS_Marten_trimmed20191204"
name_new_patient_data <- "Crossomics_DBS_Marten_trimmed20191205"



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Functions ---------------------------------------------------------------
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

fix_patient_number <- function(xls_data, old_patient_column, DBS = TRUE){
  if(DBS){
    tmp <- sapply(strsplit(xls_data[, old_patient_column],"[P.]"), function(x)
      paste0("P", sprintf("%03d",as.numeric(x[2])), ".", x[3])
    )
  } else {
    tmp <- sapply(strsplit(xls_data[, old_patient_column],"[P.]"), function(x)
      paste0("P", sprintf("%03d",as.numeric(x[2])))
    )
  }
  return(tmp)
}
get_PatientID <- function(xls_data, column1, column2, seperator = "^"){
  patientID <- paste(xls_data[,column1], xls_data[,column2], sep = seperator)
  return(patientID)
}



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Make RData file from xlsx file ------------------------------------------
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

path_of_xls <- paste0(dirname(rstudioapi::getActiveDocumentContext()$path),"/../Data/")
xls_df <- read.xlsx(paste0(path_of_xls,name_patient_data,".xlsx"), sheetIndex = 1)
xls_df <- xls_df[,colSums(is.na(xls_df))<nrow(xls_df)]
# xls_DT <- data.table(xls_df)
# xls_data <- xls_DT

xls_data <- xls_df
xls_data <- taRifx::remove.factors(xls_data)
# xls_data <- xls_data[ !grepl("Exclusion", Judith) & Gene != "NA" ]

if(sum(colnames(xls_data) == "Patient number in set" | colnames(xls_data) == "Patient.number.in.set") > 0){
  colnames(xls_data)[colnames(xls_data) == "Patient number in set"] <- "Old.patient.number"
  colnames(xls_data)[colnames(xls_data) == "Patient.number.in.set"] <- "Old.patient.number"
  colnames(xls_data)[colnames(xls_data) == "Patient number in set.Iden"] <- "Old.patient.number.Iden"
  colnames(xls_data)[colnames(xls_data) == "Patient.number.in.set.Iden"] <- "Old.patient.number.Iden"
  colnames(xls_data)[colnames(xls_data) == "Patient number in set.Iden2"] <- "Old.patient.number.Iden2"
  colnames(xls_data)[colnames(xls_data) == "Patient.number.in.set.Iden2"] <- "Old.patient.number.Iden2"
  colnames(xls_data)[colnames(xls_data) == "Patient number in set.Iden3"] <- "Old.patient.number.Iden3"
  colnames(xls_data)[colnames(xls_data) == "Patient.number.in.set.Iden3"] <- "Old.patient.number.Iden3"
}

##### 
# Get correct patient numbers and IDs

xls_data$Patient.number <- fix_patient_number(xls_data, "Old.patient.number")
xls_data$Patient.number.Iden <- fix_patient_number(xls_data, "Old.patient.number.Iden")
xls_data$Patient.number.Iden2 <- fix_patient_number(xls_data, "Old.patient.number.Iden2")
xls_data$Patient.number.Iden3 <- fix_patient_number(xls_data, "Old.patient.number.Iden3")

xls_data[xls_data == "P0NA.NA"] <- NA

xls_data$Patient <- fix_patient_number(xls_data, "Old.patient.number", DBS = FALSE)
xls_data$Patient.Iden <- fix_patient_number(xls_data, "Old.patient.number.Iden", DBS = FALSE)
xls_data$Patient.Iden2 <- fix_patient_number(xls_data, "Old.patient.number.Iden2", DBS = FALSE)
xls_data$Patient.Iden3 <- fix_patient_number(xls_data, "Old.patient.number.Iden3", DBS = FALSE)

xls_data[xls_data == "P0NA"] <- NA

xls_data$PatientID <- get_PatientID(xls_data, "Dataset", "Patient")
xls_data$PatientID.Iden <- get_PatientID(xls_data, "Dataset.Iden", "Patient.Iden")
xls_data$PatientID.Iden2 <- get_PatientID(xls_data, "Dataset.Iden2", "Patient.Iden2")
xls_data$PatientID.Iden3 <- get_PatientID(xls_data, "Dataset.Iden3", "Patient.Iden3")

xls_data[xls_data == "NA^NA"] <- NA

file_patientIDs <- paste(get_PatientID(xls_data, "Patient", "Dataset", "_"), 
                         get_PatientID(xls_data, "Patient.Iden", "Dataset.Iden", "_"),
                         get_PatientID(xls_data, "Patient.Iden2", "Dataset.Iden2", "_"),
                         get_PatientID(xls_data, "Patient.Iden3", "Dataset.Iden3", "_"), sep = ";")
xls_data$file_patientIDs <- gsub(";NA_NA", "", file_patientIDs)

save(file = paste0(path_of_xls,name_new_patient_data,".RData"), xls_data)


# 
# ##### Check disease genes and print which ones don't have a metabolite set
# dum_list <- list.files(paste0(path_of_xls, date_of_metsets))
# smal_maxrxn <- paste0("/maxrxn",sort(as.numeric(unlist(strsplit(dum_list[grep("maxrxn", dum_list)], split = "maxrxn"))))[1], "/mss_0")
# avail_genes <- unlist(strsplit(list.files(paste0(path_of_xls,date_of_metsets,smal_maxrxn)), ".RDS"))
# 
# gene_col <- grep("Gene", colnames(xls_DT), value = TRUE)
# if(length(gene_col) != 1){
#   stop("There are multiple or no columns with the term 'Gene' as name. Check whether a disease gene column is present and change its name")
# } else {
#   disease_genes <- unique(trimws(unlist(strsplit(as.character(xls_DT[[gene_col]]), split = ";"))))
#   missing_genes <- disease_genes[!disease_genes %in% avail_genes]
#   message("these genes are listed in the patient dataset, but not present as having a metabolite set:\n", paste(missing_genes, collapse = "\n"))
# }
# 
