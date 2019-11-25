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
name_patient_data <- "Crossomics_DBS_Marten_Training_Validation_updated20191125"
name_new_patient_data <- "Crossomics_DBS_Marten_TraVal_Inclusion_only_updated20191125"



##### Make RData file from xlsx file
path_of_xls <- paste0(dirname(rstudioapi::getActiveDocumentContext()$path),"/../Data/")
xls_df <- read.xlsx(paste0(path_of_xls,name_patient_data,".xlsx"), sheetIndex = 1)
xls_df <- xls_df[,colSums(is.na(xls_df))<nrow(xls_df)]
xls_DT <- data.table(xls_df)
xls_data <- xls_DT
xls_data <- taRifx::remove.factors(xls_data)
xls_data <- xls_data[ !grepl("Exclusion", Judith) & Gene != "NA" ]

save(file = paste0(path_of_xls,name_new_patient_data,".RData"), xls_data)



##### Check disease genes and print which ones don't have a metabolite set
dum_list <- list.files(paste0(path_of_xls, date_of_metsets))
smal_maxrxn <- paste0("/maxrxn",sort(as.numeric(unlist(strsplit(dum_list[grep("maxrxn", dum_list)], split = "maxrxn"))))[1], "/mss_0")
avail_genes <- unlist(strsplit(list.files(paste0(path_of_xls,date_of_metsets,smal_maxrxn)), ".RDS"))

gene_col <- grep("Gene", colnames(xls_DT), value = TRUE)
if(length(gene_col) != 1){
  stop("There are multiple or no columns with the term 'Gene' as name. Check whether a disease gene column is present and change its name")
} else {
  disease_genes <- unique(trimws(unlist(strsplit(as.character(xls_DT[[gene_col]]), split = ";"))))
  missing_genes <- disease_genes[!disease_genes %in% avail_genes]
  message("these genes are listed in the patient dataset, but not present as having a metabolite set:\n", paste(missing_genes, collapse = "\n"))
}

