# Script to calculate the total number of (non-) aberrant metabolites per patient


# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# SessionInfo -------------------------------------------------------------
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# R version:  3.6.1 (2019-07-05)
# platform:   x86_64-apple-darwin15.6.0 (64-bit)
# OS:         macOS Mojave 10.14.6
# 
# libraries:
# rstudioapi  0.10
# stringr     1.4.0
# data.table  1.12.6
# tidyr       1.0.0


# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Libraries and input variables -------------------------------------------
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

library("stringr") # string manipulation, add leading 0's
library("rstudioapi")
library("tidyr")
library("data.table")
library("dplyr")

thresholds <- "-1;1.5,-1.5;2,-3;3,-5;5"
code_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
outdir <- "Results/"

date_input <- "2019-08-12" # The date of the data/mss_0 etc. runs
date_run <- "2019-12-10" # The date of this run 

nr_mocks <- 200



# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Other variables ---------------------------------------------------------
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Subset_Of_Patients <- FALSE

thresh_df <- sapply(strsplit(unlist(strsplit(thresholds, split = ",")), ";"), `[`)

thresh_neg_list <- as.numeric(thresh_df[1,])
thresh_pos_list <- as.numeric(thresh_df[2,])


# Remove any metabolites that (for some reason) should not be present
# HMDB0002467 <- determined to be non-bodily substances, but created in the lab
bad_mets <- c("HMDB0002467")

train_data_name <- "Crossomics_DBS_Marten_trimmed20191205.RData"


# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Functions ---------------------------------------------------------------
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

getcols <- function(iteration){
  cols <- c("PatientID", "Old.patient.number")
  if(iteration == 1){
    return(cols)
  } else if(iteration == 2) {
    return(paste0(cols, ".Iden"))
  } else {
    return(paste0(cols, ".Iden", iteration - 1))
  }
}

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Load data ---------------------------------------------------------------
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

source(paste0(code_dir,"/../src_Metabolite_Mapper/Supportive/get_Z_score_matrix.R"))

load(paste0(code_dir,"/../Data/", train_data_name))

uni_dat_pat <- unique(xls_data$PatientID)

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Calculate #aberrant metabolites per patient, per threshold --------------
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


nr_mets_mat <- matrix(nrow = length(uni_dat_pat), ncol = length(thresh_neg_list)*2)
rownames(nr_mets_mat) <- uni_dat_pat
colnames(nr_mets_mat) <- unlist(lapply(paste(thresh_neg_list, thresh_pos_list, sep = ", "), function(x) paste(x, c("Aberrant", "Total"), sep = "; ")))

for(i in uni_dat_pat){
  # i <- uni_dat_pat[patient_number]
  
  # Find out whether patient is known under multiple IDs 
  ident_pat <- unique(na.omit(xls_data[grep(i, xls_data$PatientID, fixed = TRUE), "PatientID.Iden"]))
  if(length(ident_pat) > 0){
    ident_pat <- c(ident_pat, unique(na.omit(xls_data[grep(i, xls_data$PatientID, fixed = TRUE), "PatientID.Iden2"])))
    if(length(ident_pat) > 1){
      ident_pat <- c(ident_pat, unique(na.omit(xls_data[grep(i, xls_data$PatientID, fixed = TRUE), "PatientID.Iden3"])))
    }
  }
  is_id_pat <- length(ident_pat) > 0
  
  
  # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # Get patient and dataset -------------------------------------------------
  # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  tmp <- unlist(strsplit(i, split = "\\^"))
  comp_dataset <- tmp[1]
  comp_patient <- tmp[2]
  prim_dataset <- tmp[1]
  prim_patient <- tmp[2]
  patientID <- paste(comp_patient, comp_dataset, sep = "^")
  rm(tmp)
  
  if(is_id_pat){
    for(identicals in ident_pat){
      tmp <- unlist(strsplit(identicals, split = "\\^"))
      dataset.Iden <- tmp[1]
      patient.Iden <- tmp[2]
      rm(tmp)
      pID.Iden <- paste(patient.Iden, dataset.Iden, sep = "^")
      # comp_patient <- paste(patient, patient.Iden, sep = ";")
      # comp_dataset <- paste(dataset, dataset.Iden, sep = ";")
      comp_patient <- paste(comp_patient, patient.Iden, sep = ";")
      comp_dataset <- paste(comp_dataset, dataset.Iden, sep = ";")
      patientID <- paste(patientID, pID.Iden, sep = ";")
    }
  }
  
  cat("start patient:", comp_patient, "dataset(s):", comp_dataset, "\n")
  
  datasets <- unlist(strsplit(comp_dataset, split = ";"))
  data_locations <- xls_data$Location[match(datasets, xls_data$Dataset)]
  data_locations <- unlist(lapply(data_locations, function(x) paste(code_dir,"../Data", paste(strsplit(x, split = "\\\\")[[1]][c(6:8)], collapse = "/"), sep = "/")))
  
  # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # Obtain metabolite Z scores ----------------------------------------------
  # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  patientIDs <- c(i, ident_pat)
  
  for(nm in c(1:length(patientIDs))){
    cols <- getcols(nm)
    tmp <- get_Z_score_matrix(xls_data, 
                              patientID = patientIDs[nm], 
                              data_location = data_locations[nm], 
                              bad_mets, 
                              columns = cols)
    if(nm == 1){
      if(length(patientIDs) == 1){
        z_mat <- tmp
      } else {
        # tmp <- data.frame(cbind(tmp, HMDB_names = rownames(tmp)))
        tmp <- data.frame(tmp)
        z_df <- tmp
      }
    } else {
      z_df <- merge(z_df, tmp, by = "row.names", all = TRUE)
      rownames(z_df) <- z_df$Row.names
      z_df$Row.names <- NULL
      # colnames(z_mat) <- cnames
      if(nm == length(patientIDs)){
        z_mat <- as.matrix(z_df)
      }
    }
  }
  
  DBS <- ncol(z_mat)
  
  # Exclude all Internal Standard 'metabolites'
  z_mat <- z_mat[!grepl("HMDB[1-2]?[0-9]00000$", rownames(z_mat)), , drop = FALSE]
  
  # av_Z_scores <- as.matrix(rowMeans(z_val_df, na.rm = TRUE))
  av_Z_scores <- as.matrix(rowMeans(z_mat, na.rm = TRUE))
  colnames(av_Z_scores) <- paste0("av.z_", comp_patient)
  
  for (threshold in 1:length(thresh_pos_list)){
    # for (threshold in 1:2){
    thresh_F_pos <- thresh_pos_list[threshold]
    thresh_F_neg <- thresh_neg_list[threshold]
    threshs <- paste(c(thresh_F_neg, thresh_F_pos), collapse = ", ")
    
    column_ab <- paste(threshs, "Aberrant", sep = "; ")
    column_tot <- paste(threshs, "Total", sep = "; ")
    
    aberrant_mets <- sum(av_Z_scores < thresh_F_neg |  av_Z_scores > thresh_F_pos)
    total_mets <- nrow(av_Z_scores)
    
    nr_mets_mat[i, column_ab] <- aberrant_mets
    nr_mets_mat[i, column_tot] <- total_mets
    
  }
}

saveRDS(nr_mets_mat, file = paste0(code_dir,"/../", outdir,date_run,"/Aberrant_mets_per_patient.RDS"))
