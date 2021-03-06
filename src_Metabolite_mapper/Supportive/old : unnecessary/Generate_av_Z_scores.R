# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Session info ------------------------------------------------------------
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

cat(
  "
  This file has been made a function with the same name
  
  
  Created by:   Marten Kerkhofs, 2019-04-11
  Modified by:  Marten Kerkhofs, 2019-04-18

  Copied from 'GeneMetaboliteMapper_Marten_Adductsums.R, on 2019-04-11 which in turn was copied from similarly named files:
  'GeneMetaboliteMapper_ME.R' in the Metab/Metabolomics/DIMS_pipeline/R_workspace_ME/Crossomics/Crossomics_SinglePatients/src folder

  
  OS
  macOS 10.14.4
  
  Package versions:
  R version   3.5.1 (2018-07-02)
  
  bioDist       1.54.0
  Cairo         1.5-9
  XLConnect     0.2-15
  BridgeDbR     1.16.1
  rJava         0.9-10
  XLConnectJars 0.2-15
  KernSmooth    2.23-15
  BiocGenerics  0.28.0
  
  
  USE:
  This file takes adductsums metabolite data, collates them, calculates the Z-scores for all patients and HMDBs and saves it to a new xlsx file
  
  Input:
  Varying functions in the src folder
  metabolite data from a patient: adductSums_xx.RData (xx = negative or positive) files
  
  Output
  Complete adductsums file, including all Z scores and intensities in an RDS file (single r data file).
  "
)

# wet working directory for Marten's mac:
# try(setwd("/Volumes/DATA/Metabolomics/DIMS_pipeline/R_workspace_ME/Crossomics/Crossomics_SinglePatients"), silent = TRUE)
# try(setwd("/Volumes/DATA-1/Metabolomics/DIMS_pipeline/R_workspace_ME/Crossomics/Crossomics_SinglePatients"), silent = TRUE)
# try(setwd("/Volumes/DATA-2/Metabolomics/DIMS_pipeline/R_workspace_ME/Crossomics/Crossomics_SinglePatients"), silent = TRUE)
# try(setwd("/Volumes/DATA-3/Metabolomics/DIMS_pipeline/R_workspace_ME/Crossomics/Crossomics_SinglePatients"), silent = TRUE)
try(setwd("/Volumes/Metab/Metabolomics/DIMS_pipeline/R_workspace_ME/Crossomics/Crossomics_SinglePatients"), silent = TRUE)




# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Libraries ---------------------------------------------------------------
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

source("/Users/mkerkho7/DIMS2_repo/Crossomics/Scripts/Supportive/getZvalues_11042019.R")




# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Load DIMS data ----------------------------------------------------------
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

## Load DIMS input files, don't know why this is necessary, but it gives a different result than just load()
loadRData <- function(fileName){
  load(fileName)
  get(ls()[ls() != "fileName"])
}

# Added:
# outlist.neg.stats.more <- loadRData("../../HPCxls/results_RES_DBS_20180720_Run4_Diagnosis2017_1/outlist_identified_negative.RData")
# outlist.pos.stats.more <- loadRData("../../HPCxls/results_RES_DBS_20180720_Run4_Diagnosis2017_1/outlist_identified_positive.RData")
# To get the adducts_summed data
# outlist.neg.adducts.HMDB <- loadRData("../../HPCxls/results_RES_DBS_20180720_Run4_Diagnosis2017_1/adductSums_negative.RData")
# outlist.pos.adducts.HMDB <- loadRData("../../HPCxls/results_RES_DBS_20180720_Run4_Diagnosis2017_1/adductSums_positive.RData")

# project_folder <- "/Volumes/Metab/Metabolomics/Research Metabolic Diagnostics/Metabolomics Projects/Projects 2015/Project 2015_011_SinglePatients/17 SinglePatients_XVII/Bioinformatics/"
project_folder <- "/Volumes/Metab/Metabolomics/DIMS_pipeline/R_workspace_ME/HPCxls/results_RES_DBS_20180720_Run4_Diagnosis2017_1/"
outlist.neg.adducts.HMDB <- loadRData(paste0(project_folder,"adductSums_negative.RData"))
outlist.pos.adducts.HMDB <- loadRData(paste0(project_folder,"adductSums_positive.RData"))


# This is a check and ensures the column names are in the same order.
tmp <- intersect(colnames(outlist.neg.adducts.HMDB), colnames(outlist.pos.adducts.HMDB))
outlist.neg.adducts.HMDB <- outlist.neg.adducts.HMDB[,tmp]
outlist.pos.adducts.HMDB <- outlist.pos.adducts.HMDB[,tmp]
# outlist.neg.stats.more <- loadRData("Z:/Metabolomics/DIMS_pipeline/R_workspace_ME/HPCxls/results_RES_DBS_20180720_Run4_Diagnosis2017_1/outlist_identified_negative.RData")
# outlist.pos.stats.more <- loadRData("Z:/Metabolomics/DIMS_pipeline/R_workspace_ME/HPCxls/results_RES_DBS_20180720_Run4_Diagnosis2017_1/outlist_identified_positive.RData")

#WORKS outlist.neg.stats.more <- loadRData("Z:/Metabolomics/Research Metabolic Diagnostics/Metabolomics Projects/Projects 2015/Project 2015_011_SinglePatients/06 SinglePatients_VI/Bioinformatics 20180824/outlist_identified_negative.RData")
#WORKS outlist.pos.stats.more <- loadRData("Z:/Metabolomics/Research Metabolic Diagnostics/Metabolomics Projects/Projects 2015/Project 2015_011_SinglePatients/06 SinglePatients_VI/Bioinformatics 20180824/outlist_identified_positive.RData")




# Way of getting variable control numbers without specifically knowing the numbers beforehand
controls_list <- unique(strsplit(colnames(outlist.neg.adducts.HMDB)[grep("C", colnames(outlist.neg.adducts.HMDB))], "[C.]"))
for(i in c(1:length(controls_list))){
  if(i == 1){control_numbers <- as.integer(controls_list[[i]][2])}
  else{control_numbers <- c(control_numbers, as.integer(controls_list[[i]][2]))}
}
controls <- control_numbers[order(control_numbers)]
patients_list <- unique(strsplit(colnames(outlist.neg.adducts.HMDB)[grep("P", colnames(outlist.neg.adducts.HMDB))], "[P.]"))
patients_numbers <- NULL
for(i in c(1:length(patients_list))){
  if(patients_list[[i]][1] != "") next
  # if(i == 1){patients_numbers <- as.integer(patients_list[[i]][2])}
  # else{patients_numbers <- c(patients_numbers, as.integer(patients_list[[i]][2]))}
  if(!is.null(patients_numbers)) patients_numbers <- c(patients_numbers, as.integer(patients_list[[i]][2]))
  else patients_numbers <- as.integer(patients_list[[i]][2])
}
patients <- patients_numbers[order(patients_numbers)]
genes=NULL


n_patients=length(patients)
n_controls=length(controls)




# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Collate the positive and negative adducts data --------------------------
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Determine location of shared HMDBs in both matrices so they can be added together later
index.neg <- which(rownames(outlist.neg.adducts.HMDB) %in% rownames(outlist.pos.adducts.HMDB))
index.pos <- which(rownames(outlist.pos.adducts.HMDB) %in% rownames(outlist.neg.adducts.HMDB))

# First determine positive HMDBs mutual with negative, then their common names and then all other positive HMDBs
# Necessary to remove the common names from the shared matrix-rows as they can not be added together otherwise.
tmp.pos <- outlist.pos.adducts.HMDB[rownames(outlist.pos.adducts.HMDB)[index.pos], 1:(dim(outlist.pos.adducts.HMDB)[2]-1)]
tmp.hmdb_name.pos <- outlist.pos.adducts.HMDB[rownames(outlist.pos.adducts.HMDB)[index.pos], dim(outlist.pos.adducts.HMDB)[2]]
tmp.pos.left <- outlist.pos.adducts.HMDB[-index.pos,]

# First negative HMDBs, mutual with positive, then all other negative (names are already provided by positive)
tmp.neg <- outlist.neg.adducts.HMDB[rownames(outlist.pos.adducts.HMDB)[index.pos], 1:(dim(outlist.neg.adducts.HMDB)[2]-1)]
tmp.neg.left <- outlist.neg.adducts.HMDB[-index.neg,]

# Add together the shared HMDBs, paste the other HMDBs (still including common names) underneath
tmp <- apply(tmp.pos, 2,as.numeric) + apply(tmp.neg, 2,as.numeric)
rownames(tmp) <- rownames(tmp.pos)
tmp <- cbind(tmp, "HMDB_name"=tmp.hmdb_name.pos)
# adducts.neg.pos <- rbind(tmp, tmp.pos.left,tmp.neg.left) 
outlist.adducts.HMDB <- rbind(tmp, tmp.pos.left, tmp.neg.left) 

outlist.adducts.HMDB <- cbind(outlist.adducts.HMDB, "HMDB_code"=rownames(outlist.adducts.HMDB))

# Additional filter for samples that do not follow the name format of C... or P... OR starting with HMDB
outlist.adducts.HMDB <- outlist.adducts.HMDB[,grep(c("^[PC]\\d+\\.\\d|^HMDB"), colnames(outlist.adducts.HMDB))]




# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Get Z scores ------------------------------------------------------------
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

peaklist <- as.data.frame(outlist.adducts.HMDB)

# Determine which columns contain controls and which contain all intensity values
ctrl.cols <- grep("C", colnames(peaklist), fixed = TRUE)
int.cols <- c(ctrl.cols, grep("^P", colnames(peaklist)))
# int.cols <- c(ctrl.cols, grep("P", colnames(peaklist), fixed = TRUE))

# Some sort of check to include NA's in places where there is no value
# peaklist[,int.cols][peaklist[,int.cols]==0] <- NA

# calculate mean and sd for Control group
# tmp = data.matrix(peaklist[ , ctrl.cols], rownames.force = TRUE)
tmp <- outlist.adducts.HMDB[ , ctrl.cols]
peaklist$avg.ctrls <- apply(tmp, 1, function(x) mean(as.numeric(x),na.rm = TRUE))
peaklist$sd.ctrls <- apply(tmp, 1, function(x) sd(as.numeric(x),na.rm = TRUE))

cnames.z = NULL

# calculate the actual Z scores
for (i in int.cols) {
  cname <- colnames(peaklist)[i]
  cnames.z <- c(cnames.z, paste(cname, "Zscore", sep="_"))
  zscores.1col <- (as.numeric(as.vector(unlist(peaklist[ , i]))) - peaklist$avg.ctrls) / peaklist$sd.ctrls
  peaklist <- cbind(peaklist, zscores.1col)
}

# Correct the column names
colnames(peaklist)[grep("zscores.1col", colnames(peaklist))[1]:ncol(peaklist)] <- cnames.z

# Put columns in desired order (intensities and zscores last, all other info first)
peaklist <- peaklist[,c(grep("^[CP]\\d+", colnames(peaklist), value = TRUE, invert = TRUE),grep("^[CP]\\d+", colnames(peaklist), value = TRUE))]




# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Get mean Z-values for specific Patient ----------------------------------
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Dit geeft geen p.values, maar de gemiddelde Z-score van de HMDB intensities van de patient. 
adductsHMDB_z_ave_and_int <- getZvalues(peaklist = peaklist, 
                                        patients,
                                        controls
                                        )

path <- "/Users/mkerkho7/DIMS2_repo/TestResults/"
# filelist <- "metabolites"
# wb <- loadWorkbook(paste(path,"20180720_Run4_Diagnosis2017_1_Adductsums_Zscores.xls",sep=""), create = TRUE)
# createSheet(wb, name = filelist)
# writeWorksheet(wb, adductsHMDB_z_ave_and_int, sheet = filelist)
# saveWorkbook(wb)
# rm(wb)

saveRDS(object = adductsHMDB_z_ave_and_int, file = paste(path,"20180720_Run4_Diagnosis2017_1_Adductsums_Zscores.rds",sep=""))
