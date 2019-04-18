# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Session info ------------------------------------------------------------
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

cat(
  "
  Created by:   Marcel Willemse?
  Modified by:  Marten Kerkhofs, 2019-04-08

  Copied from 'GeneMetaboliteMapper_ME_Marten.R, on 29/03/2019 which in turn was copied from:
  'GeneMetaboliteMapper_ME.R' in the Metab/Metabolomics/DIMS_pipeline/R_workspace_ME/Crossomics/Crossomics_SinglePatients/src folder
  on 26/03/2019 to investigate and annotate the functionality of the crossomics pipeline.
  Disabled all save/generateExcel functions
  
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
  This file takes a patient's metabolite and rare-gene variant data and performs metabolite set enrichment analysis
  
  Input:
  Varying functions in the src folder
  metabolite data from a patient: adductSums_xx.RData (xx = negative or positive) files
  A gene set of a single patient which contains rare-gene variants in that patient
  A patient number
  
  Output

  "
)

# wet working directory for Marten's mac:
# try(setwd("/Volumes/DATA/Metabolomics/DIMS_pipeline/R_workspace_ME/Crossomics/Crossomics_SinglePatients"), silent = TRUE)
# try(setwd("/Volumes/DATA-1/Metabolomics/DIMS_pipeline/R_workspace_ME/Crossomics/Crossomics_SinglePatients"), silent = TRUE)
# try(setwd("/Volumes/DATA-2/Metabolomics/DIMS_pipeline/R_workspace_ME/Crossomics/Crossomics_SinglePatients"), silent = TRUE)
# try(setwd("/Volumes/DATA-3/Metabolomics/DIMS_pipeline/R_workspace_ME/Crossomics/Crossomics_SinglePatients"), silent = TRUE)
try(setwd("/Volumes/Metab/Metabolomics/DIMS_pipeline/R_workspace_ME/Crossomics/Crossomics_SinglePatients"), silent = TRUE)


# Test for getting the data from P01.1, having PKU
patient_folder <- "/Volumes/Metab/Metabolomics/Research Metabolic Diagnostics/Metabolomics Projects/Projects 2017/Project 2017_008 MetabolomicsDiagnosis_DBS/RES_BSP_20170629_MetabolomicsDiagnosis_RUN0/Bioinformatics_metIS/"
setwd(patient_folder)



# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Manual input ------------------------------------------------------------
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# To what distance from the reaction should metabolites be considered?
# distance_to_gene <- 1 # Staat nog niet aan, kijk naar het begin van de "then perform MSEA"
patients=c(221)




# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Libraries ---------------------------------------------------------------
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# library("R.matlab")
# library("KEGGREST")
library("bioDist")
library("Cairo")
library("XLConnect")
library("BridgeDbR")

# source("./src/Crossomics/sourceDir.R")
# sourceDir("./src/Crossomics")
source("/Users/mkerkho7/DIMS2_repo/Crossomics/Scripts/Supportive/sourceDir.R")
sourceDir("/Users/mkerkho7/DIMS2_repo/Crossomics/Scripts/Supportive")



# The following code block is moved to Generate_av_Z_scores.R to create and save 1 excel file with all averaged Z scores of all patients.

{
  # Space for the Generate_av_Z_score file, or making a function out of it...
  
  
  
}














# # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# # Load DIMS data ----------------------------------------------------------
# # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 
# ## Load DIMS input files, don't know why this is necessary, but it gives a different result than just load()
# loadRData <- function(fileName){
#   load(fileName)
#   get(ls()[ls() != "fileName"])
# }
# 
# # Added:
# # outlist.neg.stats.more <- loadRData("../../HPCxls/results_RES_DBS_20180720_Run4_Diagnosis2017_1/outlist_identified_negative.RData")
# # outlist.pos.stats.more <- loadRData("../../HPCxls/results_RES_DBS_20180720_Run4_Diagnosis2017_1/outlist_identified_positive.RData")
# # To get the adducts_summed data
# outlist.neg.adducts.HMDB <- loadRData("../../HPCxls/results_RES_DBS_20180720_Run4_Diagnosis2017_1/adductSums_negative.RData")
# outlist.pos.adducts.HMDB <- loadRData("../../HPCxls/results_RES_DBS_20180720_Run4_Diagnosis2017_1/adductSums_positive.RData")
# 
# # This is a check and ensures the column names are in the same order.
# tmp <- intersect(colnames(outlist.neg.adducts.HMDB), colnames(outlist.pos.adducts.HMDB))
# outlist.neg.adducts.HMDB <- outlist.neg.adducts.HMDB[,tmp]
# outlist.pos.adducts.HMDB <- outlist.pos.adducts.HMDB[,tmp]
# # outlist.neg.stats.more <- loadRData("Z:/Metabolomics/DIMS_pipeline/R_workspace_ME/HPCxls/results_RES_DBS_20180720_Run4_Diagnosis2017_1/outlist_identified_negative.RData")
# # outlist.pos.stats.more <- loadRData("Z:/Metabolomics/DIMS_pipeline/R_workspace_ME/HPCxls/results_RES_DBS_20180720_Run4_Diagnosis2017_1/outlist_identified_positive.RData")
# 
# #WORKS outlist.neg.stats.more <- loadRData("Z:/Metabolomics/Research Metabolic Diagnostics/Metabolomics Projects/Projects 2015/Project 2015_011_SinglePatients/06 SinglePatients_VI/Bioinformatics 20180824/outlist_identified_negative.RData")
# #WORKS outlist.pos.stats.more <- loadRData("Z:/Metabolomics/Research Metabolic Diagnostics/Metabolomics Projects/Projects 2015/Project 2015_011_SinglePatients/06 SinglePatients_VI/Bioinformatics 20180824/outlist_identified_positive.RData")
# 
# 
# 
# 
# ########################################################################################################
# ### Declare patient names 
# #patients=c(80,82)
# #controls=c(30:34,36:45)
# # If WES data available genes=NULL, genes from WES in "./db/P28_HGNC.txt"
# #genes=NULL
# ########################################################################################################
# 
# #patients=c(82)
# #controls=c(30:34,36:45)
# #genes=NULL
# 
# # patients=c(221)
# # ME: USE ALL CONTROLS in list!
# # controls=c(63:65,69,73:77,79:81,83,88:89,96:99,102:104,107,109,111:112,114,117:118,120)
# 
# # Way of getting variable control numbers without specifically knowing the numbers beforehand
# controls_list <- unique(strsplit(colnames(outlist.neg.adducts.HMDB)[grep("C", colnames(outlist.neg.adducts.HMDB))], "[C.]"))
# for(i in c(1:length(controls_list))){
#   if(i == 1){control_numbers <- as.integer(controls_list[[i]][2])}
#   else{control_numbers <- c(control_numbers, as.integer(controls_list[[i]][2]))}
# }
# controls <- control_numbers[order(control_numbers)]
# genes=NULL
# 
# 
# n_patients=length(patients)
# n_controls=length(controls)
# 
# 
# 
# 
# ###############################################################################################
# ################  Seperate modes only selected adducts ########################################
# ###############################################################################################
# 
# 
# # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# # Collate the positive and negative adducts data --------------------------
# # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 
# # Determine location of shared HMDBs in both matrices so they can be added together later
# index.neg <- which(rownames(outlist.neg.adducts.HMDB) %in% rownames(outlist.pos.adducts.HMDB))
# index.pos <- which(rownames(outlist.pos.adducts.HMDB) %in% rownames(outlist.neg.adducts.HMDB))
# 
# # First determine positive HMDBs mutual with negative, then their common names and then all other positive HMDBs
# # Necessary to remove the common names from the shared matrix-rows as they can not be added together otherwise.
# tmp.pos <- outlist.pos.adducts.HMDB[rownames(outlist.pos.adducts.HMDB)[index.pos], 1:(dim(outlist.pos.adducts.HMDB)[2]-1)]
# tmp.hmdb_name.pos <- outlist.pos.adducts.HMDB[rownames(outlist.pos.adducts.HMDB)[index.pos], dim(outlist.pos.adducts.HMDB)[2]]
# tmp.pos.left <- outlist.pos.adducts.HMDB[-index.pos,]
# 
# # First negative HMDBs, mutual with positive, then all other negative (names are already provided by positive)
# tmp.neg <- outlist.neg.adducts.HMDB[rownames(outlist.pos.adducts.HMDB)[index.pos], 1:(dim(outlist.neg.adducts.HMDB)[2]-1)]
# tmp.neg.left <- outlist.neg.adducts.HMDB[-index.neg,]
# 
# # Add together the shared HMDBs, paste the other HMDBs (still including common names) underneath
# tmp <- apply(tmp.pos, 2,as.numeric) + apply(tmp.neg, 2,as.numeric)
# rownames(tmp) <- rownames(tmp.pos)
# tmp <- cbind(tmp, "HMDB_name"=tmp.hmdb_name.pos)
# # adducts.neg.pos <- rbind(tmp, tmp.pos.left,tmp.neg.left) 
# outlist.adducts.HMDB <- rbind(tmp, tmp.pos.left, tmp.neg.left) 
# 
# # dummy.neg <- rep(NA, dim(adducts.neg.pos)[1])
# # outlist.adducts.HMDB <- cbind("mzmed.pgrp"=dummy.neg,
# #                              "fq.best"=dummy.neg,
# #                              "fq.worst"=dummy.neg,
# #                              "nrsamples"=dummy.neg,
# #                              "mzmin.pgrp"=dummy.neg,
# #                              "mzmax.pgrp"=dummy.neg,
# #                              adducts.neg.pos)
# 
# outlist.adducts.HMDB <- cbind(outlist.adducts.HMDB, "HMDB_code"=rownames(outlist.adducts.HMDB))
# 
# 
# 
# 
# # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# # Get Z scores ------------------------------------------------------------
# # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 
# peaklist <- as.data.frame(outlist.adducts.HMDB)
# 
# # Determine which columns contain controls and which contain all intensity values
# ctrl.cols <- grep("C", colnames(peaklist), fixed = TRUE)
# int.cols <- c(ctrl.cols, grep("P", colnames(peaklist), fixed = TRUE))
# 
# # Some sort of check to include NA's in places where there is no value
# # peaklist[,int.cols][peaklist[,int.cols]==0] <- NA
# 
# # calculate mean and sd for Control group
# # tmp = data.matrix(peaklist[ , ctrl.cols], rownames.force = TRUE)
# tmp <- outlist.adducts.HMDB[ , ctrl.cols]
# peaklist$avg.ctrls <- apply(tmp, 1, function(x) mean(as.numeric(x),na.rm = TRUE))
# peaklist$sd.ctrls <- apply(tmp, 1, function(x) sd(as.numeric(x),na.rm = TRUE))
# 
# cnames.z = NULL
# 
# # calculate the actual Z scores
# for (i in int.cols) {
#   cname <- colnames(peaklist)[i]
#   cnames.z <- c(cnames.z, paste(cname, "Zscore", sep="_"))
#   zscores.1col <- (as.numeric(as.vector(unlist(peaklist[ , i]))) - peaklist$avg.ctrls) / peaklist$sd.ctrls
#   peaklist <- cbind(peaklist, zscores.1col)
# }
# 
# # Correct the column names
# colnames(peaklist)[grep("zscores.1col", colnames(peaklist))[1]:ncol(peaklist)] <- cnames.z
# 
# # Put columns in desired order (intensities and zscores last, all other info first)
# peaklist <- peaklist[,c(grep("^[CP]\\d+", colnames(peaklist), value = TRUE, invert = TRUE),grep("^[CP]\\d+", colnames(peaklist), value = TRUE))]
# 
# 
# 
# 
# # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# # Get mean Z-values for specific Patient ----------------------------------
# # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 
# adductsSummed <- TRUE
# 
# # Dit geeft geen p.values, maar de gemiddelde Z-score van de HMDB intensities van de patient. 
# adductsHMDB_z_ave_and_int <- getPvalues(peaklist = peaklist, 
#                                n_patients, 
#                                n_controls, 
#                                assi.lab = "HMDB_code",
#                                patients,
#                                controls,
#                                adducts = FALSE)



file <- "/Users/mkerkho7/DIMS2_repo/TestResults/20180720_Run4_Diagnosis2017_1_Adductsums_Zscores.rds"
adductsHMDB_z_ave_and_int <- readRDS(file)




# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Perform MSEA ------------------------------------------------------------
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Filter adducts and collate positive and negative mode -------------------
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# 
# # Filter adducts
# names = lapply(rownames(p.values.assi.pos), function(x, adducts=c(1,2), adducts_long=c("[M+Na]+","[M+K]+")) {
#   # x=rownames(p.values.assi.pos)[1]
#   compounds=unlist(strsplit(as.vector(x),";"))
#   # compounds=compounds[-1]  # since DIMS 2.1
#   index=c(1:length(compounds))
#   adducts_index = grep("_", compounds, fixed=TRUE)
#   
#   if (length(adducts_index)>0) index=index[-adducts_index]
#   
#   all_adducts=compounds[adducts_index]
#   keep=NULL
#   for (i in 1:length(all_adducts)){
#     # unlist(strsplit(as.vector(all_adducts[i]),"_"))[2] %in% adducts
#     keep = c(keep, unlist(strsplit(as.vector(all_adducts[i]),"_"))[2] %in% adducts)
#   }
#   all_adducts=all_adducts[keep]
#   
#   for (i in 1:length(adducts)){
#     all_adducts=gsub(paste("_",toString(adducts[i]),sep=""), paste(" ", adducts_long[i], sep=""), all_adducts)
#   }
#   
#   compounds=compounds[index]
#   if (length(all_adducts)!=0){
#     if (length(compounds)!=0){
#       compounds=paste(compounds, all_adducts, sep=";", collapse=";")
#     } else {
#       compounds=paste(all_adducts, collapse=";")
#     }  
#   } else if (length(compounds)>1) {
#     compounds=paste(compounds, collapse=";")
#   }
#   
#   if (length(compounds)==0) compounds="" 
#   
#   return(compounds)
# })
# 
# index=which(names=="")
# if (length(index)>0) names=names[-index]
# p.values.assi.pos.filt = p.values.assi.pos[-index,] 
# rownames(p.values.assi.pos.filt) = names
# 
# # save(p.values.assi.pos.filt, file="p.values.assi.pos.filt.RDate")
# 
# p.values.assi.neg = getPvalues(peaklist = outlist.neg.stats.more, 
#                                n_patients, 
#                                n_controls, 
#                                assi.lab = "HMDB_code",
#                                patients,
#                                controls,
#                                adducts = FALSE)  
# 
# # Filter adducts 
# names = lapply(rownames(p.values.assi.neg), function(x, adducts=c(1), adducts_long=c("[M+Cl]-")) {
#   # x=rownames(p.values.assi.pos)[1]
#   compounds=unlist(strsplit(as.vector(x),";"))
#   # compounds=compounds[-1]   # since DIMS 2.1
#   index=c(1:length(compounds))
#   adducts_index = grep("_", compounds, fixed=TRUE)
#   
#   if (length(adducts_index)>0) index=index[-adducts_index]
#   
#   all_adducts=compounds[adducts_index]
#   keep=NULL
#   for (i in 1:length(all_adducts)){
#     # unlist(strsplit(as.vector(all_adducts[i]),"_"))[2] %in% adducts
#     keep = c(keep, unlist(strsplit(as.vector(all_adducts[i]),"_"))[2] %in% adducts)
#   }
#   all_adducts=all_adducts[keep]
#   
#   for (i in 1:length(adducts)){
#     all_adducts=gsub(paste("_",toString(adducts[i]),sep=""), paste(" ", adducts_long[i], sep=""), all_adducts)
#   }
#   
#   compounds=compounds[index]
#   if (length(all_adducts)!=0){
#     if (length(compounds)!=0){
#       compounds=paste(compounds, all_adducts, sep=";", collapse=";")
#     } else {
#       compounds=paste(all_adducts, collapse=";")
#     }  
#   } else if (length(compounds)>1) {
#     compounds=paste(compounds, collapse=";")
#   }
#   
#   if (length(compounds)==0) compounds="" 
#   
#   return(compounds)
# })
# 
# index=which(names=="")
# if (length(index)>0) names=names[-index]
# p.values.assi.neg.filt = p.values.assi.neg[-index,] 
# rownames(p.values.assi.neg.filt) = names
# 
# # save(p.values.assi.neg.filt, file="p.values.assi.neg.filt.RDate")
# 
# # Positive and negative mode together
# rownames(p.values.assi.pos.filt) = paste(rownames(p.values.assi.pos.filt), "pos", sep="_") 
# rownames(p.values.assi.neg.filt) = paste(rownames(p.values.assi.neg.filt), "neg", sep="_") 
# p.values.assi.all=rbind(p.values.assi.pos.filt, p.values.assi.neg.filt)
# ###############################################################################################
# ###############################################################################################
# 

###############################################################################################
# Save/load as Rdata 
########################################################################################################
# save(p.values.assi.all, file = "./results/p.values.assi.all_CROSS-OMICS_P28.RData")
# load("Z:/Metabolomics/DIMS_pipeline/R_workspace_ME/Crossomics/Crossomics_SinglePatients/results/p.valus.assi.all_CROSS-OMICS_P80,82.RData")
########################################################################################################



# ###############################################################################################
# ################ Summed adducts ###############################################################
# ###############################################################################################
# adductsSummed=TRUE
# outlist.adducts.stats=cbind("HMDB_code"=rownames(outlist.adducts.stats), outlist.adducts.stats) 
# p.values.assi.all = getPvalues(outlist.adducts.stats, n_patients, n_controls, assi.lab="HMDB_code",patients,controls, adducts=TRUE) # assi.hmdb, HMDB_code For new identification !!!!!!!!!!!!!!!!!!!!!!!!!!
# ###############################################################################################
# ###############################################################################################
# ###############################################################################################

# dir.create("./results", showWarnings = FALSE)
# path <- "./results/crossomics_Fisher_weighted"
path <- "/Users/mkerkho7/DIMS2_repo/TestResults/"
# dir.create(path, showWarnings = FALSE)
thresh_F_pos <- 1.5
thresh_F_neg <- -1
top <- 20
# id="InChI_key"
id <- "hmdb"
genes <- NULL

###############################################################################################
# adductsSummed=FALSE!!!!!!!!!!
# Warning messages:
#   1: In data.row.names(row.names, rowsi, i) :
#   some row.names duplicated: 6 --> row.names NOT used
# 2: In data.row.names(row.names, rowsi, i) :
#   some row.names duplicated: 16 --> row.names NOT used
###############################################################################################

# ## Generate SinglePatients gene lists ###############################################################################
# # patients=c(5,28,41:45)
# for (i in 2:length(patients)){
#   getSinglePatientGeneList(patients[i])  
# }
# #####################################################################################################################

# Then performe MSEA ######################################################################################

path2 <- "../Crossomics_Build_Metabolite_Set/results/metabolite_sets_step_0,1,2,3_1.0_filter_1.1/mss_1"
# path2 = paste0("../Crossomics_Build_Metabolite_Set/results/metabolite_sets_step_0,1,2,3_1.0_filter_1.1/mss_",distance_to_gene)
step <- 1 # 1:3  # ME wordt gebruikt in regel 495 -503, waarom?
path_wg <- "./results/mss_WG_step_0" # for sampling random genes   #ME enkel voor random data set, niet WES
sets <- list.files(path = path2)
overview <- NULL # ME ?
# rank <- 0  # ME wordt niet gebruikt
# p_value_sum <- 0 # ME wordt niet gebruikt

###########################################################################################################
########################### real WES data #################################################################
###########################################################################################################
Sys.time()
p.values.assi.all <- adductsHMDB_z_ave_and_int
for (i in 1:length(patients)){
  
  ######################## MSEA Revon2 neighbourhood ######################################################
  dir.create(paste(path,"/P",patients[i], sep=""), showWarnings = FALSE)
  dir.create(paste(path,"/P",patients[i],"/Recon2", sep=""), showWarnings = FALSE)
  
  # subset to 1 patient!
  p.values=p.values.assi.all[,c(grep(paste("P",patients[i],sep=""), colnames(p.values.assi.all), fixed=TRUE),
                                grep("C", colnames(p.values.assi.all), fixed=TRUE),
                                which(colnames(p.values.assi.all)=="avg.int.controls"))]
  
  if (is.null(genes)){
    mss = read.table(paste("./db/P",patients[i],"_HGNC.txt",sep=""), header = FALSE, sep="\t")
    mss = as.vector(unlist(mss))
    mss = paste(mss,"RData",sep=".")
    # ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
    # gene_map_table = getBM(attributes=c('hgnc_symbol', 'entrezgene', 'ensembl_gene_id'),
    #                        filters = 'hgnc_symbol', values = as.vector(unlist(gene_list)), mart = ensembl)
    
  } else {
    genespp=as.vector(unlist(genes[i]))
    
    mss = NULL
    for (j in 1:length(genespp)){
      mss = c(mss, paste(genespp[j],"RData",sep="."))
    }
    
    wg = list.files(path = path_wg)
    mss = c(mss, wg[sample(1:length(wg), 100-length(genespp), replace=TRUE)])
    
    selection = which(sets %in% mss)
    mss = sets[selection]
    
    # save(mss, file = paste0("./results/mss_",genes[i][[1]][1],i,".RData"))
  }
  
  metSetResult = NULL
  nMets = NULL    # list with number of metabolites per gene, not used for any calculations, but only for output excel file.
  for (j in 1:length(mss)){
    print(mss[j])
    # Skip the gene if there is no metabolite pathway data available, elsewise, load its file
    if (!file.exists(paste(path2, mss[j], sep="/"))) next
    
    load(paste(path2, mss[j], sep="/"))
    gene_in=strsplit(mss[j], split = "." , fixed=TRUE)[[1]][1]
    
    # Don't know why this is necessary
    if (step==2){
      metaboliteSet = result_mets_2
    } else if (step==1){
      metaboliteSet = result_mets_1
    } else if (step==3){
      metaboliteSet = result_mets_3
    } else {
      metaboliteSet = result_mets_0
    }
    
    # Don't know why this is necessary, metaboliteSet should have content at this point
    if (!is.null(metaboliteSet)){
      if (is.null(dim(metaboliteSet))){
        metaboliteSet=data.frame(t(metaboliteSet))  # ME matrix transpose
      }
    }
    
    ##################################################
    # temporarily work around to be fixed in findMetabolicEnvironment <- why?
    # NA's in metaboliteSet[,"hmdb"] are set to 'character(0)'
    # index = which(is.na(metaboliteSet[,"hmdb"]))  
    # if (length(index) > 0) metaboliteSet[index,"hmdb"] = "character(0)"
    
    # Set all NA's to "character (0)" in the ID columns
    metaboliteSet[,c("hmdb","kegg","chebi")][is.na(metaboliteSet[,c("hmdb","kegg","chebi")])] <- "character(0)"

    ##################################################
    
    ########### Recon 2.0 ############################
    # If hmdb codes are absent, try to translate KEGG and later ChEBI codes to hmdb.
    index = which(metaboliteSet[,"hmdb"] == "character(0)") # ME list positions of "Unkown" in hmdb column
    index.sub = which(metaboliteSet[index,"kegg"] != "character(0)") # ME Welke unknown hmdb zijn niet unknown in KEGG
    kegg_id = metaboliteSet[index[index.sub],"kegg"]
    
    # ADDED (new .bridge file)
    mapper <- BridgeDbR::loadDatabase("/Users/mkerkho7/DIMS2_repo/Crossomics/metabolites_20190207.bridge")
    # mapper = loadDatabase("./db/BridgeDB/metabolites_20150717.bridge")
    hmdb = BridgeDbR::getSystemCode("HMDB")
    kegg = BridgeDbR::getSystemCode("KEGG Compound")
    
    # Try to fill in empty HMDB IDs via the KEGG ID
    if (length(kegg_id) > 0){
      for (k in 1:length(kegg_id)){
        if (!is.null(unlist(BridgeDbR::map(mapper, kegg, kegg_id[k], hmdb)[1]))) {
          metaboliteSet[index[index.sub[k]],"hmdb"] <- unlist(BridgeDbR::map(mapper, kegg, kegg_id[k], hmdb)[1])
        }
      }
    }
    
    index = which(metaboliteSet[,"hmdb"] == "character(0)")
    index.sub = which(metaboliteSet[index,"chebi"] != "character(0)")
    chebi_id = metaboliteSet[index[index.sub],"chebi"]
    
    chebi = BridgeDbR::getSystemCode("ChEBI")
    if (length(chebi_id)>0){
      for (k in 1:length(chebi_id)){
        if (!is.null(unlist(map(mapper, chebi, chebi_id[k], hmdb)[1]))) {
          metaboliteSet[index[index.sub[k]],"hmdb"] = unlist(map(mapper, chebi, chebi_id[k], hmdb)[1])
        }
      }
    }
    
    # Get rid of any rows that still don't have an hmdb code
    index = which(metaboliteSet[,"hmdb"] == "character(0)")
    
    if (length(index)>0) metaboliteSet <- metaboliteSet[-index,,drop=FALSE]
    
    # Get the total number of unique metabolites for all genes in one vector
    # nMets=c(nMets, length(unique(metaboliteSet[,"hmdb"])))
    # Improved check for same compounds, including if either chebi, kegg or hmdb have duplicated values. Return number of metabolites
    # nMets <- c(nMets, sum(apply(!apply(metaboliteSet[,c("hmdb","chebi","kegg")], 2, duplicated, incomparables = NA), 1, all)))
    
    # Remove duplicate metabolites (added by Marten)
    metaboliteSet <- metaboliteSet[apply(!apply(metaboliteSet[,c("hmdb","chebi","kegg")], 2, duplicated, incomparables = NA), 1, all),, drop = FALSE]
    nMets <- nrow(metaboliteSet)
    
    savepoint_p.values <- p.values
    savepoint_metaboliteSet <- metaboliteSet
    
    
    
    # Peform MSEA, one single patient and rare-gene variant at a time
    # 1 = Fishers exact; 3 = Fisher weighted
    retVal = performMSEA(metaboliteSet = savepoint_metaboliteSet, 
                         av_int_and_z_values_matrix = savepoint_p.values, 
                         patient = patients[i], 
                         gene_in, 
                         n_patients, 
                         thresh_F_pos, 
                         thresh_F_neg, 
                         path, 
                         test = 1, 
                         top, 
                         id, 
                         adductsSummed = FALSE
    )
    
    # Get only the p.value and put all p.values with the patient and gene in one table
    print(retVal)
    p_value = as.numeric(retVal$p.value)
    if (length(p_value)==0){
      p_value=NA
    }
    metSetResult = rbind(metSetResult, c("p.value"=p_value, "patient"=paste0("P", patients[i]), "metabolite.set"=gene_in))
    
  }
  
  # Create a short excel file with the p.values and the number of metabolites associated with a gene
  tmp <- data.frame("HGNC"=metSetResult[,3],"p.value"=as.numeric(metSetResult[,"p.value"]), "metabolites"=nMets)
  genExcelFileShort(tmp[order(tmp[,"p.value"]),], paste(path,"/P",patients[i],"/Recon2/MSEA_results.xls",sep=""))
  
  if (!is.null(genes)){
    tmp1 = tmp[order(tmp[,"p.value"]),]
    
    dummy = c(NA,0,0)
    for (l in 1:(100-dim(tmp)[1])){
      tmp1 = rbind(tmp1,dummy)
    }
    
    overview = rbind(overview, t(tmp1))
  }
  
}
# if (!is.null(genes)) genExcelFileShort(overview, paste(path,"/MSEA_overview.xls",sep=""))
Sys.time()