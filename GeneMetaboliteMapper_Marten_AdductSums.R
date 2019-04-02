# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Session info ------------------------------------------------------------
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

cat("
    Copied from 'GeneMetaboliteMapper_ME_Marten.R, on 29/03/2019 which in turn was copied from:
    'GeneMetaboliteMapper_ME.R' in the Metab/Metabolomics/DIMS_pipeline/R_workspace_ME/Crossomics/Crossomics_SinglePatients/src folder
    on 26/03/2019 to investigate and annotate the functionality of the crossomics pipeline.
    Disabled all save/generateExcel functions
    
    package versions:
    bioDist
    Cairo
    XLConnect
    BridgeDbR 
    ")

# added:
try(setwd("/Volumes/DATA/Metabolomics/DIMS_pipeline/R_workspace_ME/Crossomics/Crossomics_SinglePatients"), silent = TRUE)
try(setwd("/Volumes/DATA-1/Metabolomics/DIMS_pipeline/R_workspace_ME/Crossomics/Crossomics_SinglePatients"), silent = TRUE)
try(setwd("/Volumes/DATA-2/Metabolomics/DIMS_pipeline/R_workspace_ME/Crossomics/Crossomics_SinglePatients"), silent = TRUE)
try(setwd("/Volumes/DATA-3/Metabolomics/DIMS_pipeline/R_workspace_ME/Crossomics/Crossomics_SinglePatients"), silent = TRUE)
try(setwd("/Volumes/Metab/Metabolomics/DIMS_pipeline/R_workspace_ME/Crossomics/Crossomics_SinglePatients"), silent = TRUE)




# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Manual input ------------------------------------------------------------
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# To what distance from the reaction should metabolites be considered?
# distance_to_gene <- 1 # Staat nog niet aan, kijk naar het begin van de "then perform MSEA"

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Libraries ---------------------------------------------------------------
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# library("R.matlab")
# library("KEGGREST")
library("bioDist")
library("Cairo")
library("XLConnect")
library("BridgeDbR")

source("./src/Crossomics/sourceDir.R")
sourceDir("./src/Crossomics")


################################################################
## Load DIMS datas
################################################################

## Load direct combined Rdata object (obsolete?)
#load("./input_DIMS/CROSS-OMICS_P80,82.RData")


## Load DIMS input files
loadRData <- function(fileName){
  load(fileName)
  get(ls()[ls() != "fileName"])
}

# Added:
# outlist.neg.stats.more <- loadRData("../../HPCxls/results_RES_DBS_20180720_Run4_Diagnosis2017_1/outlist_identified_negative.RData")
# outlist.pos.stats.more <- loadRData("../../HPCxls/results_RES_DBS_20180720_Run4_Diagnosis2017_1/outlist_identified_positive.RData")
# To get the adducts_summed data
outlist.neg.adducts.HMDB <- loadRData("../../HPCxls/results_RES_DBS_20180720_Run4_Diagnosis2017_1/adductSums_negative.RData")
outlist.pos.adducts.HMDB <- loadRData("../../HPCxls/results_RES_DBS_20180720_Run4_Diagnosis2017_1/adductSums_positive.RData")

tmp <- intersect(colnames(outlist.neg.adducts.HMDB), colnames(outlist.pos.adducts.HMDB))
outlist.neg.adducts.HMDB <- outlist.neg.adducts.HMDB[,tmp]
outlist.pos.adducts.HMDB <- outlist.pos.adducts.HMDB[,tmp]
# outlist.neg.stats.more <- loadRData("Z:/Metabolomics/DIMS_pipeline/R_workspace_ME/HPCxls/results_RES_DBS_20180720_Run4_Diagnosis2017_1/outlist_identified_negative.RData")
# outlist.pos.stats.more <- loadRData("Z:/Metabolomics/DIMS_pipeline/R_workspace_ME/HPCxls/results_RES_DBS_20180720_Run4_Diagnosis2017_1/outlist_identified_positive.RData")

#WORKS outlist.neg.stats.more <- loadRData("Z:/Metabolomics/Research Metabolic Diagnostics/Metabolomics Projects/Projects 2015/Project 2015_011_SinglePatients/06 SinglePatients_VI/Bioinformatics 20180824/outlist_identified_negative.RData")
#WORKS outlist.pos.stats.more <- loadRData("Z:/Metabolomics/Research Metabolic Diagnostics/Metabolomics Projects/Projects 2015/Project 2015_011_SinglePatients/06 SinglePatients_VI/Bioinformatics 20180824/outlist_identified_positive.RData")


########################################################################################################
### Declare patient names 
#patients=c(80,82)
#controls=c(30:34,36:45)
# If WES data available genes=NULL, genes from WES in "./db/P28_HGNC.txt"
#genes=NULL
########################################################################################################

#patients=c(82)
#controls=c(30:34,36:45)
#genes=NULL

patients=c(221)
# ME: USE ALL CONTROLS in list!
# controls=c(63:65,69,73:77,79:81,83,88:89,96:99,102:104,107,109,111:112,114,117:118,120)
controls_list <- unique(strsplit(colnames(outlist.neg.adducts.HMDB)[grep("C", colnames(outlist.neg.adducts.HMDB))], "[C.]"))
for(i in c(1:length(controls_list))){
  if(i == 1){control_numbers <- as.integer(controls_list[[i]][2])}
  else{control_numbers <- c(control_numbers, as.integer(controls_list[[i]][2]))}
}
controls <- control_numbers[order(control_numbers)]
genes=NULL


n_patients=length(patients)
n_controls=length(controls)
########################################################################################################

###############################################################################################
################  Seperate modes only selected adducts ########################################
###############################################################################################


# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Collate the positive and negative adducts data --------------------------
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
index.neg <- which(rownames(outlist.neg.adducts.HMDB) %in% rownames(outlist.pos.adducts.HMDB))
index.pos <- which(rownames(outlist.pos.adducts.HMDB) %in% rownames(outlist.neg.adducts.HMDB))

tmp.pos <- outlist.pos.adducts.HMDB[rownames(outlist.pos.adducts.HMDB)[index.pos], 1:(dim(outlist.pos.adducts.HMDB)[2]-1)]
tmp.hmdb_name.pos <- outlist.pos.adducts.HMDB[rownames(outlist.pos.adducts.HMDB)[index.pos], dim(outlist.pos.adducts.HMDB)[2]]
tmp.pos.left <- outlist.pos.adducts.HMDB[-index.pos,]

tmp.neg <- outlist.neg.adducts.HMDB[rownames(outlist.pos.adducts.HMDB)[index.pos], 1:(dim(outlist.neg.adducts.HMDB)[2]-1)]
tmp.neg.left <- outlist.neg.adducts.HMDB[-index.neg,]

tmp <- apply(tmp.pos, 2,as.numeric) + apply(tmp.neg, 2,as.numeric)
rownames(tmp) <- rownames(tmp.pos)
tmp <- cbind(tmp, "HMDB_name"=tmp.hmdb_name.pos)
# adducts.neg.pos <- rbind(tmp, tmp.pos.left,tmp.neg.left) 
outlist.adducts.HMDB <- rbind(tmp, tmp.pos.left,tmp.neg.left) 

# dummy.neg <- rep(NA, dim(adducts.neg.pos)[1])
# outlist.adducts.HMDB <- cbind("mzmed.pgrp"=dummy.neg,
#                              "fq.best"=dummy.neg,
#                              "fq.worst"=dummy.neg,
#                              "nrsamples"=dummy.neg,
#                              "mzmin.pgrp"=dummy.neg,
#                              "mzmax.pgrp"=dummy.neg,
#                              adducts.neg.pos)

outlist.adducts.HMDB <- cbind(outlist.adducts.HMDB, "HMDB_code"=rownames(outlist.adducts.HMDB))




# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Get Z scores ------------------------------------------------------------
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

peaklist <- as.data.frame(outlist.adducts.HMDB)

# calculate mean and sd for Control group
ctrl.cols <- grep("C", colnames(peaklist),fixed = TRUE)
int.cols <- c(ctrl.cols, grep("P", colnames(peaklist),fixed = TRUE))


peaklist[,int.cols][peaklist[,int.cols]==0] <- NA

# tmp = data.matrix(peaklist[ , ctrl.cols], rownames.force = TRUE)
tmp = outlist.adducts.HMDB[ , ctrl.cols]

peaklist$avg.ctrls <- apply(tmp, 1, function(x) mean(as.numeric(x),na.rm = TRUE))
peaklist$sd.ctrls <- apply(tmp, 1, function(x) sd(as.numeric(x),na.rm = TRUE))

cnames.z = NULL

# calculate the actual Z scores
for (i in int.cols) {
  cname = colnames(peaklist)[i]
  cnames.z = c(cnames.z, paste(cname, "Zscore", sep="_"))
  zscores.1col <- (as.numeric(as.vector(unlist(peaklist[ , i]))) - peaklist$avg.ctrls) / peaklist$sd.ctrls
  peaklist <- cbind(peaklist, zscores.1col)
}

colnames(peaklist)[grep("zscores.1col", colnames(peaklist))[1]:ncol(peaklist)] <- cnames.z

# Put columns in desired order (intensity and zscore last)
peaklist <- peaklist[,c(grep("^[CP]\\d+", colnames(peaklist), value = TRUE, invert = TRUE),grep("^[CP]\\d+", colnames(peaklist), value = TRUE))]




# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Get mean Z-values for specific Patient ----------------------------------
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

adductsSummed <- TRUE

# Dit geeft geen p.values, maar de gemiddelde Z-score van de HMDB intensities van de patient. 
adductsHMDB_z_ave_and_int <- getPvalues(peaklist = peaklist, 
                               n_patients, 
                               n_controls, 
                               assi.lab = "HMDB_code",
                               patients,
                               controls,
                               adducts = FALSE)




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
path="./results/crossomics_Fisher_weighted"
path <- "/Users/mkerkho7/DIMS2_repo/TestResults/"
# dir.create(path, showWarnings = FALSE)
thresh_F_pos=1.5
thresh_F_neg=-1
top=20
# id="InChI_key"
id="hmdb"

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

path2 = "../Crossomics_Build_Metabolite_Set/results/metabolite_sets_step_0,1,2,3_1.0_filter_1.1/mss_1"
# path2 = paste0("../Crossomics_Build_Metabolite_Set/results/metabolite_sets_step_0,1,2,3_1.0_filter_1.1/mss_",distance_to_gene)
step = 1 # 1:3  # ME wordt gebruikt in regel 495 -503, waarom?
path_wg = "./results/mss_WG_step_0" # for sampling random genes   #ME enkel voor random data set, niet WES
sets = list.files(path = path2)
overview = NULL # ME ?
rank = 0  # ME wordt niet gebruikt
p_value_sum = 0 # ME wordt niet gebruikt

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
  nMets = NULL    # ME list with number of metabolites per gene
  for (j in 1:length(mss)){
    if (!file.exists(paste(path2, mss[j], sep="/"))) next
    
    load(paste(path2, mss[j], sep="/"))
    gene_in=strsplit(mss[j], split = "." , fixed=TRUE)[[1]][1]
    
    # A lot of missing HMDB identifiers!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    # metaboliteSet = retVal$mets
    if (step==2){
      metaboliteSet = result_mets_2
    } else if (step==1){
      metaboliteSet = result_mets_1
    } else if (step==3){
      metaboliteSet = result_mets_3
    } else {
      metaboliteSet = result_mets_0
    }
    
    
    if (!is.null(metaboliteSet)){
      if (is.null(dim(metaboliteSet))){
        metaboliteSet=data.frame(t(metaboliteSet))  # ME matrix transpose
      }
    }
    
    ##################################################
    # temporarily work around to be fixed in findMetabolicEnvironment
    index = which(is.na(metaboliteSet[,"hmdb"]))  
    if (length(index)>0) metaboliteSet[index,"hmdb"] = "character(0)"  # ME als NA in HMDB colom van resul_mets_[n], zet op unknown
    ##################################################
    
    ########### Recon 2.0 ############################
    index = which(metaboliteSet[,"hmdb"] == "character(0)") # ME list positions of "Unkown" in hmdb column
    
    # present in BridgeDB?
    index.sub = which(metaboliteSet[index,"kegg"] != "character(0)") # ME Welke unkwon hmdb zijn niet unknown in KEGG
    kegg_id = metaboliteSet[index[index.sub],"kegg"]
    
    # ADDED (new .bridge file)
    mapper <- loadDatabase("/Users/mkerkho7/DIMS2_repo/Crossomics/metabolites_20190207.bridge")
    # mapper = loadDatabase("./db/BridgeDB/metabolites_20150717.bridge")  # ME loadDatabase function zit in BridgeDbR library 
    hmdb = getSystemCode("HMDB") # ME function zit in BridgeDbR library
    kegg = getSystemCode("KEGG Compound")# ME function zit in BridgeDbR library
    
    if (length(kegg_id)>0){
      for (k in 1:length(kegg_id)){
        if (!is.null(unlist(map(mapper, kegg, kegg_id[k], hmdb)[1]))) metaboliteSet[index[index.sub[k]],"hmdb"] = unlist(map(mapper, kegg, kegg_id[k], hmdb)[1])
      }
    }
    
    index = which(metaboliteSet[,"hmdb"] == "character(0)")
    index.sub = which(metaboliteSet[index,"chebi"] != "character(0)")
    chebi_id = metaboliteSet[index[index.sub],"chebi"]
    
    chebi = getSystemCode("ChEBI")
    if (length(chebi_id)>0){
      for (k in 1:length(chebi_id)){
        if (!is.null(unlist(map(mapper, chebi, chebi_id[k], hmdb)[1]))) metaboliteSet[index[index.sub[k]],"hmdb"] = unlist(map(mapper, chebi, chebi_id[k], hmdb)[1])
      }
    }
    
    index = which(metaboliteSet[,"hmdb"] == "character(0)")
    
    if (length(index)>0) metaboliteSet = metaboliteSet[-index,,drop=FALSE]
    
    # nMets=c(nMets, dim(metaboliteSet)[1])
    nMets=c(nMets, length(unique(metaboliteSet[,"hmdb"])))
    
    savepoint_p.values <- p.values
    
    
    
    
    # 3 = Fisher weighted
    retVal = performMSEA(metaboliteSet = metaboliteSet, 
                         p_valuesAll = p.values, 
                         patient = patients[i], 
                         gene_in, 
                         n_patients, 
                         thresh_F_pos, 
                         thresh_F_neg, 
                         path, 
                         test = 3, 
                         top, 
                         id, 
                         adductsSummed = FALSE
    )
    print(retVal)
    p_value = as.numeric(retVal$p.value)
    if (length(p_value)==0){
      p_value=NA
    }
    
    metSetResult = rbind(metSetResult, c("p.value"=p_value, "patient"=paste("P", patients[i], sep=""), "metabolite.set"=gene_in))
    
  }
  
  tmp=data.frame("HGNC"=metSetResult[,3],"p.value"=as.numeric(metSetResult[,"p.value"]), "metabolites"=nMets)
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


### End of Pipeline for WES #####











########################################################################
########################################################################
########################################################################
path="./results/crossomics_Fisher_Weighted"
# dir.create(path, showWarnings = FALSE)
overview = NULL

for (i in 1:length(patients)){
  
  ######################## MSEA Revon2 neighbourhood ######################################################
  # dir.create(paste(path,"/P",patients[i], sep=""), showWarnings = FALSE)
  # dir.create(paste(path,"/P",patients[i],"/Recon2", sep=""), showWarnings = FALSE)
  
  # message(paste("Patient", patients[i]))
  
  # subset to 1 patient!
  p.values=p.values.assi.all[,c(grep(paste("P",patients[i],sep=""), colnames(p.values.assi.all), fixed=TRUE),
                                grep("C", colnames(p.values.assi.all), fixed=TRUE),
                                which(colnames(p.values.assi.all)=="avg.int.controls"))] 
  
  load(paste0("./results/mss_",genes[i][[1]][1],i,".RData"))
  
  metSetResult = NULL
  nMets = NULL
  
  #   for (j in 1:dim(gene_map_table)[1]){
  for (j in 1:length(mss)){
    
    if (!file.exists(paste(path2, mss[j], sep="/"))) next
    
    load(paste(path2, mss[j], sep="/"))
    gene_in=strsplit(mss[j], split = "." , fixed=TRUE)[[1]][1]
    
    # A lot of missing HMDB identifiers!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    # metaboliteSet = retVal$mets
    if (step==2){
      metaboliteSet = result_mets_2
    } else if (step==1){
      metaboliteSet = result_mets_1
    } else if (step==3){
      metaboliteSet = result_mets_3
    } else {
      metaboliteSet = result_mets_0
    }
    
    if (!is.null(metaboliteSet)){
      if (is.null(dim(metaboliteSet))){
        metaboliteSet=data.frame(t(metaboliteSet))
      }
    }  
    
    ##################################################
    # temporarily work around to be fixed in findMetabolicEnvironment
    index = which(is.na(metaboliteSet[,"hmdb"]))
    if (length(index)>0) metaboliteSet[index,"hmdb"] = "character(0)"
    ##################################################
    
    ########### Recon 2.0 ############################
    index = which(metaboliteSet[,"hmdb"] == "character(0)") 
    
    # pressent in BridgeDB?
    index.sub = which(metaboliteSet[index,"kegg"] != "character(0)")  # ME filter metaboliteSet for empty Kegg 
    kegg_id = metaboliteSet[index[index.sub],"kegg"] # ME get kegg_id's for metaboliteSet
    
    mapper = loadDatabase("./db/BridgeDB/metabolites_20150717.bridge")
    hmdb = getSystemCode("HMDB")  
    kegg = getSystemCode("KEGG Compound") 
    
    if (length(kegg_id)>0){
      for (k in 1:length(kegg_id)){
        if (!is.null(unlist(map(mapper, kegg, kegg_id[k], hmdb)[1]))) metaboliteSet[index[index.sub[k]],"hmdb"] = unlist(map(mapper, kegg, kegg_id[k], hmdb)[1]) 
      }
    }  
    
    index = which(metaboliteSet[,"hmdb"] == "character(0)")
    index.sub = which(metaboliteSet[index,"chebi"] != "character(0)")
    chebi_id = metaboliteSet[index[index.sub],"chebi"]
    
    chebi = getSystemCode("ChEBI")# ME system code for chebi ??
    
    if (length(chebi_id)>0){
      for (k in 1:length(chebi_id)){
        if (!is.null(unlist(map(mapper, chebi, chebi_id[k], hmdb)[1]))) metaboliteSet[index[index.sub[k]],"hmdb"] = unlist(map(mapper, chebi, chebi_id[k], hmdb)[1]) 
      }
    }
    
    index = which(metaboliteSet[,"hmdb"] == "character(0)") 
    if (length(index)>0) metaboliteSet = metaboliteSet[-index,,drop=FALSE]
    #######################################################
    
    #     ########### Recon 2.2 ############################
    #     inchi = metaboliteSet[,"InChI_key"]
    #     chebi_id = metaboliteSet[,"chebi"]
    #     index=which((is.na(chebi_id)&is.na(inchi)))
    #     if (length(index)>0){
    #       metaboliteSet = metaboliteSet[-index,]
    #     }
    #     #######################################################
    
    # Add glycine and carnitine conjugates to CoA compounds
    #addConjugates(metaboliteSet)
    
    #     if (!is.null(metaboliteSet) && (dim(metaboliteSet)[1]!=0)){
    #       if (is.null(dim(metaboliteSet))){
    #         metaboliteSet=data.frame(t(metaboliteSet))
    #       }
    
    # nMets=c(nMets, dim(metaboliteSet)[1])
    nMets=c(nMets, length(unique(metaboliteSet[,"hmdb"])))
    # message(paste("gene: ", gene_in))
    # message(paste("metaboliteSet: ", dim(metaboliteSet)[1]))
    
    #     Warning messages:
    #       1: In data.row.names(row.names, rowsi, i) :
    #       some row.names duplicated: 6 --> row.names NOT used
    #     2: In data.row.names(row.names, rowsi, i) :
    #       some row.names duplicated: 17 --> row.names NOT used
    
    
    # which(is.na(p.values[1,]))
    
    retVal = performMSEA(metaboliteSet, p.values, patients[i], gene_in, n_patients, thresh_F_pos, thresh_F_neg, path, 3, top, id, adductsSummed)
    
    p_value = as.numeric(retVal$p.value)
    if (length(p_value)==0){
      p_value=NA
    }
    
    metSetResult = rbind(metSetResult, c("p.value"=p_value, "patient"=paste("P", patients[i], sep=""), "metabolite.set"=gene_in))
    
  }
  
  tmp=data.frame("HGNC"=metSetResult[,3],"p.value"=as.numeric(metSetResult[,"p.value"]), "metabolites"=nMets)
  # genExcelFileShort(tmp[order(as.numeric(tmp[,"p.value"])),], paste(path,"/P",patients[i],"/Recon2/MSEA_results.xls",sep=""))
  
  if (!is.null(genes)){
    tmp1 = tmp[order(tmp[,"p.value"]),]
    
    dummy = c(NA,0,0)
    for (l in 1:(100-dim(tmp)[1])){
      tmp1 = rbind(tmp1,dummy) 
    }
    
    overview = rbind(overview, t(tmp1))
  }
  
  # # rank
  # index.rank = min(which(tmp2[1,]%in%unlist(genes[i])))
  # rank = rank + index.rank 
  # p_value_sum = p_value_sum + as.numeric(tmp2[2,index.rank]) 
  ##########################################################################################################
}
# if (!is.null(genes)) genExcelFileShort(overview, paste(path,"/MSEA_overview.xls",sep=""))
# rank = rank/length(patients)
# save(rank, p_value_sum, file=paste(path,"/rank.RData",sep=""))
Sys.time()

########################################################################
########################################################################
########################################################################
path="./results/crossomics_Hyper"
# dir.create(path, showWarnings = FALSE)
overview = NULL

for (i in 1:length(patients)){
  
  ######################## MSEA Revon2 neighbourhood ######################################################
  # dir.create(paste(path,"/P",patients[i], sep=""), showWarnings = FALSE)
  # dir.create(paste(path,"/P",patients[i],"/Recon2", sep=""), showWarnings = FALSE)
  
  # message(paste("Patient", patients[i]))
  
  # subset to 1 patient!
  p.values=p.values.assi.all[,c(grep(paste("P",patients[i],sep=""), colnames(p.values.assi.all), fixed=TRUE),
                                grep("C", colnames(p.values.assi.all), fixed=TRUE),
                                which(colnames(p.values.assi.all)=="avg.int.controls"))] 
  
  load(paste0("./results/mss_",genes[i],i,".RData"))
  
  metSetResult = NULL
  nMets = NULL
  
  #   for (j in 1:dim(gene_map_table)[1]){
  for (j in 1:length(mss)){
    
    if (!file.exists(paste(path2, mss[j], sep="/"))) next
    
    load(paste(path2, mss[j], sep="/"))
    gene_in=strsplit(mss[j], split = "." , fixed=TRUE)[[1]][1]
    
    # A lot of missing HMDB identifiers!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    # metaboliteSet = retVal$mets
    if (step==2){
      metaboliteSet = result_mets_2
    } else if (step==1){
      metaboliteSet = result_mets_1
    } else if (step==3){
      metaboliteSet = result_mets_3
    } else {
      metaboliteSet = result_mets_0
    }
    
    if (!is.null(metaboliteSet)){
      if (is.null(dim(metaboliteSet))){
        metaboliteSet=data.frame(t(metaboliteSet))
      }
    }  
    
    ##################################################
    # temporarily work around to be fixed in findMetabolicEnvironment
    index = which(is.na(metaboliteSet[,"hmdb"]))
    if (length(index)>0) metaboliteSet[index,"hmdb"] = "character(0)"
    ##################################################
    
    ########### Recon 2.0 ############################
    index = which(metaboliteSet[,"hmdb"] == "character(0)") 
    
    # pressent in BridgeDB?
    index.sub = which(metaboliteSet[index,"kegg"] != "character(0)")
    kegg_id = metaboliteSet[index[index.sub],"kegg"]
    
    mapper = loadDatabase("./db/BridgeDB/metabolites_20150717.bridge")
    hmdb = getSystemCode("HMDB")
    kegg = getSystemCode("KEGG Compound")
    
    if (length(kegg_id)>0){
      for (k in 1:length(kegg_id)){
        if (!is.null(unlist(map(mapper, kegg, kegg_id[k], hmdb)[1]))) metaboliteSet[index[index.sub[k]],"hmdb"] = unlist(map(mapper, kegg, kegg_id[k], hmdb)[1]) 
      }
    }  
    
    index = which(metaboliteSet[,"hmdb"] == "character(0)")
    index.sub = which(metaboliteSet[index,"chebi"] != "character(0)")
    chebi_id = metaboliteSet[index[index.sub],"chebi"]
    
    chebi = getSystemCode("ChEBI")
    
    if (length(chebi_id)>0){
      for (k in 1:length(chebi_id)){
        if (!is.null(unlist(map(mapper, chebi, chebi_id[k], hmdb)[1]))) metaboliteSet[index[index.sub[k]],"hmdb"] = unlist(map(mapper, chebi, chebi_id[k], hmdb)[1]) 
      }
    }
    
    index = which(metaboliteSet[,"hmdb"] == "character(0)") 
    if (length(index)>0) metaboliteSet = metaboliteSet[-index,,drop=FALSE]
    #######################################################
    
    #     ########### Recon 2.2 ############################
    #     inchi = metaboliteSet[,"InChI_key"]
    #     chebi_id = metaboliteSet[,"chebi"]
    #     index=which((is.na(chebi_id)&is.na(inchi)))
    #     if (length(index)>0){
    #       metaboliteSet = metaboliteSet[-index,]
    #     }
    #     #######################################################
    
    # Add glycine and carnitine conjugates to CoA compounds
    #addConjugates(metaboliteSet)
    
    #     if (!is.null(metaboliteSet) && (dim(metaboliteSet)[1]!=0)){
    #       if (is.null(dim(metaboliteSet))){
    #         metaboliteSet=data.frame(t(metaboliteSet))
    #       }
    
    # nMets=c(nMets, dim(metaboliteSet)[1])
    nMets=c(nMets, length(unique(metaboliteSet[,"hmdb"])))
    # message(paste("gene: ", gene_in))
    # message(paste("metaboliteSet: ", dim(metaboliteSet)[1]))
    
    #     Warning messages:
    #       1: In data.row.names(row.names, rowsi, i) :
    #       some row.names duplicated: 6 --> row.names NOT used
    #     2: In data.row.names(row.names, rowsi, i) :
    #       some row.names duplicated: 17 --> row.names NOT used
    
    
    # which(is.na(p.values[1,]))
    
    retVal = performMSEAenv(metaboliteSet, p.values, patients[i], gene_in, n_patients, thresh_F_pos, thresh_F_neg, path, 1, top, id, adductsSummed)
    
    p_value = as.numeric(retVal$p.value)
    if (length(p_value)==0){
      p_value=NA
    }
    
    metSetResult = rbind(metSetResult, c("p.value"=p_value, "patient"=paste("P", patients[i], sep=""), "metabolite.set"=gene_in))
    
  }
  
  tmp=data.frame("HGNC"=metSetResult[,3],"p.value"=as.numeric(metSetResult[,"p.value"]), "metabolites"=nMets)
  # genExcelFileShort(tmp[order(as.numeric(tmp[,"p.value"])),], paste(path,"/P",patients[i],"/Recon2/MSEA_results.xls",sep=""))
  
  if (!is.null(genes)){
    tmp1 = tmp[order(tmp[,"p.value"]),]
    
    dummy = c(NA,0,0)
    for (l in 1:(100-dim(tmp)[1])){
      tmp1 = rbind(tmp1,dummy) 
    }
    
    overview = rbind(overview, t(tmp1))
  }
  
  # # rank
  # index.rank = min(which(tmp2[1,]%in%unlist(genes[i])))
  # rank = rank + index.rank 
  # p_value_sum = p_value_sum + as.numeric(tmp2[2,index.rank]) 
  ##########################################################################################################
}
# if (!is.null(genes)) genExcelFileShort(overview, paste(path,"/MSEA_overview.xls",sep=""))
# rank = rank/length(patients)
# save(rank, p_value_sum, file=paste(path,"/rank.RData",sep=""))
Sys.time()

########################################################################
########################################################################
########################################################################
path="./results/crossomics_Hyper_Weighted"
# dir.create(path, showWarnings = FALSE)
overview = NULL

for (i in 1:length(patients)){
  
  ######################## MSEA Revon2 neighbourhood ######################################################
  # dir.create(paste(path,"/P",patients[i], sep=""), showWarnings = FALSE)
  # dir.create(paste(path,"/P",patients[i],"/Recon2", sep=""), showWarnings = FALSE)
  
  # message(paste("Patient", patients[i]))
  
  # subset to 1 patient!
  p.values=p.values.assi.all[,c(grep(paste("P",patients[i],sep=""), colnames(p.values.assi.all), fixed=TRUE),
                                grep("C", colnames(p.values.assi.all), fixed=TRUE),
                                which(colnames(p.values.assi.all)=="avg.int.controls"))] 
  
  load(paste0("./results/mss_",genes[i],i,".RData"))
  
  metSetResult = NULL
  nMets = NULL
  
  #   for (j in 1:dim(gene_map_table)[1]){
  for (j in 1:length(mss)){
    
    if (!file.exists(paste(path2, mss[j], sep="/"))) next
    
    load(paste(path2, mss[j], sep="/"))
    gene_in=strsplit(mss[j], split = "." , fixed=TRUE)[[1]][1]
    
    # A lot of missing HMDB identifiers!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    # metaboliteSet = retVal$mets
    if (step==2){
      metaboliteSet = result_mets_2
    } else if (step==1){
      metaboliteSet = result_mets_1
    } else if (step==3){
      metaboliteSet = result_mets_3
    } else {
      metaboliteSet = result_mets_0
    }
    
    if (!is.null(metaboliteSet)){
      if (is.null(dim(metaboliteSet))){
        metaboliteSet=data.frame(t(metaboliteSet))
      }
    }  
    
    ##################################################
    # temporarily work around to be fixed in findMetabolicEnvironment
    index = which(is.na(metaboliteSet[,"hmdb"]))
    if (length(index)>0) metaboliteSet[index,"hmdb"] = "character(0)"
    ##################################################
    
    ########### Recon 2.0 ############################
    index = which(metaboliteSet[,"hmdb"] == "character(0)") 
    
    # pressent in BridgeDB?
    index.sub = which(metaboliteSet[index,"kegg"] != "character(0)")
    kegg_id = metaboliteSet[index[index.sub],"kegg"]
    
    mapper = loadDatabase("./db/BridgeDB/metabolites_20150717.bridge")
    hmdb = getSystemCode("HMDB")
    kegg = getSystemCode("KEGG Compound")
    
    if (length(kegg_id)>0){
      for (k in 1:length(kegg_id)){
        if (!is.null(unlist(map(mapper, kegg, kegg_id[k], hmdb)[1]))) metaboliteSet[index[index.sub[k]],"hmdb"] = unlist(map(mapper, kegg, kegg_id[k], hmdb)[1]) 
      }
    }  
    
    index = which(metaboliteSet[,"hmdb"] == "character(0)")
    index.sub = which(metaboliteSet[index,"chebi"] != "character(0)")
    chebi_id = metaboliteSet[index[index.sub],"chebi"]
    
    chebi = getSystemCode("ChEBI")
    
    if (length(chebi_id)>0){
      for (k in 1:length(chebi_id)){
        if (!is.null(unlist(map(mapper, chebi, chebi_id[k], hmdb)[1]))) metaboliteSet[index[index.sub[k]],"hmdb"] = unlist(map(mapper, chebi, chebi_id[k], hmdb)[1]) 
      }
    }
    
    index = which(metaboliteSet[,"hmdb"] == "character(0)") 
    if (length(index)>0) metaboliteSet = metaboliteSet[-index,,drop=FALSE]
    #######################################################
    
    #     ########### Recon 2.2 ############################
    #     inchi = metaboliteSet[,"InChI_key"]
    #     chebi_id = metaboliteSet[,"chebi"]
    #     index=which((is.na(chebi_id)&is.na(inchi)))
    #     if (length(index)>0){
    #       metaboliteSet = metaboliteSet[-index,]
    #     }
    #     #######################################################
    
    # Add glycine and carnitine conjugates to CoA compounds
    #addConjugates(metaboliteSet)
    
    #     if (!is.null(metaboliteSet) && (dim(metaboliteSet)[1]!=0)){
    #       if (is.null(dim(metaboliteSet))){
    #         metaboliteSet=data.frame(t(metaboliteSet))
    #       }
    
    # nMets=c(nMets, dim(metaboliteSet)[1])
    nMets=c(nMets, length(unique(metaboliteSet[,"hmdb"])))
    # message(paste("gene: ", gene_in))
    # message(paste("metaboliteSet: ", dim(metaboliteSet)[1]))
    
    #     Warning messages:
    #       1: In data.row.names(row.names, rowsi, i) :
    #       some row.names duplicated: 6 --> row.names NOT used
    #     2: In data.row.names(row.names, rowsi, i) :
    #       some row.names duplicated: 17 --> row.names NOT used
    
    
    # which(is.na(p.values[1,]))
    
    retVal = performMSEAenv(metaboliteSet, p.values, patients[i], gene_in, n_patients, thresh_F_pos, thresh_F_neg, path, 2, top, id, adductsSummed)
    
    p_value = as.numeric(retVal$p.value)
    if (length(p_value)==0){
      p_value=NA
    }
    
    metSetResult = rbind(metSetResult, c("p.value"=p_value, "patient"=paste("P", patients[i], sep=""), "metabolite.set"=gene_in))
    
  }
  
  tmp=data.frame("HGNC"=metSetResult[,3],"p.value"=as.numeric(metSetResult[,"p.value"]), "metabolites"=nMets)
  # genExcelFileShort(tmp[order(tmp[,"p.value"]),], paste(path,"/P",patients[i],"/Recon2/MSEA_results.xls",sep=""))
  
  if (!is.null(genes)){
    tmp1 = tmp[order(tmp[,"p.value"]),]
    
    dummy = c(NA,0,0)
    for (l in 1:(100-dim(tmp)[1])){
      tmp1 = rbind(tmp1,dummy) 
    }
    
    overview = rbind(overview, t(tmp1))
  }
  
  # # rank
  # index.rank = min(which(tmp2[1,]%in%unlist(genes[i])))
  # rank = rank + index.rank 
  # p_value_sum = p_value_sum + as.numeric(tmp2[2,index.rank]) 
  ##########################################################################################################
}
# if (!is.null(genes)) genExcelFileShort(overview, paste(path,"/MSEA_overview.xls",sep=""))
# rank = rank/length(patients)
# save(rank, p_value_sum, file=paste(path,"/rank.RData",sep=""))
Sys.time()

