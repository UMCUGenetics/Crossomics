# load("./input_DIMS/SinglePatients_IV_DIMS_2.1_4Crossomics.RData")
# load("./input_DIMS/SinglePatients_IV_DIMS_2.1_area_4Crossomics.RData")
# load("./input_DIMS/SinglePatients_IV_DIMS_2.1_filtered_4Crossomics.RData")
# load("./input_DIMS/BSP20160902_2,28,34-47.RData")

############### used for manuscript ############################
# load("./input_DIMS/SinglePatients_IV_DIMS_2.1_more_noise_4Crossomics.RData")
# load("./input_DIMS/SinglePatients_PAH_CBS_DIMS_2.1_more_noise_4Crossomics.RData")
# load("./input_DIMS/SinglePatients_GCDH_IVD_DIMS_2.1_more_noise_4Crossomics.RData")
# load("./input_DIMS/SinglePatients_MAT1A_MTHFR_DIMS_2.1_more_noise_4Crossomics.RData")
# load("./input_DIMS/SinglePatients_MCT1_normalized_DIMS_2.1_more_noise_4Crossomics.RData")
################################################################

# load("./input_DIMS/CROSS-OMICS_P5.RData")
# load("./input_DIMS/CROSS-OMICS_P38.RData")
load("./input_DIMS/CROSS-OMICS_P80,82.RData")
#load("./input_DIMS/CROSS-OMICSP_P15_pipeline_2.1.RData")
# load("./input_DIMS/RES_BSP_SP_IX.RData")
# load("./input_DIMS/RES_BSP_20171109_SP_X.RData")
# load("./input_DIMS/CROSS-OMICS_SP_VIII_P99.RData")
# load("./input_DIMS/CROSS-OMICS_SPVII_P95.RData")


#load()
#load()

# library("R.matlab")
# library("KEGGREST")
library("bioDist")
library("Cairo")
library("XLConnect")
library("BridgeDbR")

source("./src/Crossomics/sourceDir.R")
sourceDir("./src/Crossomics")

# Recon 2.2 metabolites have been annotated with InChI string representations of molecular structure.

# library(rsbml)
# sbml=rsbml_read("./db/MODEL1603150001.xml", consistency = FALSE)
# rdf = as.vector(unlist(sapply(species(model(sbml)), annotation)))
# length(rdf) # 4326
# 
# length(sbml@model@species) # 6047
# 
# unlist(strsplit(unlist(strsplit(rdf, "CHEBI:", fixed = TRUE))[2],"\"/>\n",fixed = TRUE))[1]

# model=readSBMLmod("./db/MODEL1603150001.xml")
# save(model, file="./db/Recon_2.2_biomodels.RData")
# load("./db/Recon_2.2_biomodels.RData")

#ls("package:sybil")
#recon2=readMat("./db/Recon2.v04.mat")
#save(recon2, file="./db/Recon2.RData")

# On HPC
######################################################################
# load("./db/Recon2.RData")
# print(recon2)
# model=recon2$modelR204[,,1]
# names(model)
######################################################################

# length(model$metCHEBIID) # 5063

# S=model@S
# dim(S) # 5324 x 7785

# head(model@react_id)
# head(model@met_id)
# 
# #  7785(rxns) X 1675(genes)  7440(rxns) X 2140(genes)
# rxnGeneMat=model@rxnGeneMat 
# dim(rxnGeneMat)

# ############################# get P-values #############################################################
# # HEXA
# patients=c(28,41:45)
# controls=c(1:17)
# # If WES data available genes=NULL, genes from WES in "./db/P28_HGNC.txt"
# genes=NULL
# ########################################################################################################

# # P5
# patients=c(5)
# controls=c(56:63)
# # If WES data available genes=NULL, genes from WES in "./db/P28_HGNC.txt"
# genes=NULL
# #########################################################################################################

# # P38
# patients=c(38)
# controls=c(1:20)
# # If WES data available genes=NULL, genes from WES in "./db/P28_HGNC.txt"
# genes=NULL
# ########################################################################################################

# P80,82
patients=c(80,82)
controls=c(30:34,36:45)
# If WES data available genes=NULL, genes from WES in "./db/P28_HGNC.txt"
genes=NULL
########################################################################################################

# # P15
# patients=c(15)
# controls=c(45:50)
# # If WES data available genes=NULL, genes from WES in "./db/P28_HGNC.txt"
# genes=NULL
# #########################################################################################################

# # P12
# patients=c(12)
# controls=c(1:20)
# # If WES data available genes=NULL, genes from WES in "./db/P28_HGNC.txt"
# genes=NULL
# #########################################################################################################

# # P137, P135, P138
# patients=c(137,135,138)
# controls=c(31:34,36:63)
# # If WES data available genes=NULL, genes from WES in "./db/P28_HGNC.txt"
# genes=NULL
# #########################################################################################################

# # P99 (18)
# patients=c(18)
# controls=c(1:30)
# # If WES data available genes=NULL, genes from WES in "./db/P28_HGNC.txt"
# genes=NULL
# #########################################################################################################

# # P95
# patients=c(95)
# controls=c(1:30)
# # If WES data available genes=NULL, genes from WES in "./db/P28_HGNC.txt"
# genes=NULL
# #########################################################################################################

# patients=c(3,17,18,21,23)
# controls=c(30,31,32,33,34,38,39,40,41,42,43,44,45,46,47,48,49,50)
# patients=c(4)
# controls=c(56,57,58,59)
# patients=c(28) #patients=c(28,41:45) #patients=c(2,28:47)
# controls=c(1:20)
# patients=c(14)
# controls=c(56:63)
# patients=c(5) #patients=c(2,28:47)
# controls=c(56:63)

# # SPIV
# # patients=c(53:71)
# patients=c(53:68,70,71)
# controls=c(30:34,36:45)
# genes = list(c("MCCC1","MCCC2"),
# c("GLDC","GCSH","AMT"),
# c("GLDC","GCSH","AMT"),
# "CYP27A1",
# "CYP27A1",
# "DHCR7",
# "DHCR7",
# "FAH",
# "ACADVL",
# "ACADM",
# "ACADM",
# "MUT",
# "MUT",
# "OTC",
# "OTC",
# c("ABCG5","ABCG8"),
# #"AGA",
# "SLC7A7",
# "CPT1A")

# # SP_PAH_CBS
# patients=c(1:5,8)
# controls=c(1:5)
# genes = c(rep("PAH",5),"CBS")

# # SP_GCDH_IVD
# patients=c(11:14)
# controls=c(21:34)
# genes = c(rep("GCDH",2),rep("IVD",2))

# # SP_MAT1A_MTHFR
# patients=c(6,7,9,10)
# controls=c(6:11)
# genes = c(rep("MAT1A",2),rep("MTHFR",2))

# # SP_MCT1
# patients=c(15)
# # remove outliers
# controls=c(1:40)
# # remove outliers
# controls=controls[-c(6,9,14,25,40)] # only if there is only 1 measurement!!!
# genes = c(rep("SLC16A1",1))

###################################################
n_patients=length(patients)
n_controls=length(controls)

# labelAdducts <- function(peaklist, adducts, adducts_long) {
#   
#   names = lapply(rownames(peaklist), function(x, adducts, adducts_long) {
#     # x=rownames(p.values.assi.pos)[1]
#     compounds=unlist(strsplit(as.vector(x),";"))
#     compounds=compounds[-1]
#     index=c(1:length(compounds))
#     adducts_index = grep("_", compounds, fixed=TRUE)
#     
#     if (length(adducts_index)>0) index=index[-adducts_index]
#     
#     all_adducts=compounds[adducts_index]
#     keep=NULL
#     for (i in 1:length(all_adducts)){
#       # unlist(strsplit(as.vector(all_adducts[i]),"_"))[2] %in% adducts
#       keep = c(keep, unlist(strsplit(as.vector(all_adducts[i]),"_"))[2] %in% adducts)
#     }
#     all_adducts=all_adducts[keep]
#     
#     for (i in 1:length(adducts)){
#       all_adducts=gsub(paste("_",toString(adducts[i]),sep=""), paste(" ", adducts_long[i], sep=""), all_adducts)
#     }
#     
#     compounds=compounds[index]
#     if (length(all_adducts)!=0){
#       if (length(compounds)!=0){
#         compounds=paste(compounds, all_adducts, sep=";", collapse=";")
#       } else {
#         compounds=paste(all_adducts, collapse=";")
#       }  
#     } else if (length(compounds)>1) {
#       compounds=paste(compounds, collapse=";")
#     }
#     
#     if (length(compounds)==0) compounds="" 
#     
#     return(compounds)
#   })
#   
#   index=which(names=="")
#   if (length(index)>0) names=names[-index]
#   peaklist.filt = peaklist[-index,] 
#   rownames(peaklist.filt) = names
#   
#   return(peaklist.filt)
# }

# # Z-scores to P-values
# if warnings => getPvalues labels P4 vs. P04


# index = grep("HMDB00062", outlist.pos.stats.more[,"HMDB_code"], fixed = TRUE)
# outlist.pos.stats.more[index,c("mzmed.pgrp","P71.1","P71.2","HMDB_code")]
# 
# index=grep("HMDB00062", rownames(p.values.assi.pos), fixed = TRUE)
# p.values.assi.pos[index, c("P71")]
# 
# p.values.assi.pos=p.values.assi.pos[index,]

###############################################################################################
################  Saperate modes only selected adducts ########################################
###############################################################################################
adductsSummed=FALSE

p.values.assi.pos = getPvalues(outlist.pos.stats.more, n_patients, n_controls, assi.lab="HMDB_code",patients,controls) # assi.hmdb, HMDB_code For new identification !!!!!!!!!!!!!!!!!!!!!!!!!!

# Filter adducts
names = lapply(rownames(p.values.assi.pos), function(x, adducts=c(1,2), adducts_long=c("[M+Na]+","[M+K]+")) {
  # x=rownames(p.values.assi.pos)[1]
  compounds=unlist(strsplit(as.vector(x),";"))
  # compounds=compounds[-1]  # since DIMS 2.1
  index=c(1:length(compounds))
  adducts_index = grep("_", compounds, fixed=TRUE)
  
  if (length(adducts_index)>0) index=index[-adducts_index]
  
  all_adducts=compounds[adducts_index]
  keep=NULL
  for (i in 1:length(all_adducts)){
    # unlist(strsplit(as.vector(all_adducts[i]),"_"))[2] %in% adducts
    keep = c(keep, unlist(strsplit(as.vector(all_adducts[i]),"_"))[2] %in% adducts)
  }
  all_adducts=all_adducts[keep]
  
  for (i in 1:length(adducts)){
    all_adducts=gsub(paste("_",toString(adducts[i]),sep=""), paste(" ", adducts_long[i], sep=""), all_adducts)
  }
  
  compounds=compounds[index]
  if (length(all_adducts)!=0){
    if (length(compounds)!=0){
      compounds=paste(compounds, all_adducts, sep=";", collapse=";")
    } else {
      compounds=paste(all_adducts, collapse=";")
    }  
  } else if (length(compounds)>1) {
    compounds=paste(compounds, collapse=";")
  }
  
  if (length(compounds)==0) compounds="" 
  
  return(compounds)
})

index=which(names=="")
if (length(index)>0) names=names[-index]
p.values.assi.pos.filt = p.values.assi.pos[-index,] 
rownames(p.values.assi.pos.filt) = names

save(p.values.assi.pos.filt, file="p.values.assi.pos.filt.RDate")

p.values.assi.neg = getPvalues(outlist.neg.stats.more, n_patients, n_controls, assi.lab="HMDB_code",patients,controls)  

# Filter adducts 
names = lapply(rownames(p.values.assi.neg), function(x, adducts=c(1), adducts_long=c("[M+Cl]-")) {
  # x=rownames(p.values.assi.pos)[1]
  compounds=unlist(strsplit(as.vector(x),";"))
  # compounds=compounds[-1]   # since DIMS 2.1
  index=c(1:length(compounds))
  adducts_index = grep("_", compounds, fixed=TRUE)
  
  if (length(adducts_index)>0) index=index[-adducts_index]
  
  all_adducts=compounds[adducts_index]
  keep=NULL
  for (i in 1:length(all_adducts)){
    # unlist(strsplit(as.vector(all_adducts[i]),"_"))[2] %in% adducts
    keep = c(keep, unlist(strsplit(as.vector(all_adducts[i]),"_"))[2] %in% adducts)
  }
  all_adducts=all_adducts[keep]
  
  for (i in 1:length(adducts)){
    all_adducts=gsub(paste("_",toString(adducts[i]),sep=""), paste(" ", adducts_long[i], sep=""), all_adducts)
  }
  
  compounds=compounds[index]
  if (length(all_adducts)!=0){
    if (length(compounds)!=0){
      compounds=paste(compounds, all_adducts, sep=";", collapse=";")
    } else {
      compounds=paste(all_adducts, collapse=";")
    }  
  } else if (length(compounds)>1) {
    compounds=paste(compounds, collapse=";")
  }
  
  if (length(compounds)==0) compounds="" 
  
  return(compounds)
})

index=which(names=="")
if (length(index)>0) names=names[-index]
p.values.assi.neg.filt = p.values.assi.neg[-index,] 
rownames(p.values.assi.neg.filt) = names

save(p.values.assi.neg.filt, file="p.values.assi.neg.filt.RDate")

# Positive and negative mode together
rownames(p.values.assi.pos.filt) = paste(rownames(p.values.assi.pos.filt), "pos", sep="_") 
rownames(p.values.assi.neg.filt) = paste(rownames(p.values.assi.neg.filt), "neg", sep="_") 
p.values.assi.all=rbind(p.values.assi.pos.filt, p.values.assi.neg.filt)
###############################################################################################
###############################################################################################
###############################################################################################

# save(p.values.assi.all, file = "./results/p.values.assi.all_CROSS-OMICS_P28.RData")
# save(p.values.assi.all, file = "./results/p.values.assi.all_RES_BSP_SP_IX.RData")
# save(p.values.assi.all, file = "./results/p.values.assi.all_RES_BSP_20171109_SP_X")
# save(p.values.assi.all, file = "./results/p.values.assi.all_CROSS-OMICS_SPVII_P95.RData")
# save(p.values.assi.all, file = "./results/p.values.assi.all_CROSS-OMICS_SPVIII_P99(18).RData")

# load("./results/p.values.assi.all_PAH_CBS_DIMS_2.1_more_noise.RData")
# load("./results/p.values.assi.all_GCDH_IVD_DIMS_2.1_more_noise.RData")
# load("E:/Metabolomics/Crossomics_SinglePatients/results/p.values.assi.all_MAT1A_MTHFR_DIMS_2.1_more_noise.RData")
# load("E:/Metabolomics/Crossomics_SinglePatients/results/p.values.assi.all_MCT1_DIMS_2.1_more_noise.RData")
# load("./results/p.values.assi.all_SPIV_2.1_more_noise.RData")
# load("E:/Metabolomics/Crossomics_SinglePatients/results/p.values.assi.all_CROSS-OMICS_P80,82.RData")
load("Z:/Metabolomics/DIMS_pipeline/R_workspace_ME/Crossomics/Crossomics_SinglePatients/results/p.values.assi.all_CROSS-OMICS_P80,82.RData")
# load("E:/Metabolomics/Crossomics_SinglePatients/results/p.values.assi.all_CROSS-OMICS_P38.RData")
# load("./results/p.values.assi.all_RES_BSP_SP_IX.RData")
# load("./results/p.values.assi.all_CROSS-OMICS_P28.RData")
# load("./results/p.values.assi.all_RES_BSP_20171109_SP_X")
# load("./results/p.values.assi.all_CROSS-OMICS_SPVII_P95.RData")

# ###############################################################################################
# ################ Summed adducts ###############################################################
# ###############################################################################################
# adductsSummed=TRUE
# outlist.adducts.stats=cbind("HMDB_code"=rownames(outlist.adducts.stats), outlist.adducts.stats) 
# p.values.assi.all = getPvalues(outlist.adducts.stats, n_patients, n_controls, assi.lab="HMDB_code",patients,controls, adducts=TRUE) # assi.hmdb, HMDB_code For new identification !!!!!!!!!!!!!!!!!!!!!!!!!!
# ###############################################################################################
# ###############################################################################################
# ###############################################################################################

dir.create("./results", showWarnings = FALSE)
path="./results/crossomics_Fisher_weighted"
dir.create(path, showWarnings = FALSE)
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

# path2 = "../Crossomics_Build_Metabolite_Set/results/metabolite_sets_step_0,1,2,3_1.0_filter_1.1/mss_0"
path2 = "../Crossomics_Build_Metabolite_Set/results/metabolite_sets_step_0,1,2,3_1.0_filter_1.1/mss_1"
# path2 = "../Crossomics_Build_Metabolite_Set/results/metabolite_sets_step_0,1,2,3_1.0_filter_1.1/mss_2"
# path2 = "../Crossomics_Build_Metabolite_Set/results/metabolite_sets_step_0,1,2,3_1.0_filter_1.1/mss_3"

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

for (i in 1:length(patients)){

  ######################## MSEA Revon2 neighbourhood ######################################################
  dir.create(paste(path,"/P",patients[i], sep=""), showWarnings = FALSE)
  dir.create(paste(path,"/P",patients[i],"/Recon2", sep=""), showWarnings = FALSE)

  # message(paste("Patient", patients[i]))

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

    save(mss, file = paste0("./results/mss_",genes[i][[1]][1],i,".RData"))
  }

  # load(paste0("./results/mss_",genes[i][[1]][1],i,".RData"))

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

    mapper = loadDatabase("./db/BridgeDB/metabolites_20150717.bridge")  # ME loadDatabase function zit in BridgeDbR library 
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

    # 3 = Fisher weighted
    retVal = performMSEA(metaboliteSet, p.values, patients[i], gene_in, n_patients, thresh_F_pos, thresh_F_neg, path, 3, top, id, adductsSummed)

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

  # # rank
  # index.rank = min(which(tmp2[1,]%in%unlist(genes[i])))
  # rank = rank + index.rank
  # p_value_sum = p_value_sum + as.numeric(tmp2[2,index.rank])
  ##########################################################################################################
}
if (!is.null(genes)) genExcelFileShort(overview, paste(path,"/MSEA_overview.xls",sep=""))
# rank = rank/length(patients)
# save(rank, p_value_sum, file=paste(path,"/rank.RData",sep=""))
Sys.time()

########################################################################
########################################################################
########################################################################
path="./results/crossomics_Fisher_Weighted"
dir.create(path, showWarnings = FALSE)
overview = NULL

for (i in 1:length(patients)){
  
  ######################## MSEA Revon2 neighbourhood ######################################################
  dir.create(paste(path,"/P",patients[i], sep=""), showWarnings = FALSE)
  dir.create(paste(path,"/P",patients[i],"/Recon2", sep=""), showWarnings = FALSE)
  
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
    hmdb = getSystemCode("HMDB")   # ME system code for HGMD??
    kegg = getSystemCode("KEGG Compound") # ME system code for KEGG??
    
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
  genExcelFileShort(tmp[order(as.numeric(tmp[,"p.value"])),], paste(path,"/P",patients[i],"/Recon2/MSEA_results.xls",sep=""))
  
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
if (!is.null(genes)) genExcelFileShort(overview, paste(path,"/MSEA_overview.xls",sep=""))
# rank = rank/length(patients)
# save(rank, p_value_sum, file=paste(path,"/rank.RData",sep=""))
Sys.time()

########################################################################
########################################################################
########################################################################
path="./results/crossomics_Hyper"
dir.create(path, showWarnings = FALSE)
overview = NULL

for (i in 1:length(patients)){
  
  ######################## MSEA Revon2 neighbourhood ######################################################
  dir.create(paste(path,"/P",patients[i], sep=""), showWarnings = FALSE)
  dir.create(paste(path,"/P",patients[i],"/Recon2", sep=""), showWarnings = FALSE)
  
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
  genExcelFileShort(tmp[order(as.numeric(tmp[,"p.value"])),], paste(path,"/P",patients[i],"/Recon2/MSEA_results.xls",sep=""))
  
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
if (!is.null(genes)) genExcelFileShort(overview, paste(path,"/MSEA_overview.xls",sep=""))
# rank = rank/length(patients)
# save(rank, p_value_sum, file=paste(path,"/rank.RData",sep=""))
Sys.time()

########################################################################
########################################################################
########################################################################
path="./results/crossomics_Hyper_Weighted"
dir.create(path, showWarnings = FALSE)
overview = NULL

for (i in 1:length(patients)){
  
  ######################## MSEA Revon2 neighbourhood ######################################################
  dir.create(paste(path,"/P",patients[i], sep=""), showWarnings = FALSE)
  dir.create(paste(path,"/P",patients[i],"/Recon2", sep=""), showWarnings = FALSE)
  
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
  genExcelFileShort(tmp[order(tmp[,"p.value"]),], paste(path,"/P",patients[i],"/Recon2/MSEA_results.xls",sep=""))
  
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
if (!is.null(genes)) genExcelFileShort(overview, paste(path,"/MSEA_overview.xls",sep=""))
# rank = rank/length(patients)
# save(rank, p_value_sum, file=paste(path,"/rank.RData",sep=""))
Sys.time()

