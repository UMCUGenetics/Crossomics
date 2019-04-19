# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Session info ------------------------------------------------------------
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

cat(
  "
  Created by:   Marcel Willemse?
  Modified by:  Marten Kerkhofs, 2019-04-19

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
patients=c(20)




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
patient <- patients[1]
adductsHMDB_z_ave_and_int <- generate_av_Z_scores(patient = patient)



# file <- "/Users/mkerkho7/DIMS2_repo/TestResults/20180720_Run4_Diagnosis2017_1_Adductsums_Zscores.rds"
# adductsHMDB_z_ave_and_int <- readRDS(file)




# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Perform MSEA ------------------------------------------------------------
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

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


# Then performe MSEA ######################################################################################

path2 <- "../Crossomics_Build_Metabolite_Set/results/metabolite_sets_step_0,1,2,3_1.0_filter_1.1/mss_1"
# path2 = paste0("../Crossomics_Build_Metabolite_Set/results/metabolite_sets_step_0,1,2,3_1.0_filter_1.1/mss_",distance_to_gene)
step <- 1 # 1:3  # ME wordt gebruikt in regel 495 -503, waarom?
path_wg <- "./results/mss_WG_step_0" # for sampling random genes   #ME enkel voor random data set, niet WES
sets <- list.files(path = path2)
overview <- NULL # ME ?


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