# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Session info ------------------------------------------------------------
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# cat(
#   "
#   Created by:   Marcel Willemse?
#   Modified by:  Marten Kerkhofs, 2019-04-26
# 
#   Copied from 'GeneMetaboliteMapper_ME_Marten.R, on 29/03/2019 which in turn was copied from:
#   'GeneMetaboliteMapper_ME.R' in the Metab/Metabolomics/DIMS_pipeline/R_workspace_ME/Crossomics/Crossomics_SinglePatients/src folder
#   on 26/03/2019 to investigate and annotate the functionality of the crossomics pipeline.
#   Disabled all save/generateExcel functions
#   
#   OS
#   macOS 10.14.4
#   
#   Package versions:
#   R version   3.5.1 (2018-07-02)
#   
#   stringr       1.4.0
#   Cairo         1.5-9
#   XLConnect     0.2-15
#   BridgeDbR     1.16.1
#   rJava         0.9-10
#   XLConnectJars 0.2-15
#   bioDist       1.54.0
#   KernSmooth    2.23-15
#   Biobase       2.42.0
#   BiocGenerics  0.28.0
#   
#   
#   USE:
#   This file takes a patient's metabolite and rare-gene variant data and performs metabolite set enrichment analysis
#   OR
#   This file takes a patient without WES data and creates a mock-gene set + diseased gene to perform metabolite set enrichment analysis
#   
#   Input:
#   Varying functions in the src folder
#   metabolite data from a patient: adductSums_xx.RData (xx = negative or positive) files
#   A gene set of a single patient which contains rare-gene variants in that patient
#   A patient number
#   ...
#   
#   Output
# 
#   "
# )




# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Libraries ---------------------------------------------------------------
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# library("R.matlab")
# library("KEGGREST")
library("bioDist")
library("Cairo")
library("XLConnect")
library("BridgeDbR")
library("stringr") # add leading 0's

# source("./src/Crossomics/sourceDir.R")
# sourceDir("./src/Crossomics")
source("/Users/mkerkho7/DIMS2_repo/Crossomics/Scripts/Supportive/sourceDir.R")
sourceDir("/Users/mkerkho7/DIMS2_repo/Crossomics/Scripts/Supportive")



# for(distance in 0:3){
  

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Manual inputs -----------------------------------------------------------
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# To what distance from the reaction should metabolites be considered? (0 to 3)
distance_to_gene <- distance

# Save mock gene set?
save_mock <- TRUE

# Get the correct data location on the basis of the patient number alone. 
patients <- 2

# Get correct dataset if patient is present in multiple (choose 0 or 1)(ignore when there is only 1 run present)
run <- 0

# Make the sampling of random genes repeatable
seed <- 313

# Number of mock genes to add
nr_mocks <- 100




# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Set working directory ---------------------------------------------------
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

try(setwd("/Volumes/Metab/Metabolomics/DIMS_pipeline/R_workspace_ME/Crossomics/Crossomics_SinglePatients"), silent = TRUE)

# Test for getting the data from P01.1, having PKU
# patient_folder <- "/Volumes/Metab/Metabolomics/Research Metabolic Diagnostics/Metabolomics Projects/Projects 2017/Project 2017_008 MetabolomicsDiagnosis_DBS/RES_BSP_20170629_MetabolomicsDiagnosis_RUN0/Bioinformatics_metIS/"
# setwd(patient_folder)

patient_file <- "/Users/mkerkho7/DIMS2_repo/Crossomics/Data/Crossomics_DBS_Marten_Training.xlsx"
wb = loadWorkbook(patient_file, create = TRUE)
data = readWorksheet(wb, sheet = 1, startRow = 0, endRow = 0, startCol = 0, endCol = 0)


# Search for patient number with and without leading 0 (when patientnumber is single digid)
if (nchar(patients) == 1){
  dummy_patient_number <- c(paste0("P", patients,"\\."), paste0("P0", patients,"\\."))
} else {
  dummy_patient_number <- paste0("P",patients)
}

patient_rows <- grep(paste(dummy_patient_number, collapse = "|"), data$Patient.number)
patient_datasets <- unique(data$Dataset[patient_rows])

# Get the correct dataset name and location
if(length(patient_datasets) > 1){
  patient_datasets <- patient_datasets[grep(paste0("RUN", run), patient_datasets)]
}

data_location <- data$Location[match(patient_datasets, data$Dataset)]

# Get disease gene for patient
dis_gene <- unique(data$Gene[grepl(paste(dummy_patient_number, collapse = "|"), data$Patient.number) & data$Location == data_location])

# Make dataset location mac-compatible
data_location <- gsub("Y:", "/Volumes/Metab", data_location)
data_location <- gsub("\\", "/", data_location, fixed = TRUE)

setwd(data_location)




# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Calculate Z scores ------------------------------------------------------
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# patient <- patients[1]
# Make sure to get a patient with at least 2 digits (incl leading 0 if necessary)
patient <- paste0("P", str_pad(patients, width = 2, pad = "0"))
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

# For getting random mock genes
# mock_genes <- read.table(file = "/Users/mkerkho7/DIMS2_repo/Crossomics/Data/All_Genes_Ensembl_apr_2019_GRCh38p12.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
mock_genes <- read.table(file = "/Users/mkerkho7/DIMS2_repo/Crossomics/Data/All_Genes_Ensembl_apr_2019_GRCh38p12_extended.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
mock_genes <- mock_genes[mock_genes$Gene.type == "protein_coding",]
mock_genes <- mock_genes[!mock_genes$HGNC.ID == "",]
mock_genes <- mock_genes[!duplicated(mock_genes$Gene.name),]
mock_genes <- mock_genes[!grepl("orf", mock_genes$Gene.name),]



mock_genes <- mock_genes$Gene.name

# make sure that the disease gene isn't included twice by removing it from the mock_genes list
mock_genes <- mock_genes[mock_genes != dis_gene]
set.seed(seed = seed)
genes <- sample(mock_genes, size = nr_mocks)
genes <- c(genes, dis_gene)





# Then perform MSEA ######################################################################################

# path2 <- "../Crossomics_Build_Metabolite_Set/results/metabolite_sets_step_0,1,2,3_1.0_filter_1.1/mss_1"
# path2 = paste0("../Crossomics_Build_Metabolite_Set/results/metabolite_sets_step_0,1,2,3_1.0_filter_1.1/mss_",distance_to_gene)
path2 <- paste0("/Volumes/Metab/Metabolomics/DIMS_pipeline/R_workspace_ME/Crossomics/Crossomics_Build_Metabolite_Set/results/metabolite_sets_step_0,1,2,3_1.0_filter_1.1/mss_", distance_to_gene)
# step <- 1 # 1:3  # ME wordt gebruikt in regel 495 -503, waarom?
step <- distance_to_gene
# path_wg <- "./results/mss_WG_step_0" # for sampling random genes   #ME enkel voor random data set, niet WES
# sets <- list.files(path = path2)
overview <- NULL # at the end


###########################################################################################################
########################### real WES data #################################################################
###########################################################################################################
Sys.time()
p.values.assi.all <- adductsHMDB_z_ave_and_int
for (i in 1:length(patients)){
  
  patient_folder <- paste(patient[i], patient_datasets, "dis", distance_to_gene, sep = "_")
  
  ######################## MSEA Revon2 neighbourhood ######################################################
  # dir.create(paste(path,"/P",patients[i], sep=""), showWarnings = FALSE)
  # dir.create(paste(path,"/P",patients[i],"/Recon2", sep=""), showWarnings = FALSE)
  # dir.create(paste(path,"/",patient[i], sep=""), showWarnings = FALSE)
  # dir.create(paste(path,"/",patient[i],"/Recon2", sep=""), showWarnings = FALSE)
  dir.create(paste(path,"/",patient_folder, sep=""), showWarnings = FALSE)
  # dir.create(paste(path,"/",patient_folder,"/Recon2", sep=""), showWarnings = FALSE)
  
  # save mock genes + real gene
  # if(save_mock) write.table(genes, file = paste0(path, "/", patient[i], "/Recon2/genes_used.txt"), row.names = FALSE, col.names = FALSE)
  # if(save_mock) write.table(genes, file = paste0(path, "/", patient_folder, "/Recon2/genes_used.txt"), row.names = FALSE, col.names = FALSE)
  if(save_mock) write.table(genes, file = paste0(path, "/", patient_folder, "/genes_used.txt"), row.names = FALSE, col.names = FALSE)
  
  # subset to 1 patient!
  # p.values <- p.values.assi.all[,c(grep(paste("P",patients[i],sep=""), colnames(p.values.assi.all), fixed=TRUE),
  #                               grep("C", colnames(p.values.assi.all), fixed=TRUE),
  #                               which(colnames(p.values.assi.all)=="avg.int.controls"))]
  p.values <- p.values.assi.all[,c(grep(patient[i], colnames(p.values.assi.all), fixed=TRUE),
                                   grep("C", colnames(p.values.assi.all), fixed=TRUE),
                                   which(colnames(p.values.assi.all)=="avg.int.controls"))]
  
  if (is.null(genes)){
    mss = read.table(paste("./db/P",patients[i],"_HGNC.txt",sep=""), header = FALSE, sep="\t")
    mss = as.vector(unlist(mss))
    mss = paste(mss,"RData",sep=".")
    # ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
    # gene_map_table = getBM(attributes=c('hgnc_symbol', 'entrezgene', 'ensembl_gene_id'),
    #                        filters = 'hgnc_symbol', values = as.vector(unlist(gene_list)), mart = ensembl)
    
    # when mock genes are supplied / WES data unavailable
  } else {
    genespp <- genes
    # genespp=as.vector(unlist(genes[i])) # old
    
    mss = NULL
    for (k in 1:length(genespp)){
      mss = c(mss, paste(genespp[k],"RData",sep="."))
      # print(mss)
    }
    
    # wg = list.files(path = path_wg)
    # mss = c(mss, wg[sample(1:length(wg), 100-length(genespp), replace=TRUE)]) # Why with replacement? 
    
    # selection = which(sets %in% mss)
    # mss = sets[selection]
    
    # save(mss, file = paste0("./results/mss_",genes[i][[1]][1],i,".RData"))
  }
  
  metSetResult = NULL
  nMets = NULL    # list with number of metabolites per gene, not used for any calculations, but only for output excel file.
  for (j in 1:length(mss)){
    # print(mss[j])
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
        metaboliteSet <- data.frame(t(metaboliteSet), stringsAsFactors = FALSE)  # ME matrix transpose
      }
    }
    
    ##################################################
    # temporarily work around to be fixed in findMetabolicEnvironment <- why?
    # NA's in metaboliteSet[,"hmdb"] are set to 'character(0)'
    # index = which(is.na(metaboliteSet[,"hmdb"]))  
    # if (length(index) > 0) metaboliteSet[index,"hmdb"] = "character(0)"
    
    # Set all NA's to "character (0)" in the ID columns
    metaboliteSet[,c("hmdb","kegg","chebi")][is.na(metaboliteSet[,c("hmdb","kegg","chebi")])] <- "character(0)"
    
    # Check for when genes have no metaboliteSet info, continue with next gene
    if (all(metaboliteSet[,c("hmdb","kegg","chebi")] == "character(0)")) next

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
    
    # replace 'new' hmdb code-format with old ones (remove two zero's)
    metaboliteSet[nchar(metaboliteSet[,"hmdb"]) == 11,"hmdb"] <- str_replace(metaboliteSet[nchar(metaboliteSet[,"hmdb"]) == 11,"hmdb"], pattern = "B00", replacement = "B")
    
    # Get rid of any rows that still don't have an hmdb code
    index = which(metaboliteSet[,"hmdb"] == "character(0)")
    
    if (length(index)>0) metaboliteSet <- metaboliteSet[-index,,drop=FALSE]
    
    # Get the total number of unique metabolites for all genes in one vector
    # nMets=c(nMets, length(unique(metaboliteSet[,"hmdb"])))
    # Improved check for same compounds, including if either chebi, kegg or hmdb have duplicated values. Return number of metabolites
    # nMets <- c(nMets, sum(apply(!apply(metaboliteSet[,c("hmdb","chebi","kegg")], 2, duplicated, incomparables = NA), 1, all)))
    
    # Remove duplicate metabolites (added by Marten)
    if(nrow(metaboliteSet) >1){
      metaboliteSet <- metaboliteSet[apply(!apply(metaboliteSet[,c("hmdb","chebi","kegg")], 2, duplicated, incomparables = c("character(0)", NA)), 1, all),, drop = FALSE]
    }
    nMets <- c(nMets,nrow(metaboliteSet))
    
    savepoint_p.values <- p.values
    savepoint_metaboliteSet <- metaboliteSet
    
    
    
    # Peform MSEA, one single patient and rare-gene variant at a time
    # 0 = Fishers exact; 3 = Fisher weighted
    # Top stands for the top genes to perform the statistical test with, this is to normalise the different metabolite set-sizes
    retVal = performMSEA(metaboliteSet = savepoint_metaboliteSet, 
                         av_int_and_z_values_matrix = savepoint_p.values, 
                         patient = patient[i], 
                         gene_in, 
                         n_patients, 
                         thresh_F_pos, 
                         thresh_F_neg, 
                         path, 
                         test = 0, 
                         top, 
                         id, 
                         patient_folder = patient_folder
    )
    print(mss[j])

    # Get only the p.value and put all p.values with the patient and gene in one table
    print(retVal)
    p_value = as.numeric(retVal$p.value)
    if (length(p_value)==0){
      p_value=NA
    }
    metSetResult <- rbind(metSetResult, c("p.value"=p_value, "patient"=paste0("P", patients[i]), "metabolite.set"=gene_in))
    
  }
  
  # Create a short excel file with the p.values and the number of metabolites associated with a gene
  tmp <- data.frame("HGNC"=metSetResult[,3],"p.value"=as.numeric(metSetResult[,"p.value"]), "metabolites"=nMets)
  # genExcelFileShort(tmp[order(tmp[,"p.value"]),], paste(path,"/",patient[i],"/Recon2/MSEA_results.xls",sep=""))
  # genExcelFileShort(tmp[order(tmp[,"p.value"]),], paste(path,"/P",patients[i],"/Recon2/MSEA_results.xls",sep=""))
  # genExcelFileShort(tmp[order(tmp[,"p.value"]),], paste(path,"/",patient_folder,"/Recon2/MSEA_results.xls",sep=""))
  genExcelFileShort(tmp[order(tmp[,"p.value"]),], paste(path,"/",patient_folder,"/MSEA_results.xls",sep=""))
  
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


# }