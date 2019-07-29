# Script to convert all kegg and chebi codes in the metabolite sets to (new) HMDB codes where possible.
# The need for this separate script is to make the metabolite mapper script independent from any library that uses Java.

library("BridgeDbR")
library("stringr")


date <- "2019-07-19" # The date of the data/mss_0 etc. runs

setwd("/Users/mkerkho7/DIMS2_repo/Crossomics/Data/")
mss <- list.files("mss")
genes <- unlist(lapply(mss, function(x) strsplit(x,"\\.")[[1]][1]))


if(Sys.getenv("RSTUDIO") == "1") {
  code_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
} else {
  code_dir <- "/Users/mkerkho7/DIMS2_repo/Crossomics/src_Metabolite_Mapper"
}

# mets2remove <- as.data.frame(readRDS("./mets2remove.RDS"))
mapper <- BridgeDbR::loadDatabase("/Users/mkerkho7/DIMS2_repo/Crossomics/Data/metabolites_20190509.bridge")
hmdb_id <- BridgeDbR::getSystemCode("HMDB")




# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Functions ---------------------------------------------------------------
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

getHMDBcode <- function(metaboliteSet, identifier){
  index <- which(metaboliteSet[,"hmdb"] == "character(0)") 
  index.sub <- which(metaboliteSet[index, identifier] != "character(0)") 
  id <- metaboliteSet[index[index.sub], identifier]
  
  if(identifier == "kegg"){
    kegg_chebi <- BridgeDbR::getSystemCode("KEGG Compound")
  } else if (identifier == "chebi"){
    kegg_chebi <- BridgeDbR::getSystemCode("ChEBI")
  }
  
  # Try to fill in empty HMDB IDs via the KEGG ID
  if (length(id) > 0){
    for (k in 1:length(id)){
      if (!is.null(unlist(BridgeDbR::map(mapper, kegg_chebi, id[k], hmdb_id)[1]))) {
        metaboliteSet[index[index.sub[k]],"hmdb"] <- unlist(BridgeDbR::map(mapper, kegg_chebi, id[k], hmdb_id)[1])
      }
    }
  }
  return(metaboliteSet)
}



# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Code --------------------------------------------------------------------
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

for (j in 1:length(mss)){
  for(max_rxn_number in c(8,10,12,15)){
    max_rxn <- paste0("maxrxn",max_rxn_number)
    for(step in c(0,1,2,3,4)){
      path2 <- paste0("./",date,"_",max_rxn,"/mss_",step)

      
      
      load(paste(path2, mss[j], sep="/"))
      
      gene_in <- genes[j]
      
      # result_mets_x is the name of how the RData file is loaded in and depends on the step size.
      if (step==1){
        metaboliteSet <- result_mets_1
      } else if (step==2){
        metaboliteSet <- result_mets_2
      } else if (step==3){
        metaboliteSet <- result_mets_3
      } else if (step==4){
        metaboliteSet <- result_mets_4
      } else {
        metaboliteSet <- result_mets_0
      }
      
      if(is.vector(metaboliteSet)) { # for when a single metabolite is present (it could be converted to a character vector)
        dimnames <- names(metaboliteSet)
        metaboliteSet <- matrix(metaboliteSet, nrow = 1)
        colnames(metaboliteSet) <- dimnames
        rm(dimnames)
      }
      
      # Strip the "chebi" column of the text "CHEBI:", which is present in some of the rows
      metaboliteSet[,"chebi"] <- gsub("CHEBI:", "", metaboliteSet[,"chebi"], fixed = TRUE)
      
      # Set all NA's to "character (0)" in the ID columns
      metaboliteSet[,c("hmdb","kegg","chebi")][is.na(metaboliteSet[,c("hmdb","kegg","chebi")])] <- "character(0)"
      
      metaboliteSet <- getHMDBcode(metaboliteSet = metaboliteSet, identifier = "kegg")
      metaboliteSet <- getHMDBcode(metaboliteSet = metaboliteSet, identifier = "chebi")
      
      metaboliteSet[nchar(metaboliteSet[,"hmdb"]) == 9,"hmdb"] <- str_replace(metaboliteSet[nchar(metaboliteSet[,"hmdb"]) == 9,"hmdb"], pattern = "HMDB", replacement = "HMDB00")
      
      dir.create(paste0(path2,"_HMDBtranslated"), showWarnings = FALSE)
      save(metaboliteSet, file = paste0(path2,"_HMDBtranslated/",mss[j]))
    }
  }
}
