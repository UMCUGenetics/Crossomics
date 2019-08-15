# Script to fix defect HMDB codes 

library("BridgeDbR")
library("stringr")


date <- "2019-08-12" # The date of the data/mss_0 etc. runs
steps <- c(0,1,2,3,4,5)
max_rxns <- c(8,10,12,15,17,19)


setwd("/Users/mkerkho7/DIMS2_repo/Crossomics/Data/")
# mss <- list.files(paste0(date, "_maxrxn8/mss_0_HMDBtranslated/"))
# if(unlist(strsplit(mss[1], split = "\\."))[2] == "RData"){
#   save_as <- "RData"
# } else {
#   save_as <- "RDS"
# }
mss <- list.files("mss")
if(date >= "2019-08-12") {
  mss <- paste0(unlist(lapply(strsplit(mss, split = "\\."), function(x) unlist(x)[1])),".RDS")
  save_as <- "RDS"
} else {
  save_as <- "RData"
}

# mss <- list.files("2019-07-19_maxrxn15/mss_4_HMDBtranslatedtest/")
genes <- unlist(lapply(mss, function(x) strsplit(x,"\\.")[[1]][1]))


if(Sys.getenv("RSTUDIO") == "1") {
  code_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
} else {
  code_dir <- "/Users/mkerkho7/DIMS2_repo/Crossomics/src_Metabolite_Mapper"
}

# mets2remove <- as.data.frame(readRDS("./mets2remove.RDS"))
mapper <- BridgeDbR::loadDatabase("/Users/mkerkho7/DIMS2_repo/Crossomics/Data/metabolites_20190509.bridge")
hmdb_id <- BridgeDbR::getSystemCode("HMDB")
kegg_id <- BridgeDbR::getSystemCode("KEGG Compound")
chebi_id <- BridgeDbR::getSystemCode("ChEBI")

# Function to translate old to new HMDB
translateHMDB <- function(hmdb_code){
  if(nchar(hmdb_code) != 11){
    if(nchar(hmdb_code) == 9){
      hmdb_code <- str_replace(hmdb_code, pattern = "HMDB", replacement = "HMDB00")
      metaboliteSet[index,"hmdb"] <- hmdb_code
    } else {
      hmdb_code <- NULL
      warning("Metabolite does not have correct HMDB code, deleting it")
    }
  }
  return(hmdb_code)
}

# loop over all parameters to get all gene files
for(j in mss){
  for(maxrxn in max_rxns){
    for(step in steps){
      if(date >= "2019-08-12") {
        filename <- paste0(date,"/maxrxn",maxrxn,"/mss_",step,"_HMDBtranslated/",j)
      } else {
        filename <- paste0(date,"_maxrxn",maxrxn,"/mss_",step,"_HMDBtranslated/",j)
      }
      
      if(save_as == "RData"){
        load(filename)
      } else {
        metaboliteSet <- as.matrix(readRDS(filename))
      }
      indices <- which(!xor(nchar(metaboliteSet[,"hmdb"]) != 11, metaboliteSet[,"hmdb"] != "character(0)"))
      # move on to next file unless a faulty HMDB code is found
      if(length(indices)==0){
        next
      } else {
        message("faulty HMDB code found in file: ", j, " metabolite(s): ", paste(metaboliteSet[indices,"hmdb"], collapse = ", "))
        metaboliteSet[indices, "hmdb"] <- "character(0)"
        for(index in indices){
          fixed_hmdb <- NULL
          id <- metaboliteSet[index,"chebi"]
          if(id != "character(0)"){
            fixed_hmdb <- unlist(BridgeDbR::map(mapper, chebi_id, id, hmdb_id)[1])
          } 
          if(!is.null(fixed_hmdb)){
            fixed_hmdb <- translateHMDB(fixed_hmdb)
          } else {
            id <- metaboliteSet[index,"kegg"]
            if(id != "character(0)"){
              fixed_hmdb <- unlist(BridgeDbR::map(mapper, kegg_id, id, hmdb_id)[1])
            }
            if(!is.null(fixed_hmdb)){
              fixed_hmdb <- translateHMDB(fixed_hmdb)
            } 
          }
          if(!is.null(fixed_hmdb)){
            metaboliteSet[index, "hmdb"] <- fixed_hmdb
          }
        }
        if(save_as == "RData"){
          save(metaboliteSet, file = filename)
        } else {
          saveRDS(metaboliteSet, file = filename)
        }
      }
    }
  }
}
