# script to find out which genes have duplicate metabolites attached to them which do have different hmdb
# codes. These genes were incorrectly used in the GeneMetabMapper script and need to be identified so
# any seeds that contained any of these genes need to be rerun


# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Session Info and loading libraries --------------------------------------
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



library(rstudioapi)
library(data.table)
library(stringr)

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# functions ---------------------------------------------------------------
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

removeMets <- function(metaboliteSet, mets2remove, identifier){
  metaboliteSet <- metaboliteSet[!metaboliteSet[,identifier] %in% mets2remove[,identifier],,drop = FALSE]
  return(metaboliteSet)
}

real_duplicated <- function(set){
  duplicated(set, incomparables = c("character(0)", NA)) | duplicated(set, fromLast = TRUE, incomparables = c("character(0)", NA))
}

make.data.table <- function(datatable, identifier){
  datatable[  , .(rxn_id = toupper(paste.unique(rxn_id)),
                  hmdb = toupper(paste.unique(hmdb)),
                  kegg = toupper(paste.unique(kegg)),
                  chebi = paste.unique(chebi),
                  pubchem = paste.unique(pubchem),
                  met_long = paste.unique(met_long)
  ), 
  by = identifier]
}

paste.unique <- function(colname){
  index <- colname != "character(0)" & colname != "NA" & !is.na(colname)
  if(sum(index) > 1){
    paste(unique(tolower(colname[index])), collapse = ", ")
  } else {
    tolower(colname[index])
  }
}


# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# General -----------------------------------------------------------------
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

code_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)

# use all genes with most lax parameters (maxrxn = 19, step = 5)
gene_dir <-  paste0(code_dir,"/../Data/2019-08-12/maxrxn19/mss_5_HMDBtranslated")
gene_list <- list.files(path = gene_dir)
mets2remove <- as.data.frame(readRDS(paste0(code_dir,"/../Data/mets2remove.RDS")))

train_data_name <- "Crossomics_DBS_Marten_TraVal_Inclusion_only_updated20191031.RData"


# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Make list of all HMDB codes present in any run file ---------------------
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

load(paste0(code_dir,"/../Data/", train_data_name))
data_locations <- unique(xls_data$Location)
data_locations <- lapply(data_locations, function(x) paste(code_dir,"../Data", paste(strsplit(x, split = "\\\\")[[1]][c(6:8)], collapse = "/"), sep = "/"))

all_HMDB_codes <- NULL
for(dl in data_locations){
  RDS_file <- list.files(path = dl, pattern = "*.RDS")
  Zint_values <- readRDS(file = paste(dl, RDS_file, sep = "/"))
  tmp_HMDB_codes <- Zint_values[,"HMDB_code"]
  all_HMDB_codes <- unique(c(all_HMDB_codes, tmp_HMDB_codes))
}
all_HMDB_codes <- str_replace(all_HMDB_codes, pattern = "HMDB", replacement = "HMDB00")


# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# looping over genes ------------------------------------------------------
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


gen_dup_list <- NULL
counter <- 0

for(ge_fi in gene_list){
  counter <- counter + 1
  if(counter %% 100 == 0) cat("gene", counter, "out of", length(gene_list), "\n")
  ge_me_set <- readRDS(paste(gene_dir,ge_fi, sep = "/"))
  ge_me_set <- ge_me_set[,c("rxn_id","hmdb","kegg","chebi","pubchem","met_long")]
  
  if(is.vector(ge_me_set)) {
    dimnames <- names(ge_me_set)
    ge_me_set <- matrix(ge_me_set, nrow = 1)
    colnames(ge_me_set) <- dimnames
    rm(dimnames)
  }
  
  if(all(ge_me_set[,c("hmdb","kegg","chebi")] == "character(0)")) next
  
  # Remove any metabolites that are non-informative
  for(identifier in c("chebi","kegg","pubchem")){
    ge_me_set <- removeMets(ge_me_set, mets2remove, identifier)
  }
  
  # Get rid of any rows that don't have an hmdb code and go to next gene if none are left
  index <- which(ge_me_set[,"hmdb"] == "character(0)")
  if (length(index)>0) ge_me_set <- ge_me_set[-index,,drop=FALSE]
  if (nrow(ge_me_set) == 0) next
  
  # # Or that do not appear in the dataset
  index <- which(!as.vector(unlist(lapply(ge_me_set[,"hmdb"], function(x) length(grep(x, all_HMDB_codes)) > 0))))
  if (length(index)>0) ge_me_set <- ge_me_set[-index,,drop=FALSE]
  if (nrow(ge_me_set) == 0) next
  
  # Manual fix for the metabolites known under "HMDB0012482, HMDB0002281" (same, but different in the hmdb dataset)
  if(length(grep("HMDB0002281", ge_me_set[,"hmdb"])) | length(grep("HMDB0012482", ge_me_set[,"hmdb"]))){
    HMDB_index <- unique(grep("HMDB0002281", ge_me_set[,"hmdb"]),grep("HMDB0012482", ge_me_set[,"hmdb"]))
    ge_me_set[HMDB_index,"hmdb"] <- "HMDB0002281"
  }
  
  # Remove duplicate metabolites, but collate their data if they are known under different ID's/names
  if(nrow(ge_me_set) >1){
    tmp <- as.data.frame(ge_me_set[apply(apply(ge_me_set[,c("hmdb","chebi","kegg")], 2, real_duplicated), 1, any),])
    setDT(tmp)
    tmp <- make.data.table(tmp, "hmdb")
    tmp <- make.data.table(tmp, "chebi")
    ge_me_set <- ge_me_set[apply(!apply(ge_me_set[,c("hmdb","chebi","kegg")], 2, real_duplicated), 1, all),]
    
    ge_me_set <- rbind(ge_me_set, as.matrix(tmp[,-1]))
    rm(tmp)
  }
  met_vec <- as.vector(unlist(lapply((strsplit(ge_me_set[,"hmdb"], split = ",")), length))) > 1
  # if(any(as.vector(unlist(lapply((strsplit(ge_me_set[,"hmdb"], split = ",")), length))) > 1)){
  if(any(met_vec)){
    gen_dup_list <- rbind(gen_dup_list, ge_me_set[met_vec,, drop = FALSE])
    # gen_dup_list <- c(gen_dup_list, ge_fi)
  }
}


gen_dup_DT <- as.data.table(gen_dup_list)

gen_dup_DT <- gen_dup_DT[, list(Gene = paste(rxn_id, collapse = "; "), Name = paste(unique(met_long), collapse = ";")), by = .(hmdb)]


