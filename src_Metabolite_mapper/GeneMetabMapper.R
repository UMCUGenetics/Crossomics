library("XLConnect") # excel connect function
library("stringr") # string manipulation, add leading 0's
library("BridgeDbR")
library("Cairo")

code_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
source(paste0(code_dir,"/Supportive/sourceDir.R"))
sourceDir(paste0(code_dir,"/Supportive"))

# source("/Users/mkerkho7/DIMS2_repo/Crossomics/src/src_Metabolite_Mapper/Supportive/sourceDir.R")
# sourceDir("/Users/mkerkho7/DIMS2_repo/Crossomics/src/src_Metabolite_Mapper/Supportive")

patient_file <- paste0(code_dir,"/../Data/Crossomics_DBS_Marten_Training.xlsx")
# patient_file <- "/Users/mkerkho7/DIMS2_repo/Crossomics/Data/Crossomics_DBS_Marten_Training.xlsx"
# wb = loadWorkbook(patient_file, create = TRUE)
wb = loadWorkbook(patient_file)
xls_data = readWorksheet(wb, sheet = 1, startRow = 0, endRow = 0, startCol = 0, endCol = 0)

dat_pat <- paste(xls_data$Dataset, xls_data$Patient.number, sep = "^")
uni_dat_pat <- unique(sapply(strsplit(dat_pat, split = "\\."), `[`, 1))

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Perform metabolite mapper on all dat_pat -------------------------------
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


for (i in uni_dat_pat){
  
  # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # Get patient and dataset -------------------------------------------------
  # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  tmp <- unlist(strsplit(i, split = "\\^"))
  dataset <- tmp[1]
  patient <- tmp[2]
  rm(tmp)
  cat(dataset, patient, "\n")
  
  # Get disease gene for patient
  dis_gene <- xls_data$Gene[grepl(patient, xls_data$Patient.number) & xls_data$Dataset == dataset][1]
  
  data_location <- xls_data$Location[match(dataset, xls_data$Dataset)]
  # Make dataset location mac-compatible
  data_location <- gsub("Y:", "/Volumes/Metab", data_location)
  data_location <- gsub("\\", "/", data_location, fixed = TRUE)
  
  setwd(data_location)
  
  
  # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # Calculate Z scores ------------------------------------------------------
  # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  # patient <- paste0("P", str_pad(patient, width = 2, pad = "0"))
  # For this function to work, you need to be in the folder where the patient data is present
  Zint_scores <- generate_av_Z_scores(patient = patient)
  rownames(Zint_scores)[nchar(rownames(Zint_scores)) == 9] <- str_replace(rownames(Zint_scores)[nchar(rownames(Zint_scores)) == 9], pattern = "HMDB", replacement = "HMDB00")
  
  tmp <- as.data.frame(apply(Zint_scores, 1, paste, collapse = ","))
  colnames(tmp) <- "values"
  tmp <- aggregate(rownames(tmp), by=tmp['values'], paste, collapse = ",")
  rownames(tmp) <- tmp$x
  tmp$x <- NULL
  
  tmp <- tidyr::separate(tmp, values, into = colnames(Zint_scores), sep = ",", convert = TRUE)
  Zint_pruned <- as.matrix(tmp[order(rownames(tmp)),])
  rm(tmp)
  
  # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # Perform MSEA on all distances -------------------------------------------
  # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  # path <- "/Users/mkerkho7/DIMS2_repo/TestResults/"
  path <- paste0(code_dir,"/../Results/")
  thresh_F_pos <- 1.5
  thresh_F_neg <- -1
  top <- 20
  id <- "hmdb"
  seed = 313
  nr_mocks = 100
  save_mock = TRUE
  
  
  # Prepare mock gene set ---------------------------------------------------
  genes <- NULL
  if (!file.exists(paste0("./db/P",patient,"_HGNC.txt"))){
    # For getting random mock genes
    mock_genes <- read.table(file = paste0(code_dir,"/../Data/All_Genes_Ensembl_apr_2019_GRCh38p12_extended.txt"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    # mock_genes <- read.table(file = "/Users/mkerkho7/DIMS2_repo/Crossomics/Data/All_Genes_Ensembl_apr_2019_GRCh38p12_extended.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
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
    mss <- paste(genes,"RData",sep=".")
    
    # save mock genes + real gene
    patient_folder <- paste(patient, dataset, sep = "_")
    dir.create(paste(path,"/",patient_folder, sep=""), showWarnings = FALSE)
    
    if(save_mock) write.table(genes, file = paste0(path, "/", patient_folder, "/mock_genes_seed",seed,".txt"), row.names = FALSE, col.names = FALSE)
  } else {
    # Real WES ----------------------------------------------------------------
    mss = read.table(paste0("./db/P",patient,"_HGNC.txt"), header = FALSE, sep="\t")
    mss = as.vector(unlist(mss))
    mss = paste(mss,"RData",sep=".")
  }

  
  
  for (step in c(0:4)){
    cat("step:", step, "\n")
    step_folder <- paste0(patient_folder,"/step_", step)
    
    path2 <- paste0(code_dir,"/../Data/mss_", step)
    # path2 <- paste0("/Users/mkerkho7/DIMS2_repo/Crossomics/Results/mss_", step)
    overview <- NULL # at the end
    
    dir.create(paste(path,"/",step_folder, sep=""), showWarnings = FALSE)
    
    metSetResult = NULL
    nMets = NULL    # list with number of metabolites per gene, not used for any calculations, but only for output excel file.

    for (j in 1:length(mss)){
      if(j%%5 == 0) cat(paste0(j, "%... "))
      # cat("gene:", mss[j], "number:", j,"\n")
      # Skip the gene if there is no metabolite pathway xls_data available, elsewise, load its file
      if (!file.exists(paste(path2, mss[j], sep="/"))) next
      
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
      
      # Strip the "chebi" column of the text "CHEBI", which is present in some of the rows
      metaboliteSet[,"chebi"] <- gsub("CHEBI:", "", metaboliteSet[,"chebi"], fixed = TRUE)
      
      
      # Set all NA's to "character (0)" in the ID columns
      metaboliteSet[,c("hmdb","kegg","chebi")][is.na(metaboliteSet[,c("hmdb","kegg","chebi")])] <- "character(0)"
      
      # Remove any faulty HMDB codes (they must have 9 characters at this point: the old numbering)
      metaboliteSet[nchar(metaboliteSet[,"hmdb"]) != 9 & metaboliteSet[,"hmdb"] != "character(0)", "hmdb"] <- "character(0)"
      
      # Add HMDB codes when possible, remove metabolites otherwise
      if (all(metaboliteSet[,c("hmdb","kegg","chebi")] == "character(0)")) next
      
      mapper <- BridgeDbR::loadDatabase(paste0(code_dir,"/../Data/metabolites_20190509.bridge"))
      # mapper <- BridgeDbR::loadDatabase("/Users/mkerkho7/DIMS2_repo/Crossomics/metabolites_20190509.bridge")
      hmdb <- BridgeDbR::getSystemCode("HMDB")
      
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
            if (!is.null(unlist(BridgeDbR::map(mapper, kegg_chebi, id[k], hmdb)[1]))) {
              metaboliteSet[index[index.sub[k]],"hmdb"] <- unlist(BridgeDbR::map(mapper, kegg_chebi, id[k], hmdb)[1])
            }
          }
        }
        return(metaboliteSet)
      }
      
      metaboliteSet <- getHMDBcode(metaboliteSet = metaboliteSet, identifier = "kegg")
      metaboliteSet <- getHMDBcode(metaboliteSet = metaboliteSet, identifier = "chebi")
      
      
      # # If hmdb codes are absent, try to translate KEGG and later ChEBI codes to hmdb.
      # index = which(metaboliteSet[,"hmdb"] == "character(0)") 
      # index.sub = which(metaboliteSet[index,"kegg"] != "character(0)") 
      # kegg_id = metaboliteSet[index[index.sub],"kegg"]
      # 
      # # ADDED (new .bridge file)
      # 
      # kegg = BridgeDbR::getSystemCode("KEGG Compound")
      # 
      # # Try to fill in empty HMDB IDs via the KEGG ID
      # if (length(kegg_id) > 0){
      #   for (k in 1:length(kegg_id)){
      #     if (!is.null(unlist(BridgeDbR::map(mapper, kegg, kegg_id[k], hmdb)[1]))) {
      #       metaboliteSet[index[index.sub[k]],"hmdb"] <- unlist(BridgeDbR::map(mapper, kegg, kegg_id[k], hmdb)[1])
      #     }
      #   }
      # }
      # 
      # index = which(metaboliteSet[,"hmdb"] == "character(0)")
      # index.sub = which(metaboliteSet[index,"chebi"] != "character(0)")
      # chebi_id = metaboliteSet[index[index.sub],"chebi"]
      # 
      # chebi = BridgeDbR::getSystemCode("ChEBI")
      # if (length(chebi_id)>0){
      #   for (k in 1:length(chebi_id)){
      #     if (!is.null(unlist(BridgeDbR::map(mapper, chebi, chebi_id[k], hmdb)[1]))) {
      #       metaboliteSet[index[index.sub[k]],"hmdb"] = unlist(BridgeDbR::map(mapper, chebi, chebi_id[k], hmdb)[1])
      #     }
      #   }
      # }
      # 
      # replace 'new' hmdb code-format with old ones (remove two zero's)
      # metaboliteSet[nchar(metaboliteSet[,"hmdb"]) == 11,"hmdb"] <- str_replace(metaboliteSet[nchar(metaboliteSet[,"hmdb"]) == 11,"hmdb"], pattern = "B00", replacement = "B")
      metaboliteSet[nchar(metaboliteSet[,"hmdb"]) == 9,"hmdb"] <- str_replace(metaboliteSet[nchar(metaboliteSet[,"hmdb"]) == 9,"hmdb"], pattern = "HMDB", replacement = "HMDB00")
      
      # Get rid of any rows that still don't have an hmdb code
      index = which(metaboliteSet[,"hmdb"] == "character(0)")
      
      if (length(index)>0) metaboliteSet <- metaboliteSet[-index,,drop=FALSE]
      if (nrow(metaboliteSet) == 0) next
      
      # Get the total number of unique metabolites for all genes in one vector
      # nMets=c(nMets, length(unique(metaboliteSet[,"hmdb"])))
      # Improved check for same compounds, including if either chebi, kegg or hmdb have duplicated values. Return number of metabolites
      # nMets <- c(nMets, sum(apply(!apply(metaboliteSet[,c("hmdb","chebi","kegg")], 2, duplicated, incomparables = NA), 1, all)))
      
      
      
      # Remove duplicate metabolites (added by Marten)
      if(nrow(metaboliteSet) >1){
        metaboliteSet <- metaboliteSet[apply(!apply(metaboliteSet[,c("hmdb","chebi","kegg")], 2, duplicated, incomparables = c("character(0)", NA)), 1, all),, drop = FALSE]
      }
      # nMets <- c(nMets,nrow(metaboliteSet))
      
      
      retVal = performMSEA(metaboliteSet = metaboliteSet, 
                           # av_int_and_z_values_matrix = Zint_scores, 
                           av_int_and_z_values_matrix = Zint_pruned,
                           patient = patient, 
                           gene_in, 
                           n_patients, 
                           thresh_F_pos, 
                           thresh_F_neg, 
                           path, 
                           test = 0, 
                           top, 
                           id, 
                           patient_folder = step_folder
      )
      # cat(retVal, "\n")
      p_value = as.numeric(retVal$p.value)
      if (length(p_value)==0){
        p_value=NA
      }
      metSetResult <- rbind(metSetResult, c("p.value"=p_value, "patient"=patient, "metabolite.set"=gene_in))
      
      nMets <- nrow(metaboliteSet)
      tmp <- data.frame("HGNC"=metSetResult[,3],"p.value"=as.numeric(metSetResult[,"p.value"]), "metabolites"=nMets)
      genExcelFileShort(tmp[order(tmp[,"p.value"]),], paste0(path,step_folder,"/MSEA_results.xls"))
      
      # if (!is.null(genes)){
      #   tmp1 = tmp[order(tmp[,"p.value"]),]
      #   
      #   dummy = c(NA,0,0)
      #   for (l in 1:(100-dim(tmp)[1])){
      #     tmp1 = rbind(tmp1,dummy)
      #   }
      #   
      #   overview = rbind(overview, t(tmp1))
      # }
    }
    cat("\n\n")
  }
}