# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# SessionInfo -------------------------------------------------------------
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# R version:  3.6.1 (2019-07-05)
# platform:   x86_64-apple-darwin15.6.0 (64-bit)
# OS:         macOS Mojave 10.14.6
# 
# libraries:
# tidyr       0.8.3
# rstudioapi  0.10
# stringr     1.4.0
# data.table  1.12.6
# tidyr       1.0.0

# ON HPC:
# R version:  3.6.0 (2019-04-26)
# Platform:   x86_64-pc-linux-gnu (64-bit)
# Running under: CentOS Linux 7 (Core)
# stringr     1.4.0     
# dplyr       0.8.3       
# data.table  1.12.2 
# tidyr       0.8.3      
# tidyselect  0.2.5  
# vctrs       0.2.0       
# crayon      1.3.4      
# backports   1.1.4  
  



# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Libraries and input variables -------------------------------------------
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# This part ensures a correct loading of libraries on my (Marten) HPC environment with my R version (3.6.0). If someone else wants to run this, 
# he/she should point to a location on the HPC with the correct libraries.
if(Sys.getenv("RSTUDIO") != "1") {
  cmd_args <- commandArgs(trailingOnly = TRUE)
  patient_number <- as.numeric(cmd_args[1])
  thresholds <- cmd_args[2]
  max_rxns <- cmd_args[3]
  max_rxns <- as.numeric(unlist(strsplit(max_rxns, split = ",")))
  steps <- cmd_args[4]
  steps <- as.numeric(unlist(strsplit(steps, split = ",")))
  code_dir <- cmd_args[5]
  seed_file <- cmd_args[6]
  seed <- as.integer(sub(".txt", "", sub(".*seed", "", seed_file)))
  R_location <- cmd_args[7]
  R_location <- sub("/bin", "/lib64", R_location)
  
  outdir <- "Results/"
  
  # R_location <- "/hpc/local/CentOS7/dbg_mz/R_libs/3.6.0/lib64" 
  
  suppressMessages(library("stringr",lib.loc = R_location)) # string manipulation, add leading 0's
  suppressMessages(library("backports",lib.loc = R_location))
  suppressMessages(library("crayon",lib.loc = R_location))
  suppressMessages(library("vctrs",lib.loc = R_location))
  suppressMessages(library("tidyselect",lib.loc = R_location))
  suppressMessages(library("tidyr",lib.loc = R_location))
  suppressMessages(library("data.table",lib.loc = R_location))
  suppressMessages(library("dplyr",lib.loc = R_location))
  mock_date <- NULL
  
} else {
  library("stringr") # string manipulation, add leading 0's
  library("rstudioapi")
  library("tidyr")
  library("data.table")
  library("dplyr")
  
  patient_number <- 1
  thresholds <- "-1;1.5,-1.5;2,-3;3,-5;5"
  # thresholds <- "-1;1.5"
  max_rxns <-"8,10,12,15,17,19"
  # max_rxns <-"19"
  max_rxns <- as.numeric(unlist(strsplit(max_rxns, split = ",")))
  steps <- "0,1,2,3,4,5"
  # steps <- "5"
  steps <- as.numeric(unlist(strsplit(steps, split = ",")))
  code_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
  seed <- 72563
  outdir <- "Resultstest/"
  mock_date <- "2019-10-22/"
}

date_input <- "2019-08-12" # The date of the data/mss_0 etc. runs
date_run <- "2019-11-21" # The date of this run 

nr_mocks <- 200



# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Other variables ---------------------------------------------------------
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Subset_Of_Patients <- FALSE

thresh_df <- sapply(strsplit(unlist(strsplit(thresholds, split = ",")), ";"), `[`)

thresh_neg_list <- as.numeric(thresh_df[1,])
thresh_pos_list <- as.numeric(thresh_df[2,])


# Remove any metabolites that (for some reason) should not be present
# HMDB0002467 <- determined to be non-bodily substances, but created in the lab
bad_mets <- c("HMDB0002467")

train_data_name <- "Crossomics_DBS_Marten_TraVal_Inclusion_only_updated20191031.RData"

redo <- NULL #(null or "_REDO")



# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Functions ---------------------------------------------------------------
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

removeMets <- function(metaboliteSet, mets2remove, identifier){
  metaboliteSet <- metaboliteSet[!metaboliteSet[,identifier] %in% mets2remove[,identifier],,drop = FALSE]
  return(metaboliteSet)
}

make.data.table <- function(datatable, identifier){
  datatable[  , .(rxn_id = toupper(paste.unique(rxn_id)),
                  step = paste.unique(step),
                  met_in = paste.unique(met_in),
                  left_right = paste.unique(left_right),
                  met_short = paste.unique(met_short),
                  met_long = paste.unique(met_long),
                  hmdb = toupper(paste.unique(hmdb)),
                  kegg = toupper(paste.unique(kegg)),
                  chebi = paste.unique(chebi),
                  pubchem = paste.unique(pubchem),
                  rxn_name= paste.unique(rxn_name),
                  rxn= paste.unique(rxn),
                  resource= paste.unique(resource),
                  rxn_formula= paste.unique(rxn_formula),
                  path= paste.unique(path)
  ), 
  by = identifier]
}

real_duplicated <- function(set){
  duplicated(set, incomparables = c("character(0)", NA)) | duplicated(set, fromLast = TRUE, incomparables = c("character(0)", NA))
}

paste.unique <- function(colname){
  index <- colname != "character(0)" & colname != "NA" & !is.na(colname)
  if(sum(index) > 1){
    paste(unique(tolower(colname[index])), collapse = ", ")
  } else {
    tolower(colname[index])
  }
}

# To concatenate unique IDs of metabolites in the tibbles
concat_unique <- function(x){paste(unique(x),  collapse=',')}



# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Load data ---------------------------------------------------------------
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

source(paste0(code_dir,"/Supportive/sourceDir.R"))
sourceDir(paste0(code_dir,"/Supportive"), trace = FALSE)


load(paste0(code_dir,"/../Data/", train_data_name))

# Load mock gene set
mss <- read.table(paste0(code_dir,"/../Results/Mock_genes/",mock_date,"mock_genes",nr_mocks,"_seed",seed,".txt"), stringsAsFactors = FALSE)[,1]

# correct naming of new training set to old format
if(sum(colnames(xls_data) == "Patient number in set" | colnames(xls_data) == "Patient.number.in.set") > 0){
  colnames(xls_data)[colnames(xls_data) == "Patient number in set"] <- "Old.patient.number"
  colnames(xls_data)[colnames(xls_data) == "Patient.number.in.set"] <- "Old.patient.number"
}

xls_data$Patient.number <- sapply(strsplit(xls_data$Old.patient.number,"[P.]"), function(x)
  paste0("P", sprintf("%03d",as.numeric(x[2])), ".", x[3])
)

xls_data$Patient <- sapply(strsplit(xls_data$Old.patient.number,"[P.]"), function(x)
  paste0("P", sprintf("%03d",as.numeric(x[2])))
)


# Load patient subset
if (!Subset_Of_Patients){
  xls_data$PatientIDs <- paste(xls_data$Dataset, xls_data$Patient, sep = "^")
  uni_dat_pat <- unique(xls_data$PatientIDs)
} else {
  uni_dat_pat <- read.table(paste0(code_dir,"/../Results/",date_input,"/patient_subset_seed",seed,".txt"), stringsAsFactors = FALSE)[,1]
}


# Set which metabolites to remove
mets2remove <- as.data.frame(readRDS(paste0(code_dir,"/../Data/mets2remove.RDS")))




# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Perform metabolite mapper on patients -----------------------------------
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

i <- uni_dat_pat[patient_number]


# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Get patient and dataset -------------------------------------------------
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

tmp <- unlist(strsplit(i, split = "\\^"))
dataset <- tmp[1]
patient <- tmp[2]
rm(tmp)

cat("start patient:", patient, "dataset:", dataset, "\n")

# Get disease gene for patient
dis_gene <- xls_data$Gene[grepl(patient, xls_data$Patient.number) & xls_data$Dataset == dataset][1]

# In the case there are multiple disease genes stated:
if(grepl("[;,]+", dis_gene)) {
  dis_gene <- unique(trimws(unlist(strsplit(dis_gene, split = "[;,]+"))))
}

# Fix name for MUT / MMUT gene to pathwaycommons version
dis_gene <- gsub(pattern = "^MUT$", replacement = "MMUT", x = dis_gene)



data_location <- xls_data$Location[match(dataset, xls_data$Dataset)]
data_location <- paste(code_dir,"../Data", paste(strsplit(data_location, split = "\\\\")[[1]][c(6:8)], collapse = "/"), sep = "/")




# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Obtain metabolite Z scores ----------------------------------------------
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

old_patient_number <- sub(xls_data[grep(i, xls_data$PatientIDs, fixed = TRUE)[1], Old.patient.number], pattern = "\\..*", replacement = "")

RDS_file <- list.files(path = data_location, pattern = "*.RDS")
Zint_values <- readRDS(file = paste(data_location, RDS_file, sep = "/"))
Zint_values <- as.data.table(Zint_values)
Zint_values[nchar(HMDB_code) == 9, HMDB_code := str_replace(HMDB_code, pattern = "HMDB", replacement = "HMDB00")]

# Remove any metabolites that (for some reason) should not be present
Zint_values <- Zint_values[!HMDB_code %in% bad_mets]

# Collate indistinguishable metabolites
columns_to_compare <- grep("[PC][0-9]+\\.[0-9]", colnames(Zint_values), value = TRUE)
tmp <- as.data.frame(apply(Zint_values[,..columns_to_compare], 1, paste, collapse = ","), stringsAsFactors = FALSE)
rownames(tmp) <- Zint_values$HMDB_code
colnames(tmp) <- "values"
tmp <- aggregate(rownames(tmp), by=tmp['values'], paste, collapse = ",")
rownames(tmp) <- tmp$x
tmp$x <- NULL
tmp <- tidyr::separate(tmp, values, into = columns_to_compare, sep = ",", convert = TRUE)
Zint_pruned <- as.matrix(tmp[order(rownames(tmp)),])
rm(tmp)

#  Get patient specific columns and average Z-scores if #DBS > 1
DBS <- xls_data[tolower(PatientIDs) == tolower(i), Old.patient.number]

if(length(DBS) == 1){
  av_Z_scores <- Zint_pruned[,grep(paste0(DBS,"_Zscore"), colnames(Zint_pruned)), drop = FALSE]
} else {
  Z_cols <- unlist(lapply(paste0(DBS,"_Zscore"), function(x) grep(x, colnames(Zint_pruned))))
  av_Z_scores <- as.matrix(rowMeans(Zint_pruned[,Z_cols]))
  rm(Z_cols)
}
colnames(av_Z_scores) <- paste0("av.z_", patient)



# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Perform MSEA on all distances -------------------------------------------
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Append mock genes with disease gene(s)
mss <- c(mss, dis_gene)

# Check if mss files are saved as RData or RDS and append file name accordingly:
tmp_dir <- paste0(code_dir,"/../Data/",date_input,"/maxrxn",max_rxns[1],"/mss_", steps[1],"_HMDBtranslated")
if (unlist(strsplit(list.files(tmp_dir)[1], split = "\\."))[2] == "RData"){
  save_as <- "RData"
  mss <- paste0(mss,".RData")
} else if (unlist(strsplit(list.files(tmp_dir)[1], split = "\\."))[2] == "RDS"){
  save_as <- "RDS"
  mss <- paste0(mss,".RDS")
}
rm(tmp_dir)

# Go through all parameter combinations
Patient_metSetResult <- data.table()

for (threshold in 1:length(thresh_pos_list)){
  # for (threshold in 1:2){
  thresh_F_pos <- thresh_pos_list[threshold]
  thresh_F_neg <- thresh_neg_list[threshold]
  threshs <- paste(c(thresh_F_pos,thresh_F_neg), collapse = ", ")
  
  for (step in steps){
    # for (step in c(steps[4],steps[5])){
    
    for (maxrxn in max_rxns){
      # for (maxrxn in c(max_rxns[5],max_rxns[6])){
      # cat("Z-value threshold:", thresh_F_neg,"/", thresh_F_pos, "Max rxn:", maxrxn, "step:", step, "\n")
      
      if (date_input >= "2019-08-12"){
        indir <- paste0(code_dir,"/../Data/",date_input,"/maxrxn",maxrxn,"/mss_", step,"_HMDBtranslated")
      } else {
        indir <- paste0(code_dir,"/../Data/",date_input,"_maxrxn",maxrxn,"/mss_", step,"_HMDBtranslated")
      }
      
      metSetResult = NULL
      nMets = NULL    # list with number of metabolites per gene, not used for any calculations, but only for output excel file.
      
      for (j in 1:length(mss)){
        # for (j in 1:70){
        # Skip the gene if there is no metabolite pathway xls_data available, elsewise, load its file
        if (!file.exists(paste(indir, mss[j], sep="/"))) next
        
        # In the new version (12-08-2019) I saved the gene files as RDS files instead of RData
        if(save_as == "RData"){
          load(paste(indir, mss[j], sep="/"))
        } else if (save_as == "RDS"){
          metaboliteSet <- as.matrix(readRDS(paste(indir, mss[j], sep="/")))
        }
        
        # for when a single metabolite is present (it is possibly converted to a character vector)
        if(is.vector(metaboliteSet)) {
          dimnames <- names(metaboliteSet)
          metaboliteSet <- matrix(metaboliteSet, nrow = 1)
          colnames(metaboliteSet) <- dimnames
          rm(dimnames)
        }
        
        #take only (somewhat) interesting columns and remove the rest
        metaboliteSet <- metaboliteSet[,c("rxn_id", "step","met_long","hmdb","kegg","chebi","pubchem")]
        
        gene_in <- strsplit(mss[j], split = "\\.")[[1]][1]
        

        
        
        # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        # Remove metabolites ------------------------------------------------------
        # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        
        # Remove any metabolites that are non-informative
        for(identifier in c("chebi","kegg","pubchem")){
          metaboliteSet <- removeMets(metaboliteSet, mets2remove, identifier)
        }
        
        # Get rid of any rows that don't have an hmdb code and go to next gene if none are left
        index <- which(metaboliteSet[,"hmdb"] == "character(0)")
        if (length(index)>0) metaboliteSet <- metaboliteSet[-index,,drop=FALSE]
        if (nrow(metaboliteSet) == 0) next
        
        # Or that do not appear in the dataset
        index <- which(!as.vector(unlist(lapply(metaboliteSet[,"hmdb"], function(x) length(grep(x, rownames(Zint_pruned))) > 0))))
        if (length(index)>0) metaboliteSet <- metaboliteSet[-index,,drop=FALSE]
        if (nrow(metaboliteSet) == 0) next
        
        # Remove duplicate metabolites
        if(nrow(metaboliteSet) >1){
          tmp <- as.data.frame(metaboliteSet[apply(apply(metaboliteSet[,c("hmdb","chebi","kegg")], 2, real_duplicated), 1, any),])
        }
        
        # Manual fix for the metabolites known under "HMDB0012482, HMDB0002281" (same, but different in the hmdb dataset)
        if(length(grep("HMDB0002281", metaboliteSet[,"hmdb"])) | length(grep("HMDB0012482", metaboliteSet[,"hmdb"]))){
          HMDB_index <- unique(grep("HMDB0002281", metaboliteSet[,"hmdb"]),grep("HMDB0012482", metaboliteSet[,"hmdb"]))
          metaboliteSet[HMDB_index,"hmdb"] <- "HMDB0002281"
        }
        
        
        # Collate rows that have indistinguishable metabolites
        metaboliteSet <- cbind(metaboliteSet, "alt_hmdb" = unlist(lapply(metaboliteSet[,"hmdb"], function(x) grep(x, rownames(Zint_pruned), value = TRUE))))

        metaboliteSet <- as_tibble(metaboliteSet)

        # metaboliteSet <- metaboliteSet %>%
        #   group_by(alt_hmdb) %>%
        #   summarise_each(funs(paste(unique(.), collapse = ",")))
        metaboliteSet <- metaboliteSet %>%
          group_by(alt_hmdb) %>%
          summarise_all(., concat_unique)
        
        metaboliteSet <- as.matrix(metaboliteSet)
        
        
        
        
        # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        # Perform MSEA ------------------------------------------------------------
        # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        
        retVal = performMSEA(metaboliteSet = metaboliteSet, 
                             patient_z_values = av_Z_scores,
                             thresh_F_pos, 
                             thresh_F_neg
        )
        

        p_value <- as.numeric(retVal$p.value)
        if (length(p_value)==0){
          p_value=NA
        }
        mets_exc_thres <- as.numeric(retVal$mets_exc_thres)
        
        nMets <- nrow(metaboliteSet)
        metSetResult <- data.frame(rbind(metSetResult, c("metabolite.set"=gene_in, 
                                                         "p.value"=signif(p_value, digits = 5), 
                                                         "mets in set"=nMets, 
                                                         "mets exc thres"=mets_exc_thres
        )
        ), stringsAsFactors = FALSE)
        
      }
      
      metSetResult <- metSetResult[order(as.numeric(metSetResult[,"p.value"])),]
      
      # Add any disease genes to the results that were completely missed 
      for(gene in dis_gene){
        tmp_index <- which(metSetResult$metabolite.set == gene)
        if(length(tmp_index) == 0){
          metSetResult <-rbind(metSetResult, c("metabolite.set"=gene, 
                                               "p.value" = 1, 
                                               "mets in set" = 0, 
                                               "mets exc thres" = 0)
          )
        }
      }
      
      rank_dis_genes <- which(metSetResult$metabolite.set %in% dis_gene)
      total_genes <- nrow(metSetResult)
      tmp_Patient_metSetResult <- as.data.table(cbind(metSetResult[rank_dis_genes,],
                                                      "Position" = rank_dis_genes,
                                                      "Last_position" = total_genes,
                                                      "DBS" = length(DBS)))
      tmp_Patient_metSetResult[, c("Z_threshold","Step","Max_rxn","Seed","PatientID") := list(
        threshs, step, maxrxn, seed, paste(dataset, patient, sep = "^")
      )]
      
      # Collate patient-disease gene data to 1 data table
      Patient_metSetResult <- rbind(Patient_metSetResult, tmp_Patient_metSetResult)
      rm(tmp_Patient_metSetResult)
    }
  }
}

# Save patient and seed specific disease gene data to 1 file
names(Patient_metSetResult)[1:4] <- c("Gene", "P.value", "Mets_in_set", "Mets_exc_thresh")
dir.create(paste0(code_dir,"/../", outdir, date_run,"/", patient, redo,"_", dataset,"/seed",seed),recursive = TRUE, showWarnings = FALSE)
save(Patient_metSetResult, file = paste0(code_dir,"/../", outdir, date_run,"/", patient, redo,"_", dataset,"/seed",seed, "/MSEA_results.RData"))



