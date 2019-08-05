# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# SessionInfo -------------------------------------------------------------
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# R version:  3.6.0 (2019-04-26)
# platform:   x86_64-apple-darwin15.6.0 (64-bit)
# OS:         macOS Mojave 10.14.6
# 
# libraries:
# heatmap3    1.1.6
# tidyr       0.8.3
# rstudioapi  0.10
# Cairo       1.5-10
# stringr     1.4.0
# data.table  1.12.2



# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Libraries and input variables -------------------------------------------
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# This part ensures a correct loading of libraries on my (Marten) HPC environment with my R version (3.6.0). If someone else wants to run this, 
# he/she should point to a location on the HPC with the correct libraries.
if(Sys.getenv("RSTUDIO") != "1") {
  R_location <- "/hpc/local/CentOS7/dbg_mz/R_libs/3.6.0/lib64" 
  
  suppressMessages(library("stringr",lib.loc = R_location)) # string manipulation, add leading 0's
  suppressMessages(library("Cairo",lib.loc = R_location))
  suppressMessages(library("backports",lib.loc = R_location))
  suppressMessages(library("crayon",lib.loc = R_location))
  suppressMessages(library("vctrs",lib.loc = R_location))
  suppressMessages(library("tidyselect",lib.loc = R_location))
  suppressMessages(library("tidyr",lib.loc = R_location))
  suppressMessages(library("heatmap3",lib.loc = R_location))
  suppressMessages(library("data.table",lib.loc = R_location))
  suppressMessages(library("dplyr",lib.loc = R_location))
  
  # Supply the maxrxn and the directory where this script is placed
  cmd_args <- commandArgs(trailingOnly = TRUE)
  threshold <- as.numeric(cmd_args[1])
  maxrxn <- as.numeric(cmd_args[2])
  step <- as.numeric(cmd_args[3])
  patient_number <- as.numeric(cmd_args[4])
  code_dir <- cmd_args[5]
  seed <- cmd_args[6]
} else {
  library("stringr") # string manipulation, add leading 0's
  library("Cairo")
  library("rstudioapi")
  library("tidyr")
  library("heatmap3")
  library("data.table")
  threshold <- 4 # possible 1, 2, 3, 4
  maxrxn <- 15 # possible: 8, 10, 12, 15
  step <- 4 # possible 1, 2, 3, 4
  patient_number <- 29 # possible 1:51
  code_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
  seed <- 6734892
}


# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Other variables ---------------------------------------------------------
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Subset_Of_Patients <- FALSE

thresh_pos_list <- c(1,1.5,2,3)
thresh_neg_list <- c(-0.5,-1,-1.5,-3)

if(Subset_Of_Patients){
  sample_number <- 4
  Specific_Patients <- c("P64","P65")
} 

top <- 20
id <- "hmdb"
date_input <- "2019-07-19" # The date of the data/mss_0 etc. runs
date_run <- "2019-08-02"
nr_mocks <- 800

# Remove any metabolites that (for some reason) should not be present
# HMDB0002467 <- determined to be non-bodily substances, but created in the lab
bad_mets <- c("HMDB0002467")



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




# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Load data ---------------------------------------------------------------
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

setwd(paste0(code_dir,"/../"))

source(paste0(code_dir,"/Supportive/sourceDir.R"))
sourceDir(paste0(code_dir,"/Supportive"), trace = FALSE)
outdir <- "Results/"

load("Data/Crossomics_DBS_Marten_Training.RData")

# Load mock gene set
mss <- read.table(paste0("./Results/",date_input,"/mock_genes",nr_mocks,"_seed",seed,".txt"), stringsAsFactors = FALSE)[,1]

# Load patient subset
uni_dat_pat <- read.table(paste0("./Results/",date_input,"/patient_subset_seed",seed,".txt"), stringsAsFactors = FALSE)[,1]

dat_pat <- paste(xls_data$Dataset, xls_data$Patient.number, sep = "^")
# uni_dat_pat <- unique(sapply(strsplit(dat_pat, split = "\\."), `[`, 1))
if(Subset_Of_Patients){
  if(length(Specific_Patients) > 0 ) {
    uni_dat_pat <- uni_dat_pat[unlist(lapply(Specific_Patients, function(x) grep(paste0(x,'$'), uni_dat_pat)))]
  } else {
    set.seed(seed = seed)
    uni_dat_pat <- sample(uni_dat_pat, sample_number)
  }
}

# Set which metabolites to remove
mets2remove <- as.data.frame(readRDS("Data/mets2remove.RDS"))




# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Perform metabolite mapper on patients -----------------------------------
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# for (i in uni_dat_pat){
  i <- uni_dat_pat[patient_number]
  
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
  
  # Fix incorrect name for MMUT gene
  if(dis_gene == "MUT") dis_gene <- "MMUT"
  
  # In the case there are multiple disease genes stated:
  if(grepl(";",dis_gene)){
    dis_gene <- unlist(strsplit(dis_gene, split = "; "))
    dis_gene <- trimws(dis_gene)
  }
  mss <- c(mss, dis_gene)
  mss <- paste0(mss,".RData")
  
  
  data_location <- xls_data$Location[match(dataset, xls_data$Dataset)]
  # Make dataset location mac-compatible
  # data_location <- gsub("Y:", "/Volumes/Metab", data_location)
  # data_location <- gsub("\\", "/", data_location, fixed = TRUE)
  data_location <- paste("Data", paste(strsplit(data_location, split = "\\\\")[[1]][c(6:8)], collapse = "/"), sep = "/")
  

  
  
  # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # Calculate Z scores ------------------------------------------------------
  # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  if(nchar(patient) == 2) {
    tmp_patient <- unlist(strsplit(patient, split = ""))
    patient <- paste0(tmp_patient[1],"0",tmp_patient[2])
  }
  
  Zint_values <- generate_av_Z_scores(patient = patient, data_location = data_location)
  rownames(Zint_values)[nchar(rownames(Zint_values)) == 9] <- str_replace(rownames(Zint_values)[nchar(rownames(Zint_values)) == 9], pattern = "HMDB", replacement = "HMDB00")
  
  # Remove any metabolites that (for some reason) should not be present
  Zint_values <- Zint_values[-grep(bad_mets, rownames(Zint_values)),]
  
  # Collate all identical rows (indistinguishable metabolites)
  tmp <- as.data.frame(apply(Zint_values, 1, paste, collapse = ","))
  colnames(tmp) <- "values"
  tmp <- aggregate(rownames(tmp), by=tmp['values'], paste, collapse = ",")
  rownames(tmp) <- tmp$x
  tmp$x <- NULL
  tmp <- tidyr::separate(tmp, values, into = colnames(Zint_values), sep = ",", convert = TRUE)
  Zint_pruned <- as.matrix(tmp[order(rownames(tmp)),])
  rm(tmp)
  
  # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # Perform MSEA on all distances -------------------------------------------
  # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  # # Prepare mock gene set ---------------------------------------------------
  # genes <- NULL
  # if (!file.exists(paste0("./db/",patient,"_HGNC.txt"))){
  #   # # For getting random mock genes
  #   # mock_genes <- read.table(file = paste0(code_dir,"/../Data/All_Genes_Ensembl_apr_2019_GRCh38p12_extended.txt"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  #   # # mock_genes <- read.table(file = "/Users/mkerkho7/DIMS2_repo/Crossomics/Data/All_Genes_Ensembl_apr_2019_GRCh38p12_extended.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  #   # mock_genes <- mock_genes[mock_genes$Gene.type == "protein_coding",]
  #   # mock_genes <- mock_genes[!mock_genes$HGNC.ID == "",]
  #   # mock_genes <- mock_genes[!duplicated(mock_genes$Gene.name),]
  #   # mock_genes <- mock_genes[!grepl("orf", mock_genes$Gene.name),]
  #   # mock_genes <- mock_genes$Gene.name
  #   # 
  #   # # make sure that the disease gene isn't included twice by removing it from the mock_genes list
  #   # mock_genes <- mock_genes[mock_genes != dis_gene]
  #   # set.seed(seed = seed)
  #   # genes <- sample(mock_genes, size = nr_mocks)
  #   mock_genes <- scan(file = paste0(outdir, "/mock_genes_seed",seed,".txt"), what = "character", quiet = TRUE)
  #   genes <- c(mock_genes, dis_gene)
  #   mss <- paste(genes,"RData",sep=".")
  #   
  #   # save mock genes + real gene
  #   patient_folder <- paste0(date,"/", patient, "_", dataset)
  #   dir.create(paste0(outdir,"/",patient_folder), showWarnings = FALSE, recursive = TRUE)
  #   
  #   # if(save_mock) write.table(genes, file = paste0(outdir, "/", patient_folder, "/mock_genes_seed",seed,".txt"), row.names = FALSE, col.names = FALSE)
  # } else {
  #   # Real WES ----------------------------------------------------------------
  #   mss = read.table(paste0("./db/",patient,"_HGNC.txt"), header = FALSE, sep="\t")
  #   mss = as.vector(unlist(mss))
  #   mss = paste(mss,"RData",sep=".")
  # }
  
  # Any file location is fine in this format as all variable combinations have the same file(names) in them
  # mss <- list.files("Data/2019-07-19_maxrxn8/mss_0_HMDBtranslated")
  
  # for (threshold in 1:length(thresh_pos_list)){
    thresh_F_pos <- thresh_pos_list[threshold]
    thresh_F_neg <- thresh_neg_list[threshold]
    cat("Z-value threshold:", thresh_F_neg,"/", thresh_F_pos, "\n")
    
    # for (step in steps){
      cat("step:", step, "\n")
      overview <- NULL # at the end
      
      indir <- paste0("Data/",date_input,"_maxrxn",maxrxn,"/mss_", step,"_HMDBtranslated")
      # patient_folder <- paste0(date_run,"/", patient, "_", dataset,"/seed",seed)
      
      # step_folder <- paste0(patient_folder,"/maxrxn",maxrxn,"_thresh_n",thresh_F_neg,"_p",thresh_F_pos,"_step_", step)
      step_folder <- paste0(date_run,"/", patient, "_", dataset,"/seed",seed,"/maxrxn",maxrxn,"_thresh_n",thresh_F_neg,"_p",thresh_F_pos,"_step_", step)
      
      dir.create(paste0(outdir,step_folder), recursive = TRUE, showWarnings = FALSE)
      
      metSetResult = NULL
      nMets = NULL    # list with number of metabolites per gene, not used for any calculations, but only for output excel file.
      
      # for (j in 1:length(mss)){
        for (j in 1:length(mss)){
        # if(j%%5 == 0) cat(paste0(j, "%... "))
        if(j%%100 == 0) cat("gene:", mss[j], "number:", j,"\n")
        # Skip the gene if there is no metabolite pathway xls_data available, elsewise, load its file
        if (!file.exists(paste(indir, mss[j], sep="/"))) next
        
        load(paste(indir, mss[j], sep="/"))
        gene_in <- strsplit(mss[j], split = "\\.")[[1]][1]
        
        # The complete metabolitesets like they are loaded here are already named 'metaboliteSet'
        # result_mets_x is the name of how the RData file is loaded in and depends on the step size.
        # if (step==1){
        #   metaboliteSet <- result_mets_1
        # } else if (step==2){
        #   metaboliteSet <- result_mets_2
        # } else if (step==3){
        #   metaboliteSet <- result_mets_3
        # } else if (step==4){
        #   metaboliteSet <- result_mets_4
        # } else {
        #   metaboliteSet <- result_mets_0
        # }
        
        if(is.vector(metaboliteSet)) { # for when a single metabolite is present (it is possibly converted to a character vector)
          dimnames <- names(metaboliteSet)
          metaboliteSet <- matrix(metaboliteSet, nrow = 1)
          colnames(metaboliteSet) <- dimnames
          rm(dimnames)
        }
        
        
        
        
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

        # Remove duplicate metabolites, but collate their data if they are known under different ID's/names
        if(nrow(metaboliteSet) >1){
          tmp <- as.data.frame(metaboliteSet[apply(apply(metaboliteSet[,c("hmdb","chebi","kegg")], 2, real_duplicated), 1, any),])
          setDT(tmp)
          tmp <- make.data.table(tmp, "hmdb")
          tmp <- make.data.table(tmp, "chebi")
          metaboliteSet <- metaboliteSet[apply(!apply(metaboliteSet[,c("hmdb","chebi","kegg")], 2, real_duplicated), 1, all),]
          metaboliteSet <- rbind(metaboliteSet, as.matrix(tmp[,-1]))
          rm(tmp)
        }
        
        # nMets <- c(nMets,nrow(metaboliteSet))
        
        # Manual fix for the metabolites known under "HMDB0012482, HMDB0002281" (same, but different in the hmdb dataset)
        if(length(grep("HMDB0002281", metaboliteSet[,"hmdb"])) | length(grep("HMDB0012482", metaboliteSet[,"hmdb"]))){
          HMDB_index <- unique(grep("HMDB0002281", metaboliteSet[,"hmdb"]),grep("HMDB0012482", metaboliteSet[,"hmdb"]))
          metaboliteSet[HMDB_index,"hmdb"] <- "HMDB0002281"
        }
        
        retVal = performMSEA(metaboliteSet = metaboliteSet, 
                             # av_int_and_z_values_matrix = Zint_values, 
                             av_int_and_z_values_matrix = Zint_pruned,
                             patient = patient, 
                             gene_in, 
                             thresh_F_pos, 
                             thresh_F_neg, 
                             path = outdir, 
                             top, 
                             id, 
                             patient_folder = step_folder,
                             plot = FALSE
        )
        # cat(retVal, "\n")
        p_value = as.numeric(retVal$p.value)
        if (length(p_value)==0){
          p_value=NA
        }
        # metSetResult <- rbind(metSetResult, c("p.value"=p_value, "patient"=patient, "metabolite.set"=gene_in))
        # 
        # nMets <- nrow(metaboliteSet)
        # tmp <- data.frame("HGNC"=metSetResult[,3],"p.value"=as.numeric(metSetResult[,"p.value"]), "metabolites"=nMets)
        nMets <- nrow(metaboliteSet)
        metSetResult <- data.frame(rbind(metSetResult, c("metabolite.set"=gene_in, 
                                                         "p.value"=signif(p_value, digits = 5), 
                                                         "mets in set"=nMets, 
                                                         "mets exc thres"=retVal$mets_exc_thres
                                                         )
                                         ), stringsAsFactors = FALSE)
        
      }
      
      # genExcelFileShort(list = as.data.frame(metSetResult[order(as.numeric(metSetResult[,"p.value"])),]), 
      #                   wbfile = paste0(outdir,"/", step_folder,"/MSEA_results.xls"))
      metSetResult <- metSetResult[order(as.numeric(metSetResult[,"p.value"])),]
      save(metSetResult, file = paste0(outdir, step_folder,"/MSEA_results.RData"))
      
      cat("\n\n")
    # }
  # }
# }
