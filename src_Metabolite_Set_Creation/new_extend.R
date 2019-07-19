run <- function(entry, outdir, src, fr_max_rxns = NULL, max_rxns = NULL){

  source(paste0(src,"/Supportive/sourceDir.R"))
  sourceDir(paste0(src,"/Supportive"))
  

  load(paste(src, "../Data/Recon3D.RData", sep="/"))
  model=Recon3D$Recon3D[,,1]
  # Fix CHEBI: naming, where some Chebi codes are noted like this: "CHEBI:00001" and some like this: "00001"
  model$metCHEBIID <- rapply(model$metCHEBIID, gsub, pattern = "CHEBI:", replacement = "", ignore.case=TRUE, how = "list")
  shortmodel <- model
  rownames(shortmodel$S) <- as.vector(unlist(shortmodel$metNames))
  shortmodel <- Matrix.utils::aggregate.Matrix(abs(shortmodel$S), rownames(shortmodel$S))
  
  # Determine maximum number of reactions before a metabolite is filtered out
  if(!is.null(fr_max_rxns)){
    max_rxns <- floor(fr_max_rxns*ncol(model$S))
  } else if (is.null(max_rxns)){
    warning("no maximum number of reactions supplied, defaulting to 10")
    max_rxns <- 10
  }
  
  recon2chebi=NULL

  
  dir.create(paste0(outdir, "/mss_0"), recursive = TRUE, showWarnings = FALSE)
  dir.create(paste0(outdir, "/mss_1"), recursive = TRUE, showWarnings = FALSE)
  dir.create(paste0(outdir, "/mss_2"), recursive = TRUE, showWarnings = FALSE)
  dir.create(paste0(outdir, "/mss_3"), recursive = TRUE, showWarnings = FALSE)
  dir.create(paste0(outdir, "/mss_4"), recursive = TRUE, showWarnings = FALSE)
  

  for(gene_file in entry){
    if(!file.exists(paste(src, "../Data/mss", gene_file, sep="/"))) next
    load(paste(src, "../Data/mss", gene_file, sep="/"))
    
    if (!is.null(rval)){
      hgnc = unlist(strsplit(gene_file, split = ".", fixed = TRUE))[1]
      cat("extending", hgnc,"\n")
      
      findMetabolicEnvironmentLocal(gene_in = hgnc, 
                                    model = model, 
                                    shortmodel = shortmodel, 
                                    recon2chebi = recon2chebi, 
                                    src = src, 
                                    outdir = outdir, 
                                    rval = rval,
                                    max_rxns = max_rxns
                                    )
    }
  }
}

message("Start")


library("rstudioapi")
library("Matrix.utils")
code_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)

# TEST CODE
load(paste0(code_dir,"/../Data/Crossomics_DBS_Marten_Training.RData"))
path <- paste0(code_dir,"/../Results")
seed <- 313
dis_genes <- as.vector(unlist(strsplit(unique(xls_data$Gene), split = "; ")))
mock_genes <- scan(file = paste0(path, "/mock_genes_seed",seed,".txt"), what = "character", quiet = TRUE)
genes <- c(mock_genes, dis_genes)
mss <- paste(genes,"RData",sep=".")
if(length(grep("MUT.RData", mss))>0) {mss[grep("MUT.RData", mss)] <- "MMUT.RData"}

today <- as.character(Sys.Date())

for(react in c(6, 8, 10, 12, 15)){
  message(Sys.time())
  cat("max # reactions:", react, "\n")
  run(entry = mss[117], 
      outdir = paste0("/Users/mkerkho7/DIMS2_repo/Crossomics/Data/",today, "_maxrxn", react),
      src = "/Users/mkerkho7/DIMS2_repo/Crossomics/src_Metabolite_Set_Creation",
      max_rxns = react)
  message(Sys.time())
}



# run(entry = as.numeric(cmd_args[1]), outdir = cmd_args[2], src = code_dir)

message("Ready")
