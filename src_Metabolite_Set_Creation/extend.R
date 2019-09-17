run <- function(entry, outdir, src, max_rxns, steps){

  source(paste0(src,"/Supportive/sourceDir.R"))
  sourceDir(paste0(src,"/Supportive"))
  

  load(paste(src, "../Data/Recon3D.RData", sep="/"))
  model=Recon3D$Recon3D[,,1]
  # Fix CHEBI: naming, where some Chebi codes are noted like this: "CHEBI:00001" and some like this: "00001"
  model$metCHEBIID <- rapply(model$metCHEBIID, gsub, pattern = "CHEBI:", replacement = "", ignore.case=TRUE, how = "list")
  shortmodel <- model
  rownames(shortmodel$S) <- as.vector(unlist(shortmodel$metNames))
  shortmodel <- Matrix.utils::aggregate.Matrix(abs(shortmodel$S), rownames(shortmodel$S))

  
  recon2chebi=NULL

  for(max_rxn in max_rxns){
    for(step in 0:steps)
    dir.create(paste0(outdir,"/maxrxn",max_rxn,"/mss_",step), recursive = TRUE, showWarnings = FALSE)
  }
  # dir.create(paste0(outdir, "/mss_0"), recursive = TRUE, showWarnings = FALSE)
  # dir.create(paste0(outdir, "/mss_1"), recursive = TRUE, showWarnings = FALSE)
  # dir.create(paste0(outdir, "/mss_2"), recursive = TRUE, showWarnings = FALSE)
  # dir.create(paste0(outdir, "/mss_3"), recursive = TRUE, showWarnings = FALSE)
  # dir.create(paste0(outdir, "/mss_4"), recursive = TRUE, showWarnings = FALSE)
  # dir.create(paste0(outdir, "/mss_5"), recursive = TRUE, showWarnings = FALSE)
  files = list.files(paste(src, "../Data/mss", sep="/"))
  
  gene_file <- files[entry]

  # for(gene_file in entry){
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
                                    max_rxns = max_rxns,
                                    steps = steps)
    }
  # }
}

message("Start")
today <- "2019-08-12"
cmd_args <- commandArgs(trailingOnly = TRUE)
R_libs <- paste0(cmd_args[4], "/lib64/")

print(R_libs)
if(Sys.getenv("RSTUDIO") != 1){
  library("Matrix", lib.loc = R_libs)
  library("Matrix.utils", lib.loc = R_libs)
}





if(Sys.getenv("RSTUDIO") == "1") {
  library("Matrix.utils")
  library("rstudioapi")
  code_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
  
  # TEST CODE
  load(paste0(code_dir,"/../Data/Crossomics_DBS_Marten_Training.RData"))
  path <- paste0(code_dir,"/../Results")
  
  cmd_args[1] <- 35
  cmd_args[2] <- paste0("/Users/mkerkho7/DIMS2_repo/Crossomics/Data/",today)
  cmd_args[3] <- "/Users/mkerkho7/DIMS2_repo/Crossomics/src_Metabolite_Set_Creation"
  # cmd_args[4] <- c(8, 10, 12, 15, 17, 19)
  # cmd_args[5] <- 5
}

# today <- as.character(Sys.Date())

# Supply 1. gene_entry number, 2. outdirectory, 3. source code directory, 4. max reactions

# for(react in c(8, 10, 12, 15)){
#   cat("max # reactions:", react, "\n")
#   # run(entry = mss[117], 
#   #     outdir = paste0("/Users/mkerkho7/DIMS2_repo/Crossomics/Data/",today, "_maxrxn", react),
#   #     src = "/Users/mkerkho7/DIMS2_repo/Crossomics/src_Metabolite_Set_Creation",
#   #     max_rxns = react)
#   run(entry = as.numeric(cmd_args[1]), 
#       outdir = paste0("/Users/mkerkho7/DIMS2_repo/Crossomics/Data/",today, "_maxrxn", react), 
#       src = "/Users/mkerkho7/DIMS2_repo/Crossomics/src_Metabolite_Set_Creation",
#       max_rxns = react)
#   message(Sys.time())
# }

run(entry = as.numeric(cmd_args[1]), 
    outdir = cmd_args[2],
    # outdir = paste0("/Users/mkerkho7/DIMS2_repo/Crossomics/Data/",today), 
    src = cmd_args[3],
    # src = "/Users/mkerkho7/DIMS2_repo/Crossomics/src_Metabolite_Set_Creation",
    max_rxns = c(8, 10, 12, 15, 17, 19),
    steps = 5)



# run(entry = as.numeric(cmd_args[1]), outdir = cmd_args[2], src = code_dir)

message("Ready")
