run <- function(entry, outdir, src){
# entry=1
# outdir="./results"
# src="./src"
  
# src="/hpc/shared/dbg_mz/marcel/Crossomics_Build_Metabolite_Set/src"
# outdir="/hpc/shared/dbg_mz/marcel/Crossomics_Build_Metabolite_Set/results"
  
  source(paste(src, "Supportive/sourceDir.R", sep="/"))
  sourceDir(paste(src, "Supportive/", sep="/"))
  
  dir.create(outdir, showWarnings = F)
  dir.create(paste(outdir, "mss", sep="/"), showWarnings = F)
  
  # entry=which(HGNC2Uniprot[,"hgnc_symbol"]=="SLC7A7")

  # dim(HGNC2Uniprot)
  # length(unique(HGNC2Uniprot[,"uniprotswissprot"]))
  # length(unique(HGNC2Uniprot[,"hgnc_symbol"]))

  # load(paste(src, "HGNC2Uniprot_prev.RData", sep="/"))
  # load(paste(src, "HGNC2Uniprot.RData", sep="/"))
  HGNC2Uniprot <- readRDS(paste(outdir, "HGNC2Uniprot20190605.RDS", sep = "/")) # New version, added 06-06-2019
  
  cat("HGNC loaded\n") # Check for where the script is
  
  id_uniprot = HGNC2Uniprot[entry,"uniprotswissprot"]
  hgnc = HGNC2Uniprot[entry,"hgnc_symbol"]
  
  cat("HGNC:", hgnc, "\n") # Check for where the script is
    
  rval = getMetsPathwayCommons(id_uniprot = id_uniprot, src = src)
  
  outdir=paste(outdir, "mss", sep="/")
  if (!is.null(rval$left) | !is.null(rval$right)){
    
    f = paste(outdir, paste0(hgnc, ".RData"), sep="/")
    n=1
    while(file.exists(f)){
      n = n+1
      f = unlist(strsplit(f,".",fixed=TRUE))
      f = f[length(f)-1]
      # if (length(grep("_",f,fixed=TRUE))>0) f=unlist(strsplit(f,"_",fixed=TRUE))[1]
      if (grepl("_(?!.*/)", f, perl =TRUE)) f=unlist(strsplit(f, "_(?=[0-9])", perl=TRUE))[1]
      f = paste0(paste(f,n,sep="_"),".RData")
    }
    save(rval, file=f)
  } else {
    cat("No metabolites/reactions found for", hgnc, "\n")
  }
}
      
message("Start")
cat("==> reading arguments:\n", sep = "")

# New edit (31-05-19), loading all required packages beforehand:{
library(SPARQL)
library(rJava)
library(XML)
library(paxtoolsr)
library(RCurl)
library(bitops)
# }

cmd_args = commandArgs(trailingOnly = TRUE)

# TEST code{
# cmd_args <- c("10736", "/Users/mkerkho7/DIMS2_repo/Crossomics/Data", "/Users/mkerkho7/DIMS2_repo/Crossomics/src_Metabolite_Set_Creation")
#}

for (arg in cmd_args) cat("  ", arg, "\n", sep="")

run(entry = as.numeric(cmd_args[1]), outdir = cmd_args[2], src = cmd_args[3])

message("Ready")




