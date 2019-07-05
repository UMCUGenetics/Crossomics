run <- function(entry, outdir, src){
# entry=1
# outdir="./results"
# src="./src"
  
# src="/hpc/shared/dbg_mz/marcel/Crossomics_Build_Metabolite_Set/src"
# outdir="/hpc/shared/dbg_mz/marcel/Crossomics_Build_Metabolite_Set/results"
  
  source(paste(src, "Crossomics/sourceDir.R", sep="/"))
  sourceDir(paste(src, "Crossomics/build_mets_set", sep="/"))
  
  dir.create(outdir, showWarnings = F)
  dir.create(paste(outdir, "mss", sep="/"), showWarnings = F)
  
  # entry=which(HGNC2Uniprot[,"hgnc_symbol"]=="SLC7A7")

  # dim(HGNC2Uniprot)
  # length(unique(HGNC2Uniprot[,"uniprotswissprot"]))
  # length(unique(HGNC2Uniprot[,"hgnc_symbol"]))

  # load(paste(src, "HGNC2Uniprot_prev.RData", sep="/"))
  load(paste(src, "HGNC2Uniprot.RData", sep="/"))
  id_uniprot = HGNC2Uniprot[entry,"uniprotswissprot"]
  hgnc = HGNC2Uniprot[entry,"hgnc_symbol"]
    
  rval = getMetsPathwayCommons(id_uniprot, src)
  
  outdir=paste(outdir, "mss", sep="/")
  if (!is.null(rval$left) | !is.null(rval$right)){
    
    f = paste(outdir, paste(hgnc, "RData", sep="."), sep="/")
    n=1
    while(file.exists(f)){
      n = n+1
      f = unlist(strsplit(f,".",fixed=TRUE))
      f = f[length(f)-1]
      if (length(grep("_",f,fixed=TRUE))>0) f=unlist(strsplit(f,"_",fixed=TRUE))[1]
      f = paste(paste(f,n,sep="_"),"RData",sep=".")
      f = paste(".",f,sep="")
    }

    save(rval, file=f)
  }
}
      
message("Start")
cat("==> reading arguments:\n", sep = "")

cmd_args = commandArgs(trailingOnly = TRUE)

for (arg in cmd_args) cat("  ", arg, "\n", sep="")

run(as.numeric(cmd_args[1]), cmd_args[2], cmd_args[3])

message("Ready")
