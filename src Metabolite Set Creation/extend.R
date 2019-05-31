run <- function(entry, outdir, src){
  # entry=1
  # outdir="./results"
  # src="./src"
  
  # src="/hpc/shared/dbg_mz/marcel/Crossomics_Build_Metabolite_Set/src"
  # outdir="/hpc/shared/dbg_mz/marcel/Crossomics_Build_Metabolite_Set/results"
  
  reconVersion=2.0

  source(paste(src, "Crossomics/sourceDir.R", sep="/"))
  sourceDir(paste(src, "Crossomics/build_mets_set", sep="/"))
  sourceDir(paste(src, "Crossomics", sep="/"))
  
  if (reconVersion==2.2){
    
    load(paste(src, "Recon_2.2_biomodels.RData", sep="/"))
    load(paste(src, "recon2chebi_MODEL1603150001.RData", sep="/"))
    rownames(recon2chebi)=recon2chebi[,1]
    recon2chebi=recon2chebi[model@met_id,]
    
  } else if (reconVersion==2.0){
    
    load(paste(src, "Recon2.RData", sep="/"))
    model=recon2$modelR204[,,1]
    recon2chebi=NULL
  }
  
  message(dim(model$S))
  
  dir.create(paste(src, "../results/mss_0", sep="/"),showWarnings = FALSE)
  dir.create(paste(src, "../results/mss_1", sep="/"),showWarnings = FALSE)
  dir.create(paste(src, "../results/mss_2", sep="/"),showWarnings = FALSE)
  dir.create(paste(src, "../results/mss_3", sep="/"),showWarnings = FALSE)
  dir.create(paste(src, "../results/mss_4", sep="/"),showWarnings = FALSE)
  
  # load(paste(src, "../results/metabolite_sets_step_0,1,2,3_1.0_filter_1.1/mss_1", files[i], sep="/"))
  # 
  # myFirstDim=NULL
  # files = list.files(paste(src, "../results/metabolite_sets_step_0,1,2,3_1.0_filter_1.1/mss_0", sep="/"))
  # for (i in 1:length(files)) {
  #   load(paste(src, "../results/metabolite_sets_step_0,1,2,3_1.0_filter_1.1/mss_0", files[i], sep="/"))
  #   # message(paste("i",i))
  #   # message(dim(result_mets_0)[1])
  #   if (is.null(dim(result_mets_0)[1])){
  #     message(i)
  #     myFirstDim=c(myFirstDim,1)
  #   }
  #   myFirstDim=c(myFirstDim,dim(result_mets_0)[1])
  # }

  files = list.files(paste(src, "../results/metabolite_sets_step_0,1,2,3_1.0_filter_1.1/mss_WG_step_0", sep="/"))
  # files = list.files(paste(src, "../results/mss_PathwayCommons", sep="/"))
  # SLC7A7: entry=18425
  # grep("SLC7A7", files)
  
  message(paste(src, "../results/mss_PathwayCommons"))
  message(files[entry])
  message(paste(src, "../results/mss_PathwayCommons", files[entry], sep="/"))

  # load(paste(src, "../results/mss_PathwayCommons", files[entry], sep="/"))
  load(paste(src, "../results/metabolite_sets_step_0,1,2,3_1.0_filter_1.1/mss_WG_step_0", files[entry], sep="/"))
  
  if (!is.null(rval)){
    hgnc = unlist(strsplit(files[entry], split = ".", fixed = TRUE))[1]
      
    findMetabolicEnvironmentLocal(hgnc, model, recon2chebi, src, rval)
  }
  
}

message("Start")
cat("==> reading arguments:\n", sep = "")

cmd_args = commandArgs(trailingOnly = TRUE)

for (arg in cmd_args) cat("  ", arg, "\n", sep="")

run(as.numeric(cmd_args[1]), cmd_args[2], cmd_args[3])

message("Ready")
