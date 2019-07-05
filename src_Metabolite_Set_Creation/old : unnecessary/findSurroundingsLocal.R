
reconVersion=2.0
path="./results"
src="./src"

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

dir.create(paste(getwd(),src, "../results/mss_0", sep="/"),showWarnings = FALSE)
dir.create(paste(getwd(),src, "../results/mss_1", sep="/"),showWarnings = FALSE)
dir.create(paste(getwd(),src, "../results/mss_2", sep="/"),showWarnings = FALSE)
dir.create(paste(getwd(),src, "../results/mss_3", sep="/"),showWarnings = FALSE)

files = list.files("./results/mss_PathwayCommons")
# HEXA i=1278
# grep("HEXA", files)

# for (i in 1:length(files)){
for (i in 1:10){
    load(paste("./results/mss_PathwayCommons", files[i], sep="/"))

  if (!is.null(rval)){
    hgnc = unlist(strsplit(files[i], split = ".", fixed = TRUE))[1]

    # outdir=paste(path, "mss_0", sep="/")
    # step = 0
    findMetabolicEnvironmentLocal(hgnc, model, recon2chebi, src, rval)
  }
}
# 
library("XLConnect")

hgnc="NGLY1"
step="3"
load(paste0("./results/metabolite_sets_step_0,1,2,3_1.1_filter_1.1/mss_", step, "/", hgnc, ".RData"))
genExcelFileShort(as.data.frame(result_mets_3), paste0("./results/MetsSets/", hgnc, "_step_", step, ".xls"))
#   


