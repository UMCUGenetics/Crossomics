# test to see if any of the disease genes is in one of the mock gene sets
nr_mocks <- 200

code_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
mock_dir <- paste0(code_dir, "/../Results/Mock_genes/2019-12-10")

mock_files <- list.files(mock_dir)

# load(paste0("/Users/mkerkho7/DIMS2_repo/Crossomics/Data/Crossomics_DBS_Marten_Training_Validation.RData"))
load(paste0(code_dir,"/../Data/Crossomics_DBS_Marten_trimmed20191205.RData"))
dis_genes <- unique(xls_data$Gene)
dis_genes <- unique(trimws(unlist(strsplit(xls_data$Gene, split = "[;,]+"))))
# dis_genes[dis_genes == "MUT"] <- "MMUT" # Fix incorrect MUT name


for(i in mock_files){
  mss <- read.table(paste(mock_dir, i, sep = "/"), stringsAsFactors = FALSE)[,1]
  if(any(dis_genes %in% mss)){
    cat("disease gene in seed:", i, "(",dis_genes[dis_genes %in% mss],")\n")
  } else {
    cat("no disease genes in seed:", i,"\n")
  }
}
