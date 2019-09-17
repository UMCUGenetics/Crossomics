# test to see if any of the disease genes is in one of the mock gene sets
seeds <- c(8372, 2528, 6140, 3880, 2771, 
           8455, 3200, 6250, 4860, 6297, 
           244, 3764, 2464, 3218, 2282, 
           5600, 2359, 8353, 6399, 2001)
nr_mocks <- 100

code_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)

# load(paste0("/Users/mkerkho7/DIMS2_repo/Crossomics/Data/Crossomics_DBS_Marten_Training_Validation.RData"))
load(paste0(code_dir,"/../Data/Crossomics_DBS_Marten_Training_Validation.RData"))
dis_genes <- unique(xls_data$Gene)
dis_genes <- unique(trimws(unlist(strsplit(xls_data$Gene, split = "[;,]+"))))
dis_genes[dis_genes == "MUT"] <- "MMUT" # Fix incorrect MUT name


for(i in 1:length(seeds)){
  seed <- seeds[i]
  mss <- read.table(paste0("./Results/Mock_genes/mock_genes",nr_mocks,"_seed",seed,".txt"), stringsAsFactors = FALSE)[,1]
  if(any(dis_genes %in% mss)){
    cat("disease gene in seed:", seed, "(",dis_genes[dis_genes %in% mss],")\n")
  } else {
    cat("no disease genes in seed:", seed,"\n")
  }
}
