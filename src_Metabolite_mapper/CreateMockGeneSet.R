library(rstudioapi)

seeds <- c(8372, 2528, 6140, 3880, 2771, 
           8455, 3200, 6250, 4860, 6297, 
           244, 3764, 2464, 3218, 2282, 
           5600, 2359, 8353, 6399, 2001)
nr_mocks <- 100

code_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
path <- paste0(code_dir,"/../Results/")

# Load dataset containing patients with their disease genes (dataset is called xls_data)
load(paste0(code_dir,"/../Data/Crossomics_DBS_Marten_Training.RData"))
dis_genes <- unique(xls_data$Gene)
dis_genes[dis_genes == "MUT"] <- "MMUT" # Fix incorrect MUT name
dis_genes <- as.vector(unlist(strsplit(dis_genes, split = ";")))
dis_genes <- trimws(dis_genes, which = "both")

# Get complete set of human genes from ensembl website (including HGNC.ID code, Gene type and Gene name)
mock_genes <- read.table(file = paste0(code_dir,"/../Data/All_Genes_Ensembl_apr_2019_GRCh38p12_extended.txt"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
mock_genes <- mock_genes[mock_genes$Gene.type == "protein_coding",]
mock_genes <- mock_genes[!mock_genes$HGNC.ID == "",]
mock_genes <- mock_genes[!duplicated(mock_genes$Gene.name),]
mock_genes <- mock_genes[!grepl("orf", mock_genes$Gene.name),]
mock_genes <- mock_genes$Gene.name

# make sure that none of the disease genes is included in the mock gene list
# mock_genes <- mock_genes[mock_genes != dis_genes]
mock_genes <- mock_genes[!mock_genes %in% dis_genes]

for (seed in seeds){
  set.seed(seed = seed)
  genes <- sample(mock_genes, size = nr_mocks)
  
  write.table(genes, file = paste0(path, "/mock_genes",nr_mocks,"_seed",seed,".txt"), row.names = FALSE, col.names = FALSE)
}