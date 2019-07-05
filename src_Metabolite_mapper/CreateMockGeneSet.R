library(rstudioapi)

seed = 313
nr_mocks = 100

code_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
path <- paste0(code_dir,"/../Results/")

# Load dataset containing patients with their disease genes (dataset is called xls_data)
load(paste0(code_dir,"/../Data/Crossomics_DBS_Marten_Training.RData"))
dis_genes <- unique(xls_data$Gene)
dis_genes <- as.vector(unlist(strsplit(dis_genes, split = ";")))

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
set.seed(seed = seed)
genes <- sample(mock_genes, size = nr_mocks)

write.table(genes, file = paste0(path, "/mock_genes_seed",seed,".txt"), row.names = FALSE, col.names = FALSE)