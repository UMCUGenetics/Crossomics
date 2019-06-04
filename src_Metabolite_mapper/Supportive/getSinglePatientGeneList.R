getSinglePatientGeneList <- function(patient){

  library("biomaRt")
  
  path=paste("./src2/crossomics/db_P",patient, sep="")
  dir.create(path)
  
  # # First get metabolite sets ###############################################################################
  for (i in 1:length(patients)){
    
    # gene_list = read.table(paste("./db/pilot_batch1_HGNC.txt",sep=""), header = FALSE, sep="\t")
    gene_list = read.table(paste("./db/P",patient,"_HGNC.txt",sep=""), header = FALSE, sep="\t")
    ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
    gene_map_table = getBM(attributes=c('hgnc_symbol', 'entrezgene'),
                           filters = 'hgnc_symbol', values = as.vector(unlist(gene_list)), mart = ensembl)
    
    for (j in 1:dim(gene_map_table)[1]){
      message(gene_map_table$entrezgene[j])
      message(gene_map_table$hgnc_symbol[j])
      
      hgnc=gene_map_table$hgnc_symbol[j]
      gene=gene_map_table$entrezgene[j]
      save(hgnc,gene,file=paste(path, "/gene.",j,".RData", sep=""))
    }  
  } 
  
}  
