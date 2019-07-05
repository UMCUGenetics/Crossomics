getRandomGeneList <- function(){
  
  require("biomaRt")
  
#   ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
#   gene_map_table = getBM(attributes=c('hgnc_symbol', 'entrezgene', 'ensembl_gene_id'), mart = ensembl)
#   entregene=gene_map_table[!is.na(gene_map_table[,2]),]
#   entregene_hgnc=entregene[-which(entregene[,1]==""),]
#   save(entregene_hgnc, file="./db/allGenesHS.RData")

#   ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
#   gene_map_table = getBM(attributes=c('hgnc_symbol', 'entrezgene', 'ensembl_gene_id'), filters = "chromosome_name", values = "X", mart = ensembl)
#   entregene=gene_map_table[!is.na(gene_map_table[,2]),]
#   entregene_hgnc=entregene[-which(entregene[,1]==""),]
#   save(entregene_hgnc, file="./db/allGenesHS_chrom_X.RData")
  
#   library("XLConnect")
#   genExcelFileShort(entregene_hgnc_X, wbfile="./db/allGenesHS_chrom_X.xls")
  
  # load("./db/allGenesHS_chrom_X.RData")
  # for (i in 1:dim(entregene_hgnc)[1]){
  #   hgnc=as.vector(entregene_hgnc[i,1])
  #   gene=as.vector(entregene_hgnc[i,2])
  #   save(hgnc,gene,file=paste("./src2/crossomics/db/gene.",i,".RData", sep=""))
  # }

  for (i in 1:15){
  
    load("./input/allGenesHS.RData")
    randomSample = entregene_hgnc[sample(1:dim(entregene_hgnc)[1], 1000, replace=F),]
  
    write.table(cbind(as.vector(randomSample[,1]),as.vector(randomSample[,2])),
                file=paste("./db/random_HGNC_",i,".txt",sep =""),sep = "\t",row.names = F,col.names = F, quote = F)
    
    dir.create(paste("./src2/crossomics/db",i,sep="_"), showWarnings = FALSE)
    
    randomSample=entregene_hgnc

    for (j in 1:dim(randomSample)[1]){
      hgnc=as.vector(randomSample[j,1])
      gene=as.vector(randomSample[j,2])
      # save(hgnc,gene,file=paste(paste("./src2/crossomics/db",i,sep="_"),"/gene.",j,".RData", sep=""))
      save(hgnc,gene,file=paste(paste("./db",hgnc,sep="/"), ".RData", sep=""))
    }
  }
}