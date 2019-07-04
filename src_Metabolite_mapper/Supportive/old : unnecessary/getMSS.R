getMSS <- function(hgnc, entrezgene, path, src) {
#   entrezgene=gene
#   path=outdir

  require("biomaRt")
  
  outdir=paste(path, "mss", sep="/")
  
  # load(paste(src, "Recon_2.2_biomodels.RData", sep="/"))
  load(paste(src, "Recon2.RData", sep="/"))
  model=recon2$modelR204[,,1]
  
  retVal = findMetabolicEnvironment(entrezgene, model, hgnc, onlyPrimary=FALSE)
###############
#   mss = list.files(path = "./results/mss_P14_primary")
#   for (j in 1:length(mss)){
#     load(paste("./results/mss_P14_primary", mss[j], sep="/"))
#     hgnc=strsplit(mss[j], split=".", fixed =TRUE)[[1]][1]
# ###############

  # retVal$mets=retVal$mets[c(1:14),]
  
    
############################################################################################
  if (!is.null(retVal$mets)){
    
    metsExtend=NULL
    
    if (is.null(dim(retVal$mets))){    
      retVal$mets=t(as.matrix(retVal$mets))
      rownames(retVal$mets)=retVal$mets[,"rxn_id"]
    }
  
    if (dim(retVal$mets)[1]>0){  
      # extend CoA2CarnitineGlycine
#       #############################################################################################
#       for (i in 1:length(retVal$mets[,"met_long"])){
#         message(retVal$mets[i,"met_long"])
#         
#         if (gregexpr(pattern ="coa", retVal$mets[i,"met_short"])[[1]][1]>1) {
#           #
#           consumed = rbind(retVal$mets[i,],retVal$mets[i,])
#           colnames(consumed)[9]="CheBI_id"
#           colnames(consumed)[8]="KEGG_id"
#           colnames(consumed)[10]="PubChem_id"
#           index_consumed = getMets2ReconID(consumed, model)
#           
#           tmp = getPreviousNext(model, index_consumed, produced=TRUE)
#           
#           #metsExtend=rbind(metsExtend,tmp$mets[which(tmp$mets[,"left_right"]=="right"),])
#           metsExtend=rbind(metsExtend,tmp$mets)
#           
#           tmp = getPreviousNext(model, index_consumed, produced=FALSE)
#           
#           #metsExtend=rbind(metsExtend,tmp$mets[which(tmp$mets[,"left_right"]=="right"),])
#           metsExtend=rbind(metsExtend,tmp$mets)
#         }
#       }
# 
#       # only keep carnitine compounds
#       metsExtend = rbind(metsExtend[grep(metsExtend[,"met_long"], pattern = "carnitine", fixed = TRUE),,drop=FALSE],
#                          metsExtend[grep(metsExtend[,"met_long"], pattern = "glycine", fixed = TRUE),,drop=FALSE])
# 
#       metsExtend = rbind(retVal$mets, removeMetsFromSet(metsExtend, model)$result_mets)
#     
#       tmp = metsExtend[,"met_long"]
#       unic = !duplicated(tmp)  ## logical vector of unique values
#       index = seq_along(tmp)[unic]  ## indices
#       retVal$mets = metsExtend[index,]
#       #############################################################################################

        # extend transcription factor    
#       #################################################################################################################  
#       controlledGenes = findControlsExpressionOf(hgnc)
#       controlledGenes = unique(controlledGenes)
#     
#       if (length(controlledGenes)>0){
#         message(paste("Bingo transcription factor found!!!",hgnc))
#         
#         ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
#         gene_map_table2 = getBM(attributes=c('hgnc_symbol', 'entrezgene', 'ensembl_gene_id'),
#                                 filters = 'hgnc_symbol', values = controlledGenes, mart = ensembl)
#         
#         for (k in 1:length(controlledGenes)){
#           
#           if (is.null(retVal$mets) || (dim(retVal$mets)[1]==0)){
#             retVal = findMetabolicEnvironment(gene_map_table2$entrezgene[k], model, gene_map_table2$hgnc_symbol[k])  
#           } else {
#             tmp = findMetabolicEnvironment(gene_map_table2$entrezgene[k], model, gene_map_table2$hgnc_symbol[k])
#             if (!is.null(tmp$mets)){
#               retVal$mets = rbind(retVal$mets, tmp$mets)  
#             }
#           }
#         }
#       }
#       #################################################################################################################  
    }
  }
  
  if (is.null(retVal$mets)) {
    message(paste("Empty set for", hgnc))
  } else {
    if (is.null(dim(retVal$mets))) {
      retVal$mets=t(as.data.frame(retVal$mets))
    }    
    if (dim(retVal$mets)[1]==0) {
      message(paste("Empty set for", hgnc))
    } else {
      save(retVal, file=paste(outdir, paste(hgnc, "RData", sep="."), sep="/"))
    }  
  }
  
# #########  
#   }
# #########  
}