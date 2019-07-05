# From gene to metabolite
getMetsRecon2 <- function(entrezgene, level, met){
#     entrezgene=59272
#     level=0
#     met=""
  
  gene_id = rapply(strsplit(genes, split = ".", fixed = T), function(x) head(x, 1))
  index=which(gene_id==as.character(entrezgene))
  if (length(index)==0) {
    message("gene not found in recon2") 
    return()
  } else {
    gene_in=genes[index]  
  } 
  
  if (length(index)>1) {
    rxn_index = rxnGeneMat[,index[2]] # <============================ To do 2 versions exist!! only second one works!!
  } else {
    rxn_index = rxnGeneMat[,index]
  }
  rxn_index = which(rxn_index!=0)
  
  S.sub = S[,rxn_index]
  if (length(dim(S.sub)>0)) {
    met_index=which(apply(S.sub, 1, function(x) !all(x==0))==T)
  } else {
    met_index=which(S.sub!=0)
  }
  
  rxn_long = as.vector(unlist(model$rxnNames[rxn_index], use.names = FALSE))
  rxn_ss = as.vector(unlist(model$subSystems[rxn_index], use.names = FALSE))
  
  ############ rxn != genes!!!!!!!!!!!!!!!! #####################################################
  #   > length(model$rxns)
  #   [1] 7440
  #   > length(model$genes)
  #   [1] 2140
  genes_out = NULL
  for (i in 1:length(rxn_index)){
    gene_index = rxnGeneMat[rxn_index[i],]
    gene_index = which(gene_index!=0)
    gene_out = NULL
    for (j in 1:length(gene_index)){
      gene_out=paste(genes[gene_index[j]], gene_out, sep=";")  
    }
    genes_out=c(genes_out, gene_out)
  }
  
  met_long = as.vector(unlist(model$metNames[met_index], use.names = FALSE)) # 5062 ????????????????????
  met_short = as.vector(unlist(model$mets[met_index], use.names = FALSE))
  
  hmdb = model$metHMDB[met_index]
  tmp = rep(NA, length(hmdb))
  for (i in 1:length(hmdb)){
    if (identical(as.character(hmdb[[i]]), character(0))) tmp[i] = NA else tmp[i] = as.character(hmdb[[i]])   
  }
  hmdb=tmp

  kegg = model$metKeggID[met_index]
  tmp = rep(NA, length(kegg))
  for (i in 1:length(kegg)){
    if (identical(as.character(kegg[[i]]), character(0))) tmp[i] = NA else tmp[i] = as.character(kegg[[i]])   
  }
  kegg=tmp

  chebi = model$metCHEBIID[met_index]
  tmp = rep(NA, length(chebi))
  for (i in 1:length(chebi)){
    if (identical(as.character(chebi[[i]]), character(0))) tmp[i] = NA else tmp[i] = as.character(chebi[[i]])   
  }
  chebi=tmp

  pubchem = model$metPubChemID[met_index]
  tmp = rep(NA, length(pubchem))
  for (i in 1:length(pubchem)){
    if (identical(as.character(pubchem[[i]]), character(0))) tmp[i] = NA else tmp[i] = as.character(pubchem[[i]])   
  }
  pubchem=tmp

  gene_in = rep(gene_in, length(met_short))
  gene_in_ = rep(gene_in[1], length(genes_out))
  level = rep(level, length(met_short))
  level_ = rep(level[1], length(genes_out))
  met_in = rep(met, length(met_short))
  met_in_ = rep(met_in[1], length(genes_out))
  
  return(list("metabolites"=cbind(gene_in, level, met_in, met_short, met_long, hmdb, kegg, chebi, pubchem), "genes"=cbind(gene_in_, level_, met_in_, genes_out, rxn_long, rxn_ss), "matrix"=S[met_index,rxn_index], "rxns"=rxn_index, "mets"=met_index))
}
