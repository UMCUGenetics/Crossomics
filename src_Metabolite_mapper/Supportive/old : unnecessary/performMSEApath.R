performMSEApath <- function(metaboliteSet, p_valuesAll, patient, gene_in){
  # p_valuesAll=p.values.assi.all
  # patient=i  
  # gene_in=gene_list$HGNC[k]  
  
  #P.adj.thresh = 0.5
  #P.adj.thresh = 0.00005
  P.adj.thresh = 0.05
  
  p.values = p_valuesAll[order(p_valuesAll[,patient]),patient]
  
  # dit moet voor verdubbelen van rijen <====================================!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  # p.values.adj = p.adjust(p.values, method = "bonferroni") 
  
  all.mets.assi = as.integer(p.values < P.adj.thresh)
  names(all.mets.assi) = names(p.values)
  
  hmdb.min.adduct.rdndnt = as.vector(unlist(lapply(names(p.values),function(x) as.vector(unlist(strsplit(x, split =" ")))[1])))
  
  index = which(hmdb.min.adduct.rdndnt %in% metaboliteSet[,"hmdb"])
  
  metsInset = all.mets.assi[index]
  
  if (length(metsInset)>50) {
  ############################################################################################  
  # top 50 on p-value of set
  ints = p_valuesAll[names(metsInset), ]
  ints.ord = ints[order(ints[,paste("p.value_P",patient,sep="")], decreasing=FALSE),]
  ints.ord.top = ints.ord[1:50, -grep("p.value_", colnames(p_valuesAll), fixed=TRUE)] 
  metsInset = metsInset[rownames(ints.ord.top)]
  ############################################################################################
  }
  
  inSetAboveThresh=length(which(metsInset==1))
  inSetBelowThresh=length(which(metsInset==0))
  notInSetAboveThresh=length(which(all.mets.assi==1)) - inSetAboveThresh
  notInSetBelowThresh=length(which(all.mets.assi==0)) - inSetBelowThresh
  
  enrichment <-
    matrix(c(inSetAboveThresh, inSetBelowThresh, notInSetAboveThresh, notInSetBelowThresh),
           nrow = 2,
           dimnames = list(c("above_thresh", "below_thresh"),c("in_set", "not_in_setn")))
  #   print(enrichment)
  #   print(fisher.test(enrichment, alternative = "greater"))
  retVal = fisher.test(enrichment, alternative = "greater")
  
  # if (patient == 1) {
    
    ints = p_valuesAll[names(metsInset), -grep("p.value_", colnames(p_valuesAll), fixed=TRUE)]
    p_values = p_valuesAll[names(metsInset), grep("p.value_", colnames(p_valuesAll), fixed=TRUE)]
    
    names = rep("", dim(ints)[1])
    ids = as.vector(unlist(lapply(rownames(ints),function(x) as.vector(unlist(strsplit(x, split =" ")))[1])))
    for (j in 1:dim(ints)[1]){
      #message(pathway$metabolites[which(pathway$metabolites[,"hmdb"]==ids[j]),"name"])
      
      #       Warning message:
      #         In names[j] = metaboliteSet[which(metaboliteSet[, "hmdb"] == ids[j]),  :
      #                                       number of items to replace is not a multiple of replacement length
      
      # metaboliteSet[,"hmdb"] is not unique!!!!!!!!!!!!!!!!!!! 
      
      names[j] = metaboliteSet[which(metaboliteSet[,"hmdb"]==ids[j]), "name"]    # "met_long" <=====================!!!!!!!!!!!!!!!!!
    }
    
    adducts.label = as.vector(unlist(lapply(rownames(ints),function(x) as.vector(unlist(strsplit(x, split =" ")))[2])))
    adducts.label[which(adducts.label=="NA")]=""
    adducts.label[which(is.na(adducts.label))]=""
    names = paste(names, adducts.label)
    
    negpos.label = as.vector(unlist(lapply(rownames(ints),function(x) as.vector(unlist(strsplit(x, split =" ")))[3])))
    negpos.label[which(negpos.label=="NA")]=""
    negpos.label[which(is.na(negpos.label))]=""
    names = paste(names, negpos.label)
    
    rownames(ints) = names
    rownames(p_values) = names
    
    #CairoPNG(filename=paste("./results/crossomics/", gene_in, "/KEGG/P", patient,".png",sep=""), width, height) #, width, height
    CairoPNG(filename=paste("./results/crossomics2/P", patient, "/KEGG/", gene_in,".png",sep=""), width, height) #, width, height

    #hm = heatmap(as.matrix(ints), scale="row", Rowv = NA, distfun = spearman.dist, col=colorRampPalette(c("white","red"))(100), margins=c(6,6), labRow = rownames(ints), main=paste("Area around gene:", gene_in)) # Rowv = NA,  labCol = labcol, 
    hm = heatmap(as.matrix(ints), scale="row", distfun = spearman.dist,
                 col=colorRampPalette(c("white","red"))(100), margins=c(5,10),
                 labRow = rownames(ints), xlab=paste("Area around gene:", gene_in, ",\np-value: ",
                                                     round(retVal$p.value, digits = 3),", P",patient, sep = "")) # Rowv = NA,  labCol = labcol, 
    dev.off()
    
    ints = ints[hm$rowInd,hm$colInd]
    p_values = p_values[hm$rowInd,]
    ints = cbind("compound"=names[hm$rowInd], ints, p_values)  
    
    #genExcelFileShort(as.data.frame(ints), paste("./results/crossomics/", gene_in, "/KEGG/P", patient,".xls",sep=""))
    genExcelFileShort(as.data.frame(ints), paste("./results/crossomics2/P", patient, "/KEGG/", gene_in,".xls",sep=""))
    
  # }
  
  return(list("p.value"=retVal$p.value))
  
}
