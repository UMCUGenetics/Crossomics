performMSEAenv <- function(metaboliteSet, p_valuesAll, patient, gene_in, n_patients, thresh_F, top = 20){
# p_valuesAll=p.values.assi.all
# patient=patients[i]  

  width=1024
  height=768
  
  label=paste("p.value_P",patient,sep = "")
  p.values = p_valuesAll[order(abs(p_valuesAll[,label])),label]
  
  
  # p.values[grep(names(p.values), pattern = "HMDB01533",fixed = TRUE)]
  
  # dit moet voor verdubbelen van rijen <====================================!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  # p.values.adj = p.adjust(p.values, method = "bonferroni") 
  
  # p-values
  #all.mets.assi = as.integer(p.values < thresh_F)
  # z-scores
  all.mets.assi = as.integer(abs(p.values) > thresh_F)
  names(all.mets.assi) = names(p.values)
  
  # hmdb.min.adduct.rdndnt = 
    
  index = as.vector(unlist(lapply(names(p.values),function(x){
      any(as.vector(unlist(lapply(strsplit(as.vector(unlist(strsplit(as.vector(unlist(strsplit(x, split ="_")))[1], split =";")))[-1], split =" "), function(y) y[1]))) %in% metaboliteSet[,"hmdb"])
    })))

  metsInset = all.mets.assi[index]
  
  x=names(metsInset)[1]
  
  metNames = as.vector(unlist(lapply(names(metsInset), function(x){
    mode=as.vector(unlist(strsplit(x, split ="_")))[2]
    hmdb=unlist(lapply(strsplit(as.vector(unlist(strsplit(as.vector(unlist(strsplit(x, split ="_")))[1], split =";")))[-1], split =" "), function(y) y[1]))
    adduct=unlist(lapply(strsplit(as.vector(unlist(strsplit(as.vector(unlist(strsplit(x, split ="_")))[1], split =";")))[-1], split =" "), function(y) y[2]))
    index=which(hmdb %in% metaboliteSet[,"hmdb"] )
    hmdb=hmdb[index]
    adduct=adduct[index]
    paste(paste(paste(unlist(lapply(hmdb, function(y){unique(metaboliteSet[which(metaboliteSet[,"hmdb"]==y),"met_long"])})), adduct), collapse = ";"),mode,sep="_")})))
    
  metsInset=cbind(metsInset, "names"=metNames)
  
  if (dim(metsInset)[1]>top) {
############################################################################################  
  # top 50 on p-value of set
  ints = p_valuesAll[rownames(metsInset), ]
  # p.value
  # ints.ord = ints[order(ints[,paste("p.value_P",patient,sep="")], decreasing=FALSE),]
  # Z-score
  ints.ord = ints[order(abs(ints[,paste("p.value_P",patient,sep="")]), decreasing=TRUE),]
  ints.ord.top = ints.ord[1:top, -grep("p.value_", colnames(p_valuesAll), fixed=TRUE)] 
  metsInset = metsInset[rownames(ints.ord.top),]
############################################################################################
  }
  
#   index=which(names(CBSmetSet2) %in% names(PAHmetSet2))
#   metsInset = CBSmetSet2[-index]
 
#   # head tail 30 on intensity of set
#   ints = p_valuesAll[names(metsInset), -grep("p.value_", colnames(p_valuesAll), fixed=TRUE)]
#   ints.ord = ints[order(ints[,"P6"], decreasing=TRUE),]
#   metsInset = metsInset[c(rownames(ints.ord)[1:30],rownames(ints.ord)[(dim(ints.ord)[1]-30):dim(ints.ord)[1]])]

  inSetAboveThresh=length(which(metsInset==1))
  inSetBelowThresh=length(which(metsInset==0))
  notInSetAboveThresh=length(which(all.mets.assi==1)) - inSetAboveThresh
  notInSetBelowThresh=length(which(all.mets.assi==0)) - inSetBelowThresh
  
  message(paste(toString(gene_in), ", inSetAboveThresh: ", inSetAboveThresh, sep = ""))
  message(paste(toString(gene_in), ", inSetBelowThresh: ", inSetBelowThresh, sep = ""))
  message(paste(toString(gene_in), ", notInSetAboveThresh: ", notInSetAboveThresh, sep = ""))
  message(paste(toString(gene_in), ", notInSetBelowThresh: ", notInSetBelowThresh, sep = ""))
  
  enrichment <-
    matrix(c(inSetAboveThresh, inSetBelowThresh, notInSetAboveThresh, notInSetBelowThresh),
           nrow = 2,
           dimnames = list(c("above_thresh", "below_thresh"),c("in_set", "not_in_setn")))
  #   print(enrichment)
  #   print(fisher.test(enrichment, alternative = "greater"))
  retVal = fisher.test(enrichment, alternative = "greater")

  ints = p_valuesAll[rownames(metsInset), -grep("p.value_", colnames(p_valuesAll), fixed=TRUE),drop=FALSE]
  p_values = p_valuesAll[rownames(metsInset), grep("p.value_", colnames(p_valuesAll), fixed=TRUE),drop=FALSE]
  
#     names = rep("", dim(ints)[1])
#     ids = as.vector(unlist(lapply(rownames(ints),function(x) as.vector(unlist(strsplit(x, split =" ")))[1])))
#     for (j in 1:dim(ints)[1]){
#       names[j] = metaboliteSet[which(metaboliteSet[,"hmdb"]==ids[j])[1], "met_long"]    # "met_long" <=====================!!!!!!!!!!!!!!!!!
#     }
#     
#     adducts.label = as.vector(unlist(lapply(rownames(ints),function(x) as.vector(unlist(strsplit(x, split =" ")))[2])))
#     adducts.label[which(adducts.label=="NA")]=""
#     adducts.label[which(is.na(adducts.label))]=""
#     names = paste(names, adducts.label)
#     
#     negpos.label = as.vector(unlist(lapply(rownames(ints),function(x) as.vector(unlist(strsplit(x, split =" ")))[3])))
#     negpos.label[which(negpos.label=="NA")]=""
#     negpos.label[which(is.na(negpos.label))]=""
#     names = paste(names, negpos.label)
#     

  rownames(ints) = metsInset[,"names"]
  rownames(p_values) = metsInset[,"names"]  

  if (dim(metsInset)[1] > 1) {
    # CairoPNG(filename=paste("./results/crossomics/", gene_in, "/Recon2/P", patient,".png",sep=""), width, height) #, width, height
    CairoPNG(filename=paste("./results/crossomics/P", patient, "/Recon2/", gene_in,".png",sep=""), width, height) #, width, height
#    CairoPNG(filename=paste("./results/crossomics2/P1/Recon2/CBSnotinPAH.png",sep=""), width, height) #, width, height
    
    #hm = heatmap(as.matrix(ints), scale="row", Rowv = NA, distfun = spearman.dist, col=colorRampPalette(c("white","red"))(100), margins=c(6,6), labRow = rownames(ints), main=paste("Area around gene:", gene_in)) # Rowv = NA,  labCol = labcol, 
    hm = heatmap(as.matrix(ints), scale="row", distfun = spearman.dist,
                 col=colorRampPalette(c("white","red"))(100), margins=c(5,10),
                 labRow = rownames(ints), xlab=paste("Area around gene:", gene_in, ",\np-value: ",
                                                     round(retVal$p.value, digits = 3),", P",patient, sep = "")) # Rowv = NA,  labCol = labcol, 
    dev.off()

    ints = ints[hm$rowInd,hm$colInd]
    labels = colnames(p_values)
    p_values = as.data.frame(p_values[hm$rowInd,])
    colnames(p_values)=labels
  }
  
  ints = data.frame("compound"=rownames(ints), ints, "p.value"=as.numeric(p_values[,grep(colnames(p_values),pattern=toString(patient),fixed=TRUE)]))

  # genExcelFileShort(as.data.frame(ints), paste("./results/crossomics/", gene_in, "/Recon2/P", patient,".xls",sep=""))
  genExcelFileShort(as.data.frame(ints), paste("./results/crossomics/P", patient, "/Recon2/", gene_in,".xls",sep=""))
  
  
  return(list("p.value"=retVal$p.value))
#     
#   } else {
#     
#     return(NULL)
#     
#   }
  
}
