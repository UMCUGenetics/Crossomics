# claculate fold change
# two-tailed Student’s t-test with an alpha value of 0.01 to identify significant features

# probability of finding X >k significant genes in a particular gene set


# In een klas van 24 personen wordt door loting een groep van 4 personen samengesteld. Deze vier personen krijgen elk een andere taak.
# Op hoeveel manieren kan dit als deze vier personen pas na de loting hun taken onderling verdelen?
# Antwoord
# 
# Nu is de volgorde in de groep die wordt geloot niet van belang: ze verdelen pas na de loting onderling hun taken.
# 
# Het gaat nu dus om het aantal combinaties van 4 uit 24.
# 
# Er zijn daarom (24 over 4)=24!/4!⋅(24-4)!=24!/4!⋅20!=24⋅23⋅22⋅214!=23⋅22⋅21=10626 mogelijkheden

performMSEAenv <- function(metaboliteSet, p_valuesAll, patient, gene_in, n_patients, thresh_F_pos, thresh_F_neg, path, top = 20, id="hmdb", adductsSummed=FALSE){
  # p_valuesAll=p.values
  # patient=patients[i]
  
  width=1024
  height=768
  
  label=paste("p.value_P",patient,sep = "")
  p.values = p_valuesAll[,label]
  
  p.values = p.values[!is.na(p.values)]
  
  # z-scores
  all.mets.assi.pos = as.integer(p.values > thresh_F_pos)
  all.mets.assi.neg = as.integer(p.values < thresh_F_neg)
  all.mets.assi = all.mets.assi.pos + all.mets.assi.neg
  names(all.mets.assi) = names(p.values)
  
  # dit moet voor verdubbelen van rijen <====================================!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  # p.values.adj = p.adjust(p.values, method = "bonferroni") 
  
  #   # recon 2.2 ######################################################################  
  #   id = "InChI_key"  
  #   
  #   index1 = as.vector(unlist(lapply(names(p.values),function(x){
  #     # x=names(p.values)[3]
  #     tmp = as.vector(unlist(lapply(
  #         as.vector(unlist(strsplit(x, split =";")))[-1], function(y){
  #           # y=as.vector(unlist(strsplit(x, split =";")))[-1][1]
  #           as.vector(unlist(strsplit(y, split ="_")))[1]
  #         })))
  #     if (length(which(is.na(tmp)))>0) tmp = tmp[-which(is.na(tmp))]
  #     return(any(tmp %in% metaboliteSet[,id]))
  #   })))
  # 
  #   id = "chebi"  
  #   
  #   index2 = as.vector(unlist(lapply(names(p.values),function(x){
  #     # x=names(p.values)[3]
  #     tmp = as.vector(unlist(lapply(
  #       as.vector(unlist(strsplit(x, split =";")))[-1], function(y){
  #         # y=as.vector(unlist(strsplit(x, split =";")))[-1][2]
  #         as.vector(unlist(strsplit(y, split ="_")))[2]
  #       })))
  #     if (length(which(is.na(tmp)))>0) tmp = tmp[-which(is.na(tmp))]
  #     return(any(tmp %in% metaboliteSet[,id]))
  #   })))
  #   
  #   index = union(which(index1==TRUE),which(index2==TRUE))
  #   
  #   metsInset = all.mets.assi[index]
  #   metsInset = cbind(metsInset, "names"=names(metsInset))
  #   #####################################################################################
  
  
  
  # x=names(metsInset)[1]
  
  if (adductsSummed) {
    
    index = which(names(p.values ) %in% metaboliteSet[,"hmdb"])
    metsInset = all.mets.assi[index]
    
    metNames = rep("", length(metsInset))
    for (i in 1:length(metsInset)){
      metNames[i] = metaboliteSet[which(metaboliteSet[,id]==names(metsInset)[i])[1],"met_long"] 
    }
    
    if (length(metsInset)==0) return(NULL)
    
  } else{
    
    index = as.vector(unlist(lapply(names(p.values),function(x){
      #     any(as.vector(unlist(lapply(strsplit(as.vector(unlist(strsplit(as.vector(unlist(strsplit(x, split ="_")))[1], split =";")))[-1], split =" "), function(y) y[1]))) %in% metaboliteSet[,"hmdb"])
      any(as.vector(unlist(lapply(strsplit(as.vector(unlist(strsplit(as.vector(unlist(strsplit(x, split ="_")))[1], split =";"))), split =" "), function(y) y[1]))) %in% metaboliteSet[,"hmdb"])
    })))
    
    metsInset = all.mets.assi[index]
    
    # x=names(metsInset)[62]
    # grep("HMDB00062", names(metsInset), fixed =TRUE)
    
    metNames = as.vector(unlist(lapply(names(metsInset), function(x){
      mode=as.vector(unlist(strsplit(x, split ="_")))[2]
      #       hmdb=unlist(lapply(strsplit(as.vector(unlist(strsplit(as.vector(unlist(strsplit(x, split ="_")))[1], split =";")))[-1], split =" "), function(y) y[1]))
      hmdb=unlist(lapply(strsplit(as.vector(unlist(strsplit(as.vector(unlist(strsplit(x, split ="_")))[1], split =";"))), split =" "), function(y) y[1]))
      #       adduct=unlist(lapply(strsplit(as.vector(unlist(strsplit(as.vector(unlist(strsplit(x, split ="_")))[1], split =";")))[-1], split =" "), function(y) y[2]))
      adduct=unlist(lapply(strsplit(as.vector(unlist(strsplit(as.vector(unlist(strsplit(x, split ="_")))[1], split =";"))), split =" "), function(y) y[2]))
      index=which(hmdb %in% metaboliteSet[,id] )
      hmdb=hmdb[index]
      adduct=adduct[index]
      paste(paste(paste(unlist(lapply(hmdb, function(y){unique(metaboliteSet[which(metaboliteSet[,"hmdb"]==y),"met_long"])[1]})), adduct), collapse = ";"),mode,sep="_")})))
  }
  
  metsInset=cbind(metsInset, "names"= metNames)
  
  # Discard double rows
  metsInset = unique(metsInset, drop=FALSE)
  
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
  
  inSetAboveThresh=length(which(metsInset==1))
  inSetBelowThresh=length(which(metsInset==0))
  notInSetAboveThresh=length(which(all.mets.assi==1)) - inSetAboveThresh
  notInSetBelowThresh=length(which(all.mets.assi==0)) - inSetBelowThresh
  
  # Total number of mets
  N = inSetAboveThresh + inSetBelowThresh + notInSetAboveThresh + notInSetBelowThresh
  
  # Number of significant mets
  K = inSetAboveThresh + notInSetAboveThresh
  
  # Number of mets in set
  n = inSetAboveThresh + inSetBelowThresh
  
  # Number of significant mets in set
  k = inSetAboveThresh
  
  sum = 0
  for (i in 0:k){
    
    sum=sum +((choose(K,i) * choose((N-K),(n-i)))/choose(N,n)) 
    
  }
  p = 1-sum
  
  
  # message(paste(toString(gene_in), ", inSetAboveThresh: ", inSetAboveThresh, sep = ""))
  # message(paste(toString(gene_in), ", inSetBelowThresh: ", inSetBelowThresh, sep = ""))
  # message(paste(toString(gene_in), ", notInSetAboveThresh: ", notInSetAboveThresh, sep = ""))
  # message(paste(toString(gene_in), ", notInSetBelowThresh: ", notInSetBelowThresh, sep = ""))
  
  enrichment =
    matrix(c(inSetAboveThresh, inSetBelowThresh, notInSetAboveThresh, notInSetBelowThresh),
           nrow = 2,
           dimnames = list(c("above_thresh", "below_thresh"),c("in_set", "not_in_setn")))
  #   print(enrichment)
  #   print(fisher.test(enrichment, alternative = "greater"))
  retVal = fisher.test(enrichment, alternative = "greater")
  
  ints = p_valuesAll[rownames(metsInset), -grep("p.value_", colnames(p_valuesAll), fixed=TRUE),drop=FALSE]
  p_values = p_valuesAll[rownames(metsInset), grep("p.value_", colnames(p_valuesAll), fixed=TRUE),drop=FALSE]
  
  rownames(ints) = metsInset[,"names"]
  rownames(p_values) = metsInset[,"names"]
  
  #   # Discard double rows
  #   tmp = cbind(ints,p_values)
  #   tmp = unique(tmp, drop = FALSE)
  #
  #   p_values = tmp[,dim(tmp)[2], drop = FALSE]
  #   ints = tmp[,-dim(tmp)[2], drop = FALSE]
  
  # Remove negative and NA values
  ints[ints<0] = NA
  ints[is.na(ints)] = 0
  
  # Whole column zeros
  remove = which(apply(ints, 2, sum)==0)
  if (length(remove)>0) ints = ints[,-as.numeric(remove), drop=FALSE]
  
  if (length(p_values) > 1) {
    # CairoPNG(filename=paste("./results/crossomics/", gene_in, "/Recon2/P", patient,".png",sep=""), width, height) #, width, height
    CairoPNG(filename=paste(path, "/P", patient, "/Recon2/", gene_in,".png",sep=""), width, height) #, width, height
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
    metsInset = metsInset[hm$rowInd,]
  }
  
  ints = data.frame("compound"=rownames(ints), "HMDB"=rownames(metsInset), "Z score" = as.numeric(p_values[,grep(colnames(p_values),pattern=toString(patient),fixed=TRUE)]), ints)
  
  # genExcelFileShort(as.data.frame(ints), paste("./results/crossomics/", gene_in, "/Recon2/P", patient,".xls",sep=""))
  genExcelFileShort(as.data.frame(ints), paste(path, "/P", patient, "/Recon2/", gene_in,".xls",sep=""))
  
  
  return(list("p.value"=retVal$p.value))
  #     
  #   } else {
  #     
  #     return(NULL)
  #     
  #   }
  
}

