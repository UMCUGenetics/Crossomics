performMSEA <- function(metaboliteSet, p_valuesAll, patient, gene_in, n_patients, thresh_F_pos, thresh_F_neg, path, test, top = 20, id="hmdb", adductsSummed=FALSE){
# p_valuesAll=p.values
# patient=patients[i]
# test=3

  # require(gmp)
  width=1024
  height=768
  # a=0.5 # if 0 only crossing the thresh affects enrichment (x^0=1)
  # a=0
  # a=0.25
  # a=3
  a=1.2
  
  label=paste("p.value_P",patient,sep = "")
  p.values = p_valuesAll[,label]
  
  p.values = p.values[!is.na(p.values)]
  
  colnames(p_valuesAll)[grep("C", colnames(p_valuesAll), fixed = TRUE)]
  
  avg.int.c = apply(p_valuesAll[,grep("C", colnames(p_valuesAll), fixed = TRUE)], 1, mean)
  p_valuesAll = cbind(p_valuesAll,"avg.int.controls"=avg.int.c) 
  
  # z-scores  ????ME: make matrix for values significant (1) and not-significant (0) ??
  all.mets.assi.pos = as.integer(p.values > thresh_F_pos)
  all.mets.assi.neg = as.integer(p.values < thresh_F_neg)
  all.mets.assi = all.mets.assi.pos + all.mets.assi.neg
  #print(all.mets.assi)
  
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
    #print("jrrpa")
    #print(metsInset)
    # x=names(metsInset)[62]
    # grep("HMDB00062", names(metsInset), fixed =TRUE)s
    
    metNames = as.vector(unlist(lapply(names(metsInset), function(x){
      mode=as.vector(unlist(strsplit(x, split ="_")))[2]
#       hmdb=unlist(lapply(strsplit(as.vector(unlist(strsplit(as.vector(unlist(strsplit(x, split ="_")))[1], split =";")))[-1], split =" "), function(y) y[1]))
      hmdb=unlist(lapply(strsplit(as.vector(unlist(strsplit(as.vector(unlist(strsplit(x, split ="_")))[1], split =";"))), split =" "), function(y) y[1]))
#       adduct=unlist(lapply(strsplit(as.vector(unlist(strsplit(as.vector(unlist(strsplit(x, split ="_")))[1], split =";")))[-1], split =" "), function(y) y[2]))
      adduct=unlist(lapply(strsplit(as.vector(unlist(strsplit(as.vector(unlist(strsplit(x, split ="_")))[1], split =";"))), split =" "), function(y) y[2]))
      index=which(hmdb %in% metaboliteSet[,id] )
      hmdb=unique(hmdb[index])
      adduct=unique(adduct[index])
      paste(paste(paste(unlist(lapply(hmdb, function(y){unique(metaboliteSet[which(metaboliteSet[,"hmdb"]==y),"met_long"])[1]})), adduct), collapse = "; "),mode,sep="_")
      })))
    paths = as.vector(unlist(lapply(names(metsInset), function(x){
      # x=names(metsInset)[6]
      hmdb=unlist(lapply(strsplit(as.vector(unlist(strsplit(as.vector(unlist(strsplit(x, split ="_")))[1], split =";"))), split =" "), function(y) y[1]))
      index=which(hmdb %in% metaboliteSet[,id] )
      hmdb=unique(hmdb[index])
      paste(unlist(lapply(hmdb, function(y){unique(metaboliteSet[which(metaboliteSet[,"hmdb"]==y),"path"])[1]})), collapse = "; ")
      })))
    hmdb_set = as.vector(unlist(lapply(names(metsInset), function(x){
      # x=names(metsInset)[3]
      hmdb_pkgrp=unlist(lapply(strsplit(as.vector(unlist(strsplit(as.vector(unlist(strsplit(x, split ="_")))[1], split =";"))), split =" "), function(y) y[1]))
      hmdb_pkgrp_adduct=as.vector(unlist(strsplit(as.vector(unlist(strsplit(x, split ="_")))[1], split =";")))
      index=which(hmdb_pkgrp %in% metaboliteSet[,id] )
      hmdb_metSet=paste(unique(hmdb_pkgrp_adduct[index]), collapse = ", ")
    })))
  #print(metNames)  
  }
  foldChange = p_valuesAll[,2]/avg.int.c
  g_i = round(abs(log2(foldChange))^a)
  # g_i = round(abs(p.values)^a)
  metsInset=cbind(metsInset, "names"= metNames, "hmdb_set"= hmdb_set, "g_i"=g_i[index], "fc"=foldChange[index], "z-score"=p.values[names(metsInset)], "path"=paths)
  # Discard double rows
  metsInset = unique(metsInset, drop=FALSE)
  # tmp = metsInset
  # rownames(tmp) = NULL
  # print(tmp[,c("g_i","fc","z-score")])

#   if (dim(metsInset)[1]>top) {
# ############################################################################################
#   # top 50 on p-value of set
#   ints = p_valuesAll[rownames(metsInset), ]
#   # p.value
#   # ints.ord = ints[order(ints[,paste("p.value_P",patient,sep="")], decreasing=FALSE),]
#   # Z-score
#   ints.ord = ints[order(abs(ints[,paste("p.value_P",patient,sep="")]), decreasing=TRUE),]
#   ints.ord.top = ints.ord[1:top, -grep("p.value_", colnames(p_valuesAll), fixed=TRUE)]
#   metsInset = metsInset[rownames(ints.ord.top),]
# ############################################################################################
#   }
  inSetAboveThresh=length(which(metsInset[,"metsInset"]==1))
  inSetBelowThresh=length(which(metsInset[,"metsInset"]==0))
  notInSetAboveThresh=length(which(all.mets.assi==1)) - inSetAboveThresh
  notInSetBelowThresh=length(which(all.mets.assi==0)) - inSetBelowThresh
  # message(paste(toString(gene_in), ", inSetAboveThresh: ", inSetAboveThresh, sep = ""))
  # message(paste(toString(gene_in), ", inSetBelowThresh: ", inSetBelowThresh, sep = ""))
  # message(paste(toString(gene_in), ", notInSetAboveThresh: ", notInSetAboveThresh, sep = ""))
  # message(paste(toString(gene_in), ", notInSetBelowThresh: ", notInSetBelowThresh, sep = ""))
  
  if (test==0){
    # Fishers exact test ######################################################################
    ###########################################################################################
    enrichment =
      matrix(c(inSetAboveThresh, inSetBelowThresh, notInSetAboveThresh, notInSetBelowThresh),
             nrow = 2,
             dimnames = list(c("above_thresh", "below_thresh"),c("in_set", "not_in_setn")))
    #   print(enrichment)
    #   print(fisher.test(enrichment, alternative = "greater"))
    retVal = fisher.test(enrichment, alternative = "greater")
    p = retVal$p.value
    ###########################################################################################
  
  } else if (test==1){
    # Hypergeometric test #####################################################################
    ###########################################################################################
    # Total number of mets
    N = inSetAboveThresh + inSetBelowThresh + notInSetAboveThresh + notInSetBelowThresh

    # Number of significant mets
    m = inSetAboveThresh + notInSetAboveThresh

    # Number of mets in set
    n = inSetAboveThresh + inSetBelowThresh

    # Number of significant mets in set
    k = inSetAboveThresh

    if (n==0){
      p=1
    } else {
      
      sum = 0
      for (i in 0:k){
  
        sum=sum +((choose(m,i) * choose((N-m),(n-i)))/choose(N,n))
  
      }
      p = 1-sum
    }
    ###########################################################################################
    
  } else if (test==2){
    # Weighted hypergeometric test ############################################################
    ###########################################################################################
    # Number of mets in set
    n = inSetAboveThresh + inSetBelowThresh

    if (n==0 | dim(metsInset)[1]==0){
      p=1
    } else {

      Q = max(as.numeric(metsInset[,"g_i"]))
    
      # Pseudoset
      pseudoSet = NULL
      for (i in 1:dim(metsInset)[1]){
        name = rownames(metsInset)[i]
        for (j in 1:as.numeric(metsInset[i,"g_i"])){
          pseudoSet = rbind(pseudoSet, metsInset[i,])
          rownames(pseudoSet)[dim(pseudoSet)[1]] = paste(name, j, sep="_")
        }
      }
    
      k_pseudo = length(which(pseudoSet[,"metsInset"]==1))
      # k_pseudo = dim(pseudoSet)[1]
      n_pseudo = dim(pseudoSet)[1]
      
      # Total number of mets
      N = inSetAboveThresh + inSetBelowThresh + notInSetAboveThresh + notInSetBelowThresh
  
      sigma = 0
      for (i in 0:k_pseudo){
    
        sigma=sigma +((chooseZ(N,i) * chooseZ((Q*N-N),(Q*n_pseudo-i)))/chooseZ(Q*N,Q*n_pseudo))
    
      }
      p = 1-as.double(sigma)
    }
    ###########################################################################################
  } else if (test==3){
    # Weighted Fishers exact test ############################################################
    ###########################################################################################
    
    if (dim(metsInset)[1]==0){
      p=1
    } else {
      # Pseudoset
      pseudoSet = NULL
      for (i in 1:dim(metsInset)[1]){
        name = rownames(metsInset)[i]

        
        #print(metsInset[i,"g_i"])  ## ME
        
        for (j in 1:as.numeric(metsInset[i,"g_i"])){
          pseudoSet = rbind(pseudoSet, metsInset[i,])
          rownames(pseudoSet)[dim(pseudoSet)[1]] = paste(name, j, sep="_")
        }
      }
      index=as.vector(which(pseudoSet[,"g_i"]==0))
      if (length(index)!=0) pseudoSet = pseudoSet[-index,]

      if (is.null(dim(pseudoSet))){
        f = 1/top
        
        if (pseudoSet["metsInset"]==1) {
          inSetAboveThresh=1
          inSetBelowThresh=0
        } else {
          inSetAboveThresh=0
          inSetBelowThresh=1
        }
      } else {
        f = dim(pseudoSet)[1]/top
        
        inSetAboveThresh=length(which(pseudoSet[,"metsInset"]==1))
        inSetBelowThresh=length(which(pseudoSet[,"metsInset"]==0))
      }

      notInSetAboveThresh=round(f*notInSetAboveThresh)
      notInSetBelowThresh=round(f*notInSetBelowThresh)
      
      enrichment =
        matrix(c(inSetAboveThresh, inSetBelowThresh, notInSetAboveThresh, notInSetBelowThresh),
               nrow = 2,
               dimnames = list(c("above_thresh", "below_thresh"),c("in_set", "not_in_setn")))
      #   print(enrichment)
      #   print(fisher.test(enrichment, alternative = "greater"))
      retVal = fisher.test(enrichment, alternative = "greater")
      p = retVal$p.value
      
      genExcelFileShort(as.data.frame(pseudoSet), paste(path, "/P", patient, "/Recon2/", gene_in,"_pseudo.xls",sep=""))
    }  
  }
  tmp1=NULL
  index=which(as.numeric(metsInset[,"z-score"])<thresh_F_neg)
  if (length(index)>0) tmp1=metsInset[index,,drop=FALSE]
  tmp2=NULL
  index=which(as.numeric(metsInset[,"z-score"])>thresh_F_pos)
  if (length(index)>0) tmp2=metsInset[index,,drop=FALSE]
  metsInset=rbind(tmp1,tmp2)
  
  ints = p_valuesAll[rownames(metsInset), -grep("p.value_", colnames(p_valuesAll), fixed=TRUE),drop=FALSE]
  p_values = p_valuesAll[rownames(metsInset), grep("p.value_", colnames(p_valuesAll), fixed=TRUE),drop=FALSE]

  rownames(ints) = metsInset[,"names"]
  rownames(p_values) = metsInset[,"names"]
  
  ints = ints[, -grep("avg.int.controls", colnames(ints), fixed=TRUE),drop=FALSE]

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

      # colnames(ints)[1] = "P22"
    
      hm = heatmap(as.matrix(ints),
                   scale="row",
                   #distfun = "euclidean", #spearman.dist,
                   col=colorRampPalette(c("yellow","blue"))(100),
                   margins=c(2.5,26),
                   labRow = metsInset[,"hmdb_set"])
    
    # Rowv = NA,  labCol = labcol, , xlab=paste("Enriched in metabolic environment of:", gene_in, ",\np-value: ",  round(p, digits = 8), sep = "")
    dev.off()

    ints = ints[hm$rowInd,hm$colInd]
    labels = colnames(p_values)
    p_values = as.data.frame(p_values[hm$rowInd,])
    colnames(p_values)=labels
    metsInset = metsInset[hm$rowInd,]
  }

  ints = data.frame("compound"=rownames(ints),
                    "HMDB_set"=metsInset[,"hmdb_set"], 
                    "Z score" = as.numeric(p_values[,grep(colnames(p_values),pattern=toString(patient),fixed=TRUE)]),
                    "path"=metsInset[,"path"],
                    "HMDB"=rownames(metsInset),
                    ints)

  ints=ints[order(ints[,"Z.score"]),]
  
  # genExcelFileShort(as.data.frame(ints), paste("./results/crossomics/", gene_in, "/Recon2/P", patient,".xls",sep=""))
  genExcelFileShort(as.data.frame(ints), paste(path, "/P", patient, "/Recon2/", gene_in,".xls",sep=""))
  
  
  return(list("p.value"=p))
#     
#   } else {
#     
#     return(NULL)
#     
#   }
  
}
