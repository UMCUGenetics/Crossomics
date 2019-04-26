performMSEA <- function(metaboliteSet, av_int_and_z_values_matrix, patient, gene_in, n_patients, thresh_F_pos, thresh_F_neg, path, test, top = 20, id="hmdb", patient_folder){
  # av_int_and_z_values_matrix=patient_z_values_vector
  # patient=patients[i]
  # test=3
  
  # require(gmp)
  # set the dimensions of the pictures for later
  width=1024
  height=768
  # Fold Change is put to the power logFC_weight (weighted_logFC = round(abs(log2(foldChange))^logFC_weight))
  # logFC_weight=0.5 # if 0 only crossing the thresh affects enrichment (x^0=1)
  # logFC_weight=0
  # logFC_weight=0.25
  # logFC_weight=3
  logFC_weight=1.2
  
  # label=paste("z.average_P",patient,sep = "")
  label <- paste0("av.z_", patient)
  # label=paste("p.value_P",patient,sep = "")
  patient_z_values_vector = av_int_and_z_values_matrix[,label]
  
  # This is logFC_weight check, I think there shouldn't be any is.na(p)
  patient_z_values_vector = patient_z_values_vector[!is.na(patient_z_values_vector)]
  
  # average of the control intensity values and paste it to the end of the previous matrix
  avg.int.c <- apply(av_int_and_z_values_matrix[,grep("C", colnames(av_int_and_z_values_matrix), fixed = TRUE)], 1, mean)
  av_int_and_z_values_matrix <- cbind(av_int_and_z_values_matrix,"avg.int.controls"=avg.int.c) 
  
  # Z-scores that exceed the positive/negative threshold
  # allMetsExcedingThres.pos = as.integer(patient_z_values_vector > thresh_F_pos)
  # allMetsExcedingThres.neg = as.integer(patient_z_values_vector < thresh_F_neg)
  # allMetsExcedingThres = allMetsExcedingThres.pos + allMetsExcedingThres.neg
  # 
  # names(allMetsExcedingThres) = names(patient_z_values_vector)
  
  # This one replaces the 4 lines above: Get boolean of all HMDB codes indicating exceding either one of the threshold values
  allMetsExcedingThres <- patient_z_values_vector > thresh_F_pos | patient_z_values_vector < thresh_F_neg
  
  # dit moet voor verdubbelen van rijen <====================================!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  # patient_z_values_vector.adj = p.adjust(patient_z_values_vector, method = "bonferroni") 
  
  #   # recon 2.2 ######################################################################  
  #   id = "InChI_key"  
  #   
  #   index1 = as.vector(unlist(lapply(names(patient_z_values_vector),function(x){
  #     # x=names(patient_z_values_vector)[3]
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
  #   index2 = as.vector(unlist(lapply(names(patient_z_values_vector),function(x){
  #     # x=names(patient_z_values_vector)[3]
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
  #   metsInset = allMetsExcedingThres[index]
  #   metsInset = cbind(metsInset, "names"=names(metsInset))
  #   #####################################################################################
  
  # 1. split the rownames of patient_z_values_vector on " ", ";" and "_" and keep only the HMDB identifiers
  # 2. get logFC_weight boolean vector of length nrow(patient_z_values_vector) where the identifiers are present in the metaboliteSet[,"hmdb"]
  # index <- unlist(lapply(names(patient_z_values_vector),function(x){
  #   #     any(as.vector(unlist(lapply(strsplit(as.vector(unlist(strsplit(as.vector(unlist(strsplit(x, split ="_")))[1], split =";")))[-1], split =" "), function(y) y[1]))) %in% metaboliteSet[,"hmdb"])
  #   any(unlist(lapply(strsplit(unlist(strsplit(unlist(strsplit(x, split ="_"))[1], split =";")), split =" "), function(y) y[1])) %in% metaboliteSet[,"hmdb"])
  # }))
  # This is the replacement of the above function, get the location of all patient_z_values_vector-names present in the metaboliteSet
  index <- names(patient_z_values_vector) %in% metaboliteSet[,"hmdb"]
  # x=names(metsInset)[1]
  metsInset <- data.frame("InSetAboveThres" = allMetsExcedingThres[index])
  # metsInset <- allMetsExcedingThres
  if(nrow(metsInset) == 0)   return(list("p.value"= 1))
  
  
  # if (adductsSummed) {
  #   # NEW ONE, AS IN THE adductsSummed == FALSE, BELOW, ALSO ABOVE THIS IF-ELSE STATEMENT
  #   # index = as.vector(unlist(lapply(names(patient_z_values_vector),function(x){
  #   #   #     any(as.vector(unlist(lapply(strsplit(as.vector(unlist(strsplit(as.vector(unlist(strsplit(x, split ="_")))[1], split =";")))[-1], split =" "), function(y) y[1]))) %in% metaboliteSet[,"hmdb"])
  #   #   any(as.vector(unlist(lapply(strsplit(as.vector(unlist(strsplit(as.vector(unlist(strsplit(x, split ="_")))[1], split =";"))), split =" "), function(y) y[1]))) %in% metaboliteSet[,"hmdb"])
  #   # })))
  #   # OLD ONE, NOT WORKING
  #   # index = which(names(patient_z_values_vector ) %in% metaboliteSet[,"hmdb"])
  #   # metsInset = allMetsExcedingThres[index]
  #   
  #   metNames = rep("", length(metsInset))
  #   for (i in 1:length(metsInset)){
  #     # metNames[i] = metaboliteSet[which(metaboliteSet[,id]==names(metsInset)[i])[1],"met_long"]
  #     lapply(metaboliteSet[,id], function(x) grep(x, names(metsInset)))
  #     
  #     
  #   }
  #   
  #   if (length(metsInset)==0) return(NULL)
  #   
  # } else{
  #   # index = as.vector(unlist(lapply(names(patient_z_values_vector),function(x){
  #   #   #     any(as.vector(unlist(lapply(strsplit(as.vector(unlist(strsplit(as.vector(unlist(strsplit(x, split ="_")))[1], split =";")))[-1], split =" "), function(y) y[1]))) %in% metaboliteSet[,"hmdb"])
  #   #   any(as.vector(unlist(lapply(strsplit(as.vector(unlist(strsplit(as.vector(unlist(strsplit(x, split ="_")))[1], split =";"))), split =" "), function(y) y[1]))) %in% metaboliteSet[,"hmdb"])
  #   # })))
  #   # metsInset = allMetsExcedingThres[index]
  #   #print("jrrpa")
  #   #print(metsInset)
  #   # x=names(metsInset)[62]
  #   # grep("HMDB00062", names(metsInset), fixed =TRUE)s
  #   

  # Get long names of the hmdb's in the set.
  # metNames <- as.vector(unlist(lapply(names(metsInset), function(x) metaboliteSet[metaboliteSet[, "hmdb"] == x,"met_long"])))
  
  # Because of change of class = matrix to class = dataframe of metsInset:
  # metNames <- as.vector(unlist(lapply(rownames(metsInset), function(x) metaboliteSet[metaboliteSet[, "hmdb"] == x,"met_long"])))
  
  # Automatically takes the first hmdb hit, so duplicates will be removed (first occurence will be kept)
  rownames(metaboliteSet) <- metaboliteSet[,"hmdb"]
  metNames <- metaboliteSet[rownames(metsInset),"met_long"]

  # Previously: get long names + adduct charge (_pos/_neg) + the adduct code ([M+Na]+ or something similar)
    # metNames = unlist(lapply(names(metsInset), function(x){
    #   # mode=as.vector(unlist(strsplit(x, split ="_")))[2] #old, but as.vector is unnecessary
    #   # mode: creates pos or neg
    #   mode=unlist(strsplit(x, split ="_"))[2]
    #   #       hmdb=unlist(lapply(strsplit(as.vector(unlist(strsplit(as.vector(unlist(strsplit(x, split ="_")))[1], split =";")))[-1], split =" "), function(y) y[1]))
    #   # hmdb=unlist(lapply(strsplit(as.vector(unlist(strsplit(as.vector(unlist(strsplit(x, split ="_")))[1], split =";"))), split =" "), function(y) y[1])) #old, but as.vector seems unnecessary
    #   # hmdb: created the hmdb codes
    #   hmdb=unlist(lapply(strsplit(unlist(strsplit(unlist(strsplit(x, split ="_"))[1], split =";")), split =" "), function(y) y[1]))
    #   #       adduct=unlist(lapply(strsplit(as.vector(unlist(strsplit(as.vector(unlist(strsplit(x, split ="_")))[1], split =";")))[-1], split =" "), function(y) y[2]))
    #   # adduct=unlist(lapply(strsplit(as.vector(unlist(strsplit(as.vector(unlist(strsplit(x, split ="_")))[1], split =";"))), split =" "), function(y) y[2])) #old, but as.vector seems unnecessary
    #   # adduct: creates the "[M+Na]+" "[M+K]+" adducts
    #   adduct=unlist(lapply(strsplit(unlist(strsplit(unlist(strsplit(x, split ="_"))[1], split =";")), split =" "), function(y) y[2]))
    #   index=which(hmdb %in% metaboliteSet[,id] )
    #   hmdb=unique(hmdb[index])
    #   adduct=unique(adduct[index])
    #   print(paste("mode:", mode, " hmdb:", hmdb, " adduct:", adduct))
    #   # paste the unique, long name of the metabolites in the hmdb column (of metaboliteSet) to the adduct, collapse it with logFC_weight ; and paste with the mode
    #   paste(paste(paste(unlist(lapply(hmdb, function(y){unique(metaboliteSet[which(metaboliteSet[,"hmdb"]==y),"met_long"])[1]})), adduct), collapse = "; "),mode,sep="_")
    # })) 
  
  # All reaction paths present for the hmdb codes
    # paths = unlist(lapply(names(metsInset), function(x){
    paths <- unlist(lapply(rownames(metsInset), function(x){
      # x=names(metsInset)[6]
      # hmdb=unlist(lapply(strsplit(unlist(strsplit(unlist(strsplit(x, split ="_"))[1], split =";")), split =" "), function(y) y[1]))
      hmdb <- x
      index=which(hmdb %in% metaboliteSet[,id] )
      hmdb=unique(hmdb[index])
      paste(unlist(lapply(hmdb, function(y){unique(metaboliteSet[which(metaboliteSet[,"hmdb"]==y),"path"])[1]})), collapse = "; ")
    }))
    # All hmdbs present for this 
    # hmdb_set = unlist(lapply(names(metsInset), function(x){
    hmdb_set = unlist(lapply(rownames(metsInset), function(x){
      # x=names(metsInset)[3]
      hmdb_pkgrp=unlist(lapply(strsplit(unlist(strsplit(unlist(strsplit(x, split ="_"))[1], split =";")), split =" "), function(y) y[1]))
      hmdb_pkgrp_adduct=unlist(strsplit(unlist(strsplit(x, split ="_"))[1], split =";"))
      index=which(hmdb_pkgrp %in% metaboliteSet[,id] )
      hmdb_metSet=paste(unique(hmdb_pkgrp_adduct[index]), collapse = ", ")
    }))
    #print(metNames)  
  # }
    
  # Used for weighted Fishers test
  # foldChange = av_int_and_z_values_matrix[,2]/avg.int.c # Dangerous, logFC_weight value of 1 means no change
  # weighted_logFC = round(abs(log2(foldChange))^logFC_weight)
  # weighted_logFC = round(abs(patient_z_values_vector)^logFC_weight)
  # metsInset=cbind(metsInset, "names"= metNames, "hmdb_set"= hmdb_set, "weighted_logFC"=weighted_logFC[index], "FC"=foldChange[index], "z-score"=patient_z_values_vector[names(metsInset)], "path"=paths)
  metsInset=cbind(metsInset, 
                  "names"= metNames, 
                  "hmdb_set"= hmdb_set, 
                  # "weighted_logFC"=weighted_logFC[index], 
                  # "FC"=foldChange[index], 
                  "z-score"=patient_z_values_vector[rownames(metsInset)], 
                  "path"=paths)
  
  # Discard double rows
  metsInset = unique(metsInset, drop=FALSE)
  # tmp = metsInset
  # rownames(tmp) = NULL
  # print(tmp[,c("weighted_logFC","FC","z-score")])
  
  #   if (nrow(metsInset)>top) {
  # ############################################################################################
  #   # top 50 on p-value of set
  #   ints = av_int_and_z_values_matrix[rownames(metsInset), ]
  #   # p.value
  #   # ints.ord = ints[order(ints[,paste("p.value_P",patient,sep="")], decreasing=FALSE),]
  #   # Z-score
  #   ints.ord = ints[order(abs(ints[,paste("p.value_P",patient,sep="")]), decreasing=TRUE),]
  #   ints.ord.top = ints.ord[1:top, -grep("p.value_", colnames(av_int_and_z_values_matrix), fixed=TRUE)]
  #   metsInset = metsInset[rownames(ints.ord.top),]
  # ############################################################################################
  #   }
  # inSetAboveThresh=length(which(metsInset[,"metsInset"]==1))
  # inSetBelowThresh=length(which(metsInset[,"metsInset"]==0))
  
  # Replaces the two lines above, because allMetsExcedingThres is not 0 / 1 anymore, but FALSE / TRUE, and it could be shorter
  # inSetAboveThresh=length(which(metsInset[,"metsInset"]==TRUE))
  # inSetBelowThresh=length(which(metsInset[,"metsInset"]==FALSE))
  inSetAboveThresh=sum(metsInset[,"InSetAboveThres"])
  inSetBelowThresh=sum(!metsInset[,"InSetAboveThres"])
  
  # notInSetAboveThresh=length(which(allMetsExcedingThres==1)) - inSetAboveThresh
  # notInSetBelowThresh=length(which(allMetsExcedingThres==0)) - inSetBelowThresh
  
  # Replaces the two lines above, because allMetsExcedingThres is not 0 / 1 anymore, but FALSE / TRUE, and it could be shorter
  notInSetAboveThresh <- sum(allMetsExcedingThres) - inSetAboveThresh
  notInSetBelowThresh <- sum(!allMetsExcedingThres) - inSetBelowThresh
  
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
             dimnames = list(c("above_thresh", "below_thresh"),c("in_set", "not_in_set")))
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
    
    if (n==0 | nrow(metsInset)==0){
      p=1
    } else {
      
      Q = max(as.numeric(metsInset[,"weighted_logFC"]))
      
      # Pseudoset
      pseudoSet = NULL
      for (i in 1:nrow(metsInset)){
        name = rownames(metsInset)[i]
        for (j in 1:as.numeric(metsInset[i,"weighted_logFC"])){
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
    
    # If no metabolites left (rows of metsInset), p value = 1
    if (nrow(metsInset)==0){
      p=1
    } else {
      # Pseudoset
      pseudoSet <- NULL
      for (i in 1:nrow(metsInset)){
        name <- rownames(metsInset)[i]
        
        # To inflate the value of metabolites according to its intensity (small intensity elevation vs large elevation), 
        # Metabolites are duplicated to be included multiple times in the fishers test on the basis of their fold change
        # Possible problem: if the weighted logFC = 0, genes will also be present twice in the pseudoset.
        # This method isn't used much... references date from 2015 and 2016
        for (j in 1:as.numeric(metsInset[i,"weighted_logFC"])){
          if (j == 0) next
          pseudoSet = rbind(pseudoSet, metsInset[i,])
          rownames(pseudoSet)[nrow(pseudoSet)] = paste(name, j, sep="_")
        }
      }
      # index=as.vector(which(pseudoSet[,"weighted_logFC"]==0))
      # if (length(index)!=0) pseudoSet <- pseudoSet[-index,]
      
      # if (is.null(nrow(pseudoSet))){
        if (sum(pseudoSet$InSetAboveThres == 0)){
          f = 1/top
        
        # f = 1/top

        # This metsInset column or whatever was a vector of 1 and 0 for whether the metabolite was in the set
        # if (pseudoSet["metsInset"]==1) {
        #   inSetAboveThresh=1
        #   inSetBelowThresh=0
        # } else {
        #   inSetAboveThresh=0
        #   inSetBelowThresh=1
        # }
      } else {
        f <- sum(pseudoSet$InSetAboveThres) / top
        # f = nrow(pseudoSet)/top
# 
#         inSetAboveThresh=length(which(pseudoSet[,"metsInset"]==1))
#         inSetBelowThresh=length(which(pseudoSet[,"metsInset"]==0))
      }
      
  inSetAboveThresh <- sum(pseudoSet$InSetAboveThres)
  inSetBelowThresh <- sum(!pseudoSet$InSetAboveThres)
      # If there are few genes in the set, notInSetAboveThresh will be made smaller (f = genesInSet / top (=20))
      # Thus: few genes in set --> inflate importance for those genes
      # The effect of this method seems to be very small though. 
      # On what basis would you determine what top should be?
      notInSetAboveThresh=round(f*notInSetAboveThresh)
      notInSetBelowThresh=round(f*notInSetBelowThresh)
      
      enrichment =
        matrix(c(inSetAboveThresh, inSetBelowThresh, notInSetAboveThresh, notInSetBelowThresh),
               nrow = 2,
               dimnames = list(c("above_thresh", "below_thresh"),c("in_set", "not_in_set")))
      #   print(enrichment)
      #   print(fisher.test(enrichment, alternative = "greater"))
      retVal = fisher.test(enrichment, alternative = "greater")
      p = retVal$p.value
      
      # genExcelFileShort(as.data.frame(pseudoSet), paste(path, "/P", patient, "/Recon2/", gene_in,"_pseudo.xls",sep=""))
      # genExcelFileShort(as.data.frame(pseudoSet), paste(path, "/", patient_folder, "/Recon2/", gene_in,"_pseudo.xls",sep=""))
      genExcelFileShort(as.data.frame(pseudoSet), paste(path, "/", patient_folder, "/", gene_in,"_pseudo.xls",sep=""))
    }  
  }
  
  # # Code for creating just the metsInset exceding the thresholds
  # tmp1=NULL
  # # index=which(as.numeric(metsInset[,"z-score"])<thresh_F_neg)
  # index=which(metsInset[,"z-score"]<thresh_F_neg) 
  # if (length(index)>0) tmp1=metsInset[index,,drop=FALSE]
  # tmp2=NULL
  # # index=which(as.numeric(metsInset[,"z-score"])>thresh_F_pos)
  # index=which(metsInset[,"z-score"]>thresh_F_pos)
  # if (length(index)>0) tmp2=metsInset[index,,drop=FALSE]
  # metsInset=rbind(tmp1,tmp2)
  
  # This line replaces the chuck above this
  metsInset <- metsInset[metsInset$InSetAboveThres,]
  
  # Changed name of p.value_ to z.average_
  # ints = av_int_and_z_values_matrix[rownames(metsInset), -grep("p.value_", colnames(av_int_and_z_values_matrix), fixed=TRUE),drop=FALSE]
  # z_values = av_int_and_z_values_matrix[rownames(metsInset), grep("p.value_", colnames(av_int_and_z_values_matrix), fixed=TRUE),drop=FALSE]
  # ints <- av_int_and_z_values_matrix[rownames(metsInset), -grep("z.average_", colnames(av_int_and_z_values_matrix), fixed=TRUE),drop=FALSE]
  # z_values <- av_int_and_z_values_matrix[rownames(metsInset), grep("z.average_", colnames(av_int_and_z_values_matrix), fixed=TRUE),drop=FALSE]
  ints <- av_int_and_z_values_matrix[rownames(metsInset), -grep("av.z", colnames(av_int_and_z_values_matrix), fixed=TRUE),drop=FALSE]
  z_values <- av_int_and_z_values_matrix[rownames(metsInset), grep("av.z", colnames(av_int_and_z_values_matrix), fixed=TRUE),drop=FALSE]
  
  
  rownames(ints) = metsInset[,"names"]
  rownames(z_values) = metsInset[,"names"]
  
  ints <- ints[, -grep("avg.int.controls", colnames(ints), fixed=TRUE),drop=FALSE]
  
  #   # Discard double rows
  #   tmp = cbind(ints,z_values)
  #   tmp = unique(tmp, drop = FALSE)
  #
  #   z_values = tmp[,dim(tmp)[2], drop = FALSE]
  #   ints = tmp[,-dim(tmp)[2], drop = FALSE]
  
  # Remove negative and NA values
  ints[ints<0] <- NA
  ints[is.na(ints)] <- 0
  
  # Whole column zeros
  remove <- which(colSums(ints)==0)
  if (length(remove)>0) ints = ints[,-as.numeric(remove), drop=FALSE]
  
  if (length(z_values) > 1) {
    # CairoPNG(filename=paste("./results/crossomics/", gene_in, "/Recon2/P", patient,".png",sep=""), width, height) #, width, height
    # CairoPNG(filename=paste(path, "/P", patient, "/Recon2/", gene_in,".png",sep=""), width, height) #, width, height
    # CairoPNG(filename=paste(path, "/", patient_folder, "/Recon2/", gene_in,".png",sep=""), width, height)
    CairoPNG(filename=paste(path, "/", patient_folder, "/", gene_in,".png",sep=""), width, height)
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
    labels = colnames(z_values)
    z_values = as.data.frame(z_values[hm$rowInd,])
    colnames(z_values)=labels
    metsInset = metsInset[hm$rowInd,]
  }
  
  ints = data.frame("compound"=rownames(ints),
                    "HMDB_set"=metsInset[,"hmdb_set"], 
                    "Z score" = as.numeric(z_values[,grep(colnames(z_values),pattern=toString(patient),fixed=TRUE)]),
                    "path"=metsInset[,"path"],
                    "HMDB"=rownames(metsInset),
                    ints)
  
  ints=ints[order(ints[,"Z.score"]),]
  
  # genExcelFileShort(as.data.frame(ints), paste("./results/crossomics/", gene_in, "/Recon2/P", patient,".xls",sep=""))
  # genExcelFileShort(as.data.frame(ints), paste(path, "/P", patient, "/Recon2/", gene_in,".xls",sep=""))
  # genExcelFileShort(as.data.frame(ints), paste(path, "/", patient, "/Recon2/", gene_in,".xls",sep=""))
  # genExcelFileShort(as.data.frame(ints), paste(path, "/", patient_folder, "/Recon2/", gene_in,".xls",sep=""))
  genExcelFileShort(as.data.frame(ints), paste(path, "/", patient_folder, "/", gene_in,".xls",sep=""))
  
  return(list("p.value"=p))
  #     
  #   } else {
  #     
  #     return(NULL)
  #     
  #   }
  
}
