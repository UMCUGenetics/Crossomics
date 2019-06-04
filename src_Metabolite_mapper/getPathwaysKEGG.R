library("BridgeDbR")

#location = getDatabase("Homo sapiens",location=getwd())
# manually donwnloaded
mapper = loadDatabase("./db/BridgeDB/metabolites_20150717.bridge")
hmdb = getSystemCode("HMDB")
#getMatchingSources("C00010")
kegg = getSystemCode("KEGG Compound")

pathways = keggLink("hsa", "pathway")
length(pathways)

# all pathways: from path:hsa00010 to path:hsa05416

for (i in 1:5416) {
#for (i in 1:100) {  
  
#  i=20 # TCA
  pathwayID = paste("path:hsa", sprintf("%05d", i), sep="")
  #message(pathwayID)
#  pathwayID ="path:hsa01100"
  index = which(names(pathways)==pathwayID) 
  
  if (length(index)>0){
    
    query = keggGet(pathwayID)
    
    if (is.null(query[[1]]$COMPOUND)) next
#     query[[1]]$NAME
#     query[[1]]$COMPOUND
#     query[[1]]$DISEASE

    allMets4pathway = cbind("KEGG_id"=names(query[[1]]$COMPOUND), "name"=as.vector(query[[1]]$COMPOUND), "hmdb"=as.vector(unlist(lapply(names(query[[1]]$COMPOUND), function(x) {
      if (length(map(mapper, kegg, x, hmdb))==0) NA else map(mapper, kegg, x, hmdb)[1]}))))
    
    
    # Look in Recon2
    index.na = which(is.na(allMets4pathway[,"hmdb"]))
    if (length(index.na)==1){
      notFoundInBridgeDB = allMets4pathway[index.na,]
      notFoundInBridgeDB = c(notFoundInBridgeDB, "CheBI_id"=rep(NA,length(index.na)))
      notFoundInBridgeDB = rbind(notFoundInBridgeDB, notFoundInBridgeDB) # dirty work around
    } else if (length(index.na)>1){
      notFoundInBridgeDB = allMets4pathway[index.na,]
      notFoundInBridgeDB = cbind(notFoundInBridgeDB, "CheBI_id"=rep(NA,length(index.na)))
    }

    tmp.index = getMets2ReconID(notFoundInBridgeDB)
    if (length(tmp.index)>0){
      # message("Found in recon2") 
      for (j in 1:length(tmp.index)){
        if (length(as.vector(unlist(model$metHMDB[tmp.index[j]])))>0) {
          allMets4pathway[which(allMets4pathway %in% as.vector(unlist(model$metKeggID[tmp.index[j]]))),"hmdb"]=as.vector(unlist(model$metHMDB[tmp.index[j]]))
        }
      }
    }
    
    #   
    #   # replace empty entries by NA
    #   for (j in 1:length(tmp.index)){
    #     if (length(as.vector(unlist(model$metHMDB[tmp.index[j]])))==0) model$metHMDB[tmp.index[j]][[1]]=list(matrix(NA))
    #     if (length(as.vector(unlist(model$metKeggID[tmp.index[j]])))==0) model$metKeggID[tmp.index[j]][[1]]=list(matrix(NA))
    #   }
    # 
    #   kegg2hmdb = cbind("KEGG"=as.vector(unlist(model$metKeggID[tmp.index])),
    #                     "HMDB"=as.vector(unlist(model$metHMDB[tmp.index])))
    #   kegg2hmdb = unique(kegg2hmdb)
    #   
    #   allMets4pathway = cbind(allMets4pathway, "hmdb"=rep(NA, dim(allMets4pathway)[1]))
    #   allMets4pathway = cbind(allMets4pathway, "check"=rep(NA, dim(allMets4pathway)[1]))
    #   
    #   index.1=which(kegg2hmdb[,"KEGG"] %in% allMets4pathway[,"KEGG_id"])
    #   kegg2hmdbInPathway = kegg2hmdb[index.1,]
    #   rownames(kegg2hmdbInPathway)=kegg2hmdbInPathway[,"KEGG"]
    #   
    #   index.2=which(allMets4pathway[,"KEGG_id"] %in% kegg2hmdbInPathway[,"KEGG"])
    #   allMets4pathway[index.2,]
    #   
    #   allMets4pathway[index.2,"HMDB"]=kegg2hmdbInPathway[allMets4pathway[index.2,"KEGG_id"],"HMDB"]
    #   allMets4pathway[index.2,"check"]=kegg2hmdbInPathway[allMets4pathway[index.2,"KEGG_id"],"KEGG"]
    
    pathway=list("name"=query[[1]]$NAME, "metabolites"=allMets4pathway, "diseases"=query[[1]]$DISEASE)
    save(pathway, file=paste("./db/KEGG/", gsub(":", "_", pathwayID), ".RData", sep = ""))
  }

}



