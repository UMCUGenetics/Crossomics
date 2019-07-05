# Used if reaction gene mapping doesn't work in Recon2
getMetsKEGG <- function(entrezgene, consumedMets){
# entrezgene=5053
# consumedMets = TRUE

  retVal=NULL
  
  rxns4gene = keggLink("reaction", paste("hsa", entrezgene, sep=":"))

  if (length(rxns4gene)!=0){
  
    n = ceiling(length(rxns4gene)/10)
    for (h in 1:n){
      end=h*10
      start=end-9
      if (end>length(rxns4gene)) end=length(rxns4gene)
      rxns4geneSub = rxns4gene[start:end]
      
      if (length(rxns4geneSub)>0){
        query = keggGet(rxns4geneSub)
        
        for (i in 1:length(rxns4geneSub)){
          message(query[[i]]$DEFINITION)
          if (consumedMets) {
            tmp = strsplit(query[[i]]$DEFINITION, " <=> ", fixed = T)[[1]][1]
          } else {
            tmp = strsplit(query[[i]]$DEFINITION, " <=> ", fixed = T)[[1]][2]
          }
          cmpnds_consumed = strsplit(tmp, " + ", fixed = T)[[1]][]
          
          cmpnd = keggLink("compound", rxns4geneSub[i])
          query2 = keggGet(cmpnd)
          
          for (j in 1:length(cmpnd)){
            
            synnonyms = query2[[j]]$NAME
            # remove semicolon
            for (k in 1:length(synnonyms)){
              if (substr(synnonyms[k], nchar(synnonyms[k]), nchar(synnonyms[k]))==";"){
                synnonyms[k] = substr(synnonyms[k], 1, nchar(synnonyms[k])-1)     
              }
            }
            
            if (length(which(synnonyms %in% cmpnds_consumed))>0) {
              
              grep("PubChem", as.vector(query2[[j]]$DBLINKS))
              
              retVal = rbind(retVal, c(synnonyms[1], query2[[j]]$ENTRY, query2[[j]]$DBLINKS[grep("ChEBI", as.vector(query2[[j]]$DBLINKS))], query2[[j]]$DBLINKS[grep("PubChem", as.vector(query2[[j]]$DBLINKS))]))
            }
          }
        }  
        
        colnames(retVal) = c("Compound","KEGG_id","CheBI_id", "PubChem_id")
        
      } else {
        message("No reactions found in KEGG!") 
      }
    }
  
  } else {
    message("No reactions found in KEGG!")
  }
  
  return(retVal)   
  
}    
