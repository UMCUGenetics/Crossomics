renameIdentifiers <- function(x){
# x=p.values.assi.pos

  assi = rownames(x)
  labels = rep("", dim(x)[1])
  col.names = colnames(x)
  
  retVal = NULL
  
  # one row for each assignment
  for (i in 1:length(assi)){
    
    # for each met get HMDB identifier
    tmp = as.vector(unlist(strsplit(assi[i], split =";")))
#     hmdb = substr(tmp[-which(tmp=="")], start=1,stop=9) For new identification !!!!!!!!!!!!!!!!!!!!!!!!!!
    hmdb = tmp[-which(tmp=="")]
    
#     # select the one occuring in gene set or else the first
#     hmdb.min.adduct = as.vector(unlist(lapply(hmdb,function(x) as.vector(unlist(strsplit(x, split =" ")))[1])))
    
    # multiplay rows, one for each assignment
    tmp2 = as.vector(x[i,])
    result = NULL

    for (i in 1:length(hmdb)){
      result=rbind(result,tmp2)
    }
    rownames(result)=hmdb
    
    retVal = rbind(retVal,result)

#     index = which(hmdb.min.adduct %in% set)
# 
#     if (length(index)>0){
#       message("Compound not first!")
#       if (length(index)>1){
#         message(i)
#         break
#       }  
#       labels[i] = hmdb[index[1]]
#     } else {
#       labels[i] = hmdb[1]
#     }
  }  
  
  # rownames(x) = labels
  colnames(retVal) = col.names
  
  return(retVal)
}
