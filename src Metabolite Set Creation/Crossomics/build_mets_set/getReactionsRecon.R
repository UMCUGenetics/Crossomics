getReactionsRecon <- function(index, model){
# index=tmp$mets

  thresh = 25
  labels = NULL
  result = NULL
  
  # leave this message it avoids a bug!!!!!!!!!!!!!!!!
  # message(dim(model$S))
  index2do=index
  
  index2do=removeMetsFromSet(index, model)
  
  for (i in 1:length(index)){
    # i=1
    
    labels = c(labels, convert(model$metNames[index[i]]))

    if (!index[i]%in%index2do) {
      result = c(result,paste("Filtered out"))
      # message("next")
      next
    }

    rxns_index = which(model$S[index[i],]!=0)
    
    rxns = NULL
    truncate=FALSE
    end = length(rxns_index)
    if (length(rxns_index)>thresh){
      truncate=TRUE
      end = thresh
    }

    # reactions
    for (j in 1:end){
      left_index = which(model$S[,rxns_index[j]]<0)
      stoichiometry = abs(model$S[left_index, rxns_index[j]])
      left = paste(paste(stoichiometry, convert(model$metNames[left_index])), collapse =" + ")
      left = gsub("1 ","",left)

      right_index = which(model$S[,rxns_index[j]]>0)
      stoichiometry = abs(model$S[right_index, rxns_index[j]])
      right = paste(paste(stoichiometry, convert(model$metNames[right_index])), collapse =" + ")
      right = gsub("1 ","",right)

      rxns = c(rxns, paste(left, right, sep = " => "))

    }
    
    if (truncate) {
      result = c(result, paste(paste(rxns, collapse = "; "),";...<tuncated>",sep=""))
    } else {
      result = c(result, paste(rxns, collapse = "; "))
    }
  }
  
  # names(result) = labels
  
  return(cbind("met"=labels,"rxns"=result))
}