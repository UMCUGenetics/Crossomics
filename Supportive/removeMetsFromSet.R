removeMetsFromSet <- function(result_mets, model, index_consumed=NULL, index_produced=NULL){
#result_mets =  metsExtend
  
  mets2remove =  c(
      "O2",
      "Water",
      "ATP","ADP","AMP",
      "GTP","GDP","GMP",
      "CTP","CDP","CMP",
      "TTP","TDP","TMP",
      "UTP","UDP","UMP",
      "carbon dioxide",
      "hydrogenphosphate",
      "Diphosphate",
      "hydrosulfide",
      "Nicotinamide adenine dinucleotide phosphate",
      "Nicotinamide adenine dinucleotide","proton",
      "Nicotinamide adenine dinucleotide - reduced",
      "Nicotinamide adenine dinucleotide phosphate - reduced",
      "Nicotinamide adenine dinucleotide",
      "proton",
      "Flavin adenine dinucleotide oxidized",
      "Flavin adenine dinucleotide reduced",
      "Coenzyme A")
  
  #candidates:
#   "Sodium",
#   "Chloride",
#   "Hydrogen peroxide",
#   "potassium",
#   "calcium(2+)"
  
  for (i in 1:length(mets2remove)){

    if (is.null(dim(result_mets))){
      index = grep(mets2remove[i], result_mets["met_long"])
    } else {
      index = grep(mets2remove[i], result_mets[,"met_long"])
    }
    
    if (length(index)>0){

      result_mets = result_mets[-index,,drop=FALSE]
      
      index = grep(mets2remove[i], as.vector(unlist(model$metNames[index_produced])))
      if (length(index)>0){
        index_produced = index_produced[-index]
      }

      index = grep(mets2remove[i], as.vector(unlist(model$metNames[index_consumed])))
      if (length(index)>0){
        index_consumed = index_consumed[-index]
      }
    }
  }
  
  return(list("result_mets"=result_mets, "index_produced"=index_produced, "index_consumed"=index_consumed))
  
}