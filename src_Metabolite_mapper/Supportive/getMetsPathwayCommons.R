getMetsPathwayCommons <- function(entrezgene){
#  entrezgene=gene_in
  
  require("paxtoolsr")
  
  id_uniprot = idMapping(toString(entrezgene))
  
  searchResults = NULL

  try({searchResults = searchPc(q = toString(id_uniprot), organism = "homo sapiens", type = "Catalysis")}, silent = TRUE)
  path = "Catalysis/controlled"
  if (is.null(searchResults)){
    message("No search results for Catalysis, try Control!!!")
    try({searchResults = searchPc(q = toString(id_uniprot), organism = "homo sapiens", type = "Control")}, silent = TRUE)
    path = "Control/controlled"
    if (is.null(searchResults)){
      return(NULL)
    }  
  }
  
  uris = xpathSApply(searchResults, "/searchResponse/searchHit/uri", xmlValue)

#   owl = getPc("http://pathwaycommons.org/pc2/Control_11b881ab1d6a32b7c22733df937a8e48")
#   saveXML(owl,file="Control_11b881ab1d6a32b7c22733df937a8e48.owl")
  
  
  # Dit werkt niet, finally wordt altijd uitgevoerd!
  # First try at ones
  result = tryCatch(
    { xml = traverse(uri = uris, path = path)
      if (!is.null(xml)){
        uris2 = xpathSApply(xml, "//value/text()", xmlValue)
      }
    }
    , error = function(e) {
      message(paste("CATCHED", e))
    }
    , finally = {
      # If error try in loop
      uris2 = NULL
      xml=NULL
      for (i in 1:length(uris)){
        try({xml = traverse(uri = uris[i], path = path)}, silent = TRUE)
        if (!is.null(xml)){
          u = xpathSApply(xml, "//value/text()", xmlValue)
          uris2 = c(uris2,u)
        }
      }
    })

  # first try biochemical reaction
  xml = NULL
  try({xml = traverse(uri = uris2, path = "BiochemicalReaction/left:SmallMolecule")}, silent = TRUE)
  
  if (is.null(xml)){
    
    message("No Biochemical reaction, try Transport...")
    
    try({xml = traverse(uri = uris2, path = "Transport/left:SmallMolecule")}, silent = TRUE)
    
    if (is.null(xml)){
      message("No BiochemicalReaction nor Transport, stop searching")
      return(NULL)
      
    } else {
      uris3 = xpathSApply(xml, "//value/text()", xmlValue)
      leftMets = getLeftOrRight(uris3, entrezgene, "left")
      # For transport right and left mets are the same!
      rightMets = leftMets
      
#       try({xml = traverse(uri = uris, path = "Transport/right:SmallMolecule")}, silent = TRUE)
#       if (is.null(xml)){
#         message("No rigth small molecules for Transport!!!")
#         #return(NULL)
#       } else {
#         uris3 = xpathSApply(xml, "//value/text()", xmlValue)
#         rightMets = getLeftOrRight(uris3, entrezgene)
#       }
    }

  } else {
    uris3 = xpathSApply(xml, "//value/text()", xmlValue)
    leftMets = getLeftOrRight(uris3, entrezgene, "left")
    
    try({xml = traverse(uri = uris2, path = "BiochemicalReaction/right:SmallMolecule")}, silent = TRUE)
    if (is.null(xml)){
      message("No rigth small molecules for BiochemicalReaction!!!")
      #return(NULL)
    } else {
      uris4 = xpathSApply(xml, "//value/text()", xmlValue)
      rightMets = getLeftOrRight(uris4, entrezgene,"right")
    }
  }

  return(list("left"=leftMets, "right"=rightMets))
  
}