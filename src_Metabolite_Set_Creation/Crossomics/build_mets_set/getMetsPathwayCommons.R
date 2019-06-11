getMetsPathwayCommons <- function(id_uniprot, src){
# id_uniprot="Q9UM01"

  # require("paxtoolsr")
  # require("SPARQL")
  
  # library("biomaRt")
  # ls("package:paxtoolsr")

  # ensembl = tryCatch(
  #   { useMart("ensembl",dataset="hsapiens_gene_ensembl") }
  #   , error = function(e) {
  #     message(paste("CATCHED", e))
  #     message("Use other server!")
  #     useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = "grch37.ensembl.org")
  #   })

  # listDatasets(mart = ensembl)
  # listAttributes(ensembl)

  # id_uniprot = idMapping(toString(entrezgene), verbose = TRUE)
  # id <- idMapping("HEXA", verbose = TRUE)
  # searchResults <- searchPc(q = "glycolysis", type = "pathway")
  # id_uniprot ="P10636"

  searchResults = NULL
#####################################################################################  
  
  try({searchResults = searchPc(q = toString(id_uniprot), organism = "homo sapiens", type = "Catalysis")}, silent = TRUE)
  path = "Catalysis/controlled"
  if (grepl('numHits=\"0\"', as(searchResults, "character"))){
  # if (is.null(searchResults)){ #never TRUE, searchResults will always contain an xml header
    message("No search results for Catalysis, try Control!!!")
    try({searchResults = searchPc(q = toString(id_uniprot), organism = "homo sapiens", type = "Control")}, silent = TRUE)
    path = "Control/controlled"
    if (is.null(searchResults)){
      return(NULL)
    }  
  }
  
  uris = xpathSApply(searchResults, "/searchResponse/searchHit/uri", xmlValue)
  
  leftMets=NULL
  rightMets=NULL
  xml = NULL
  
  try({xml = traverse(uri = uris, path = path)}, silent = TRUE)
  if(is.null(xml)){return(NULL)}
  uris2 = xpathSApply(xml, "//value/text()", xmlValue)

  # Transport
  try({xml2 = traverse(uri = uris2, path = "Transport/left:SmallMolecule")}, silent = TRUE)
  uris3 = xpathSApply(xml2, "//value/text()", xmlValue)
  if (!is.null(uris3)){
    leftMets = rbind(leftMets, getLeftOrRight(uris3, entrezgene, "left"))
  }
  try({xml2 = traverse(uri = uris2, path = "Transport/right:SmallMolecule")}, silent = TRUE)
  uris3 = xpathSApply(xml2, "//value/text()", xmlValue)
  if (!is.null(uris3)){
    rightMets = rbind(rightMets, getLeftOrRight(uris3, entrezgene, "right"))
  }
  
  # BiochemicalReaction
  try({xml2 = traverse(uri = uris2, path = "BiochemicalReaction/left:SmallMolecule")}, silent = TRUE)
  uris3 = xpathSApply(xml2, "//value/text()", xmlValue)
  if (!is.null(uris3)){
    leftMets = rbind(leftMets, getLeftOrRight(uris3, entrezgene, "left"))
  }
  try({xml2 = traverse(uri = uris2, path = "BiochemicalReaction/right:SmallMolecule")}, silent = TRUE)
  uris3 = xpathSApply(xml2, "//value/text()", xmlValue)
  if (!is.null(uris3)){
    rightMets = rbind(rightMets, getLeftOrRight(uris3, entrezgene, "right"))
  }

  # # owl = getPc(uris[2])
  # # saveXML(owl,file="reactome.owl")
  # 
  # # Chebi parent to child relaties
  # chebi = read.csv("compounds.tsv",header=TRUE, sep="\t")
  # 
  # entry = which(chebi[,"CHEBI_ACCESSION"]=="CHEBI:15705")
  # chebi[which(chebi[,"PARENT_ID"]==chebi[entry,"ID"]), ]
  # 
  # entry = which(chebi[,"CHEBI_ACCESSION"]=="CHEBI:16467")
  # chebi[which(chebi[,"PARENT_ID"]==chebi[entry,"ID"]), ]
  
  
  # , ns=c("time","<http://www.w3.org/2006/time#>")
  

  # http://purl.obolibrary.org/obo/chebi.owl
  
  # endpoint="http://chebi.bio2rdf.org/sparql"
  # 
  # 
  # endpoint="http://dbpedia.org/sparql"
  # d = SPARQL(url=endpoint,
  #            query="SELECT DISTINCT ?concept WHERE {?subject a ?concept .} LIMIT 50")
  # 
  # 
  # d = SPARQL(url=endpoint,
  #   query="PREFIX chebi1: <http://purl.obolibrary.org/obo/chebi#2>
  #          PREFIX chebi2: <http://purl.obolibrary.org/obo/chebi#3>
  #          PREFIX chebi3: <http://purl.obolibrary.org/obo/chebi#1>
  #          PREFIX chebi: <http://purl.obolibrary.org/obo/chebi#>
  #          PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
  #          PREFIX owl: <http://www.w3.org/2002/07/owl#>
  #          PREFIX oboInOwl: <http://www.geneontology.org/formats/oboInOwl#>
  #          PREFIX xml: <http://www.w3.org/XML/1998/namespace>
  #          PREFIX xsd: <http://www.w3.org/2001/XMLSchema#>
  #          PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
  #          PREFIX obo: <http://purl.obolibrary.org/obo/>
  #          SELECT * WHERE { ?IAO_0000115 oboInOwl:id CHEBI:15705.} LIMIT 10")
  # is.data.frame(d$results)
                   

  resources = unique(c(leftMets[,"Data Source"],rightMets[,"Data Source"]))
  rxns = NULL
  for (i in 1:length(resources)){
    rxn = paste(paste(leftMets[which(leftMets[,"Data Source"]==resources[i]), "Compound"], collapse = " + "),
                paste(rightMets[which(rightMets[,"Data Source"]==resources[i]), "Compound"], collapse = " + "),sep=" <=> ")
    names(rxn) = resources[i]
    rxns = c(rxns,rxn)
  }
  
  return(list("left"=leftMets, "right"=rightMets, "rxns"=rxns))
}