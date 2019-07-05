library("BridgeDbR")
library("XLConnect")

genExcelFileShort <- function(list, wbfile) { 
  #   list=result_from_consumed_plus_produced$mets
  #   wbfile=paste("./results/mets",gene_in,"xls", sep=".")
  
  filelist <- "metabolites"
  wb <- loadWorkbook(wbfile, create = TRUE)
  createSheet(wb, name = filelist)
  writeWorksheet(wb, list, sheet = filelist)
  saveWorkbook(wb)
  rm(wb)
}


# A lot of missing HMDB identifiers!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
metaboliteSet = mss_2
if (!is.null(metaboliteSet)){
  if (is.null(dim(metaboliteSet))){
    metaboliteSet=data.frame(t(metaboliteSet))
  }
}  

##################################################
# temporarily work around to be fixed in findMetabolicEnvironment
index = which(is.na(metaboliteSet[,"hmdb"]))
if (length(index)>0) metaboliteSet[index,"hmdb"] = "character(0)"
##################################################

########### Recon 2.0 ############################
index = which(metaboliteSet[,"hmdb"] == "character(0)") 

# pressent in BridgeDB?
index.sub = which(metaboliteSet[index,"kegg"] != "character(0)")
kegg_id = metaboliteSet[index[index.sub],"kegg"]

mapper = loadDatabase("./db/BridgeDB/metabolites_20150717.bridge")
hmdb = getSystemCode("HMDB")
kegg = getSystemCode("KEGG Compound")

if (length(kegg_id)>0){
  for (k in 1:length(kegg_id)){
    if (!is.null(unlist(map(mapper, kegg, kegg_id[k], hmdb)[1]))) metaboliteSet[index[index.sub[k]],"hmdb"] = unlist(map(mapper, kegg, kegg_id[k], hmdb)[1]) 
  }
}  

index = which(metaboliteSet[,"hmdb"] == "character(0)")
index.sub = which(metaboliteSet[index,"chebi"] != "character(0)")
chebi_id = metaboliteSet[index[index.sub],"chebi"]

chebi = getSystemCode("ChEBI")

if (length(chebi_id)>0){
  for (k in 1:length(chebi_id)){
    if (!is.null(unlist(map(mapper, chebi, chebi_id[k], hmdb)[1]))) metaboliteSet[index[index.sub[k]],"hmdb"] = unlist(map(mapper, chebi, chebi_id[k], hmdb)[1]) 
  }
}

index = which(metaboliteSet[,"hmdb"] == "character(0)") 
if (length(index)>0) metaboliteSet = metaboliteSet[-index,,drop=FALSE]
#######################################################

#     ########### Recon 2.2 ############################
#     inchi = metaboliteSet[,"InChI_key"]
#     chebi_id = metaboliteSet[,"chebi"]
#     index=which((is.na(chebi_id)&is.na(inchi)))
#     if (length(index)>0){
#       metaboliteSet = metaboliteSet[-index,]
#     }
#     #######################################################

# Add glycine and carnitine conjugates to CoA compounds
#addConjugates(metaboliteSet)

#     if (!is.null(metaboliteSet) && (dim(metaboliteSet)[1]!=0)){
#       if (is.null(dim(metaboliteSet))){
#         metaboliteSet=data.frame(t(metaboliteSet))
#       }

nMets=c(nMets, dim(metaboliteSet)[1])
message(paste("gene: ", gene_in))
message(paste("metaboliteSet: ", dim(metaboliteSet)[1]))

genExcelFileShort(metaboliteSet, "GLS_metabolite_set.xls")
