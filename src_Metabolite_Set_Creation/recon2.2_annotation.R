
library(rsbml)
# doc <- rsbml_read(system.file("sbml", "GlycolysisLayout.xml", package = "rsbml"), dom = FALSE)
#doc <- rsbml_read("./src/MODEL1603150001.xml", dom = TRUE)
# doc@model@species$M_10fthf6glu_c@annotation

# doc=model
# model=NULL
# save(doc, file=paste(src, "MODEL1603150001_rsbml.RData", sep="/"))

load("./src/MODEL1603150001_rsbml.RData")

mets=sapply(species(model(doc)),id)
annot=sapply(species(model(doc)),annotation)

# head(as.vector(mets))
# length(as.vector(annot))

chebi=NULL
met_id=NULL
InChI=NULL

for (i in 1:length(as.vector(mets))){

  tmp=unlist(strsplit(as.vector(mets[i]), "_", fixed = TRUE))
  l=length(tmp)
  localization=tmp[l]
  tmp=tmp[-c(1,l)]
  tmp=paste(tmp,collapse = "_")
  # tmp=gsub("_", "", tmp)
  tmp=paste(tmp,"[",localization,"]",sep="")
  met_id=c(met_id,tmp)

  rdf=as.vector(unlist(annot[i]))

  if (length(rdf)==0){
    chebi=c(chebi,NA)
    InChI=c(InChI,NA)
  } else {

    tmp=unlist(strsplit(unlist(strsplit(rdf, "http://identifiers.org/chebi/CHEBI:", fixed = TRUE))[2], "\"", fixed = TRUE))[1]
    if (is.na(tmp)){
      chebi=c(chebi,NA)
    } else {
      chebi=c(chebi,tmp)
    }

    tmp=unlist(strsplit(unlist(strsplit(rdf, "http://identifiers.org/inchi/", fixed = TRUE))[2], "\"", fixed = TRUE))[1]
    if (is.na(tmp)){
      InChI=c(InChI,NA)
    } else {
      InChI=c(InChI,tmp)
    }
  } 
}

recon2chebi = cbind("species_id"=met_id, "CHEBI"=chebi, "InChI"=InChI)
save(recon2chebi, file="./src/recon2chebi_MODEL1603150001.RData")

# write.table(cbind("species_id"=as.vector(mets), "CHEBI"=chebi), file="./src/recon2chebi_MODEL1603150001.txt", sep="\t")




# require("paxtoolsr")


