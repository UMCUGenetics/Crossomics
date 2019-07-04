##########################################################
# source("https://bioconductor.org/biocLite.R")
# biocLite("paxtoolsr")

ls("package:paxtoolsr")
help(getPc)

library(paxtoolsr)
library(igraph)

# Unambiguously maps, e.g., HGNC gene symbols, NCBI Gene, RefSeq, ENS*, and secondary UniProt identifiers
# to the primary UniProt accessions, or
# ChEBI and PubChem IDs to primary ChEBI. You can mix different standard ID types in one query.
id_uniprot = idMapping(toString(gene_list$HGNC[1]))
id_uniprot = idMapping(48)
id_uniprot = "P00439"
 # or
#id_uniprot = idMapping(gene_list$Gene_ID[7])

searchResults = searchPc(q = toString(id_uniprot), organism = "homo sapiens", type = "Catalysis")
#saveXML(searchResults, file="PAH.xml")

# protein = xpathSApply(searchResults, "/searchResponse/searchHit/name", xmlValue)
dataSource = xpathSApply(searchResults, "/searchResponse/searchHit/dataSource", xmlValue)
uris = xpathSApply(searchResults, "/searchResponse/searchHit/uri", xmlValue)
name = xpathSApply(searchResults, "/searchResponse/searchHit/name", xmlValue)

write.table(as.data.frame(cbind(dataSource, uris, name)), "SLC16A1_gene_protein.txt", sep="\t")

# humancyc ##########################################################################################
# uri = uris[grep("humancyc", dataSource, fixed = TRUE)]
# xml = getPc(uri)
# saveXML(xml, file="SLC16A1_humancyc.OWL")

# reactome ##########################################################################################
uri = uris[grep("reactome_human", dataSource, fixed = TRUE)]
xml = getPc(uri[3])
saveXML(xml, file="SLC16A1_reactome_human_3.OWL")

# 1 ###############
xml = traverse(uri = uri[1], path = "Protein/xref:UnificationXref")

# Extract all the URIs
value_uri = xpathSApply(xml, "//value/text()", xmlValue)
# xml = getPc(uri)
xml = traverse(uri = value_uri, path = "UnificationXref/id")

reactome_rxn = xpathSApply(xml, "//value/text()", xmlValue)
searchResults = searchPc(q = reactome_rxn, organism = "homo sapiens")

paste("http://identifiers.org/reactome",reactome_rxn,sep="/")
# ToDo!!!!!!!!!!!!!


searchResults = searchPc(q = kegg_rxn, organism = "homo sapiens", type = "pathway")

dataSource = xpathSApply(searchResults, "/searchResponse/searchHit/dataSource", xmlValue)
uris = xpathSApply(searchResults, "/searchResponse/searchHit/uri", xmlValue)
name = xpathSApply(searchResults, "/searchResponse/searchHit/name", xmlValue)

cbind(dataSource, uris, name)
write.table(as.data.frame(cbind(dataSource, uris, name)), "PAH_pathway.txt", sep="\t")

xml = getPc(uris[1])
saveXML(xml, file="PAH_KEGG_phenylalanine_tyrosine_tryptophan_biosynthesis.OWL")


searchResults = searchPc(q = kegg_rxn, organism = "homo sapiens", type = "BiochemicalReaction")



# eerst in KEGG kijken ##############################################################################
uri = uris[grep("kegg", dataSource, fixed = TRUE)]

# saveXML(getPc(uri), file="PAH_KEGG.OWL")

xml = traverse(uri = uri, path = "Protein/xref:RelationshipXref")

# Extract all the URIs
value_uris = xpathSApply(xml, "//value/text()", xmlValue)

uri = value_uris[grep("KEGG_Reaction", value_uris, fixed = TRUE)]

# xml = getPc(uri)

xml = traverse(uri = uri, path = "RelationshipXref/id")

kegg_rxn = xpathSApply(xml, "//value/text()", xmlValue)

searchResults = searchPc(q = kegg_rxn, organism = "homo sapiens", type = "pathway")

dataSource = xpathSApply(searchResults, "/searchResponse/searchHit/dataSource", xmlValue)
uris = xpathSApply(searchResults, "/searchResponse/searchHit/uri", xmlValue)
name = xpathSApply(searchResults, "/searchResponse/searchHit/name", xmlValue)

cbind(dataSource, uris, name)
write.table(as.data.frame(cbind(dataSource, uris, name)), "PAH_pathway.txt", sep="\t")

xml = getPc(uris[1])
saveXML(xml, file="PAH_KEGG_phenylalanine_tyrosine_tryptophan_biosynthesis.OWL")


searchResults = searchPc(q = kegg_rxn, organism = "homo sapiens", type = "BiochemicalReaction")


#"ModificationFeature/featureLocation:SequenceSite/sequencePosition"

# visualization ##########################################################################################

sif <- toSif("PAH_KEGG_phenylalanine_metabolism.OWL")

# graph.edgelist requires a matrix
g <- graph.edgelist(as.matrix(sif[, c(1, 3)]), directed = FALSE)
plot(g, layout = layout.fruchterman.reingold)

help("system.file")

# neighbourhood ####################################################################################

gene <- toString(gene_list$HGNC[1])
t1 <- results <- graphPc(source = gene, kind = "neighborhood", format = "BINARY_SIF")
t2 <- t1[which(t1[, 2] == "controls-state-change-of"), ]
t2 <- t1[which(t1[, 2] == "reacts-with"), ]
t2 <- t1[which(t1[, 2] == "chemical-affects"), ]
t2 <- t1[which(t1[, 2] == "consumption-controlled-by"), ]
t2 <- t1[which(t1[, 2] == "neighbor-of"), ]
t2 <- t1[which(t1[, 2] == "used-to-produce"), ]
t2 <- t1[which(t1[, 2] == "controls-expression-of"), ]
t2 <- t1[which(t1[, 2] == "controls-production-of"), ]
t2 <- t1[which(t1[, 2] == "controls-phosphorylation-of"), ]
t2 <- t1[which(t1[, 2] == "interacts-with"), ]

g <- graph.edgelist(as.matrix(t2[, c(1, 3)]), directed = FALSE)
plot(g, layout = layout.fruchterman.reingold)

# neighbourhood ####################################################################################

t1 <- results <- graphPc(source = "http://purl.org/pc2/7/Protein_42cba79fa531c997961ac505c6c4ccfd", kind = "neighborhood", format = "BINARY_SIF")
t2 <- t1[which(t1[, 2] == "controls-production-of"), ]
t2 <- t1[which(t1[, 2] == "used-to-produce"), ]
t2 <- t1[which(t1[, 2] == "consumption-controlled-by"), ]

g <- graph.edgelist(as.matrix(t2[, c(1, 3)]), directed = FALSE)
plot(g, layout = layout.fruchterman.reingold)

saveXML(getPc("http://purl.org/pc2/7/Protein_42cba79fa531c997961ac505c6c4ccfd"), file="PAH_KEGG_Protein.OWL")

for (i in 1:dim(gene_list)[1]){
  left=getMetsPathwayCommons(gene_list$Gene_ID[i], TRUE)
  write.table(left, paste(gene_list$HGNC[i], "left.txt",sep="_"), sep="\t")
  right=getMetsPathwayCommons(gene_list$Gene_ID[i], FALSE)
  write.table(right, paste(gene_list$HGNC[i], "right.txt",sep="_"), sep="\t")
}
