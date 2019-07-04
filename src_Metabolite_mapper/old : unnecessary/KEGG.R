source("https://bioconductor.org/biocLite.R")
biocLite("KEGGREST")

library("KEGGREST")


# pathways gene 48 is in
pathways = keggLink("pathway", "hsa:48") 
rxns = keggList("reaction", pathways[1])


pathways[1] # "path:hsa00020"
# option = c("aaseq", "ntseq", "mol", "kcf", "image", "kgml"))
keggGet(pathways[1]) 
  
# reactions in pathway
keggLink("reaction", "hsa:00020") #"path:map00010"

pathway = keggLink("pathway", "hsa:48") 


pathways = keggLink("pathway", "hsa:48")
pathway=keggGet(pathways[1])

rxns = keggLink("reaction", "hsa")
rxns = keggLink("reaction", "hsa:48")
rxns = keggLink("reaction", "hsa:197322")
query = keggGet(rxns)
query[[1]]$ENTRY
query[[1]]$NAME 
query[[1]]$DEFINITION
query[[1]]$ENZYME
#keggConv(query[[1]]$ENTRY, "")

cmpnd = keggLink("compound", rxns[1])
query = keggGet(cmpnd)
query[[1]]$NAME
query[[2]]$NAME

query[[1]]$DBLINKS[3]

query[[1]]$ENTRY

############################################
ls("package:BiGGR")
library("BiGGR")
library("sybilSBML")
data("Recon2")
sbml.model = buildSBMLFromGenes(gene_list$Gene_ID[1], Recon2)
createLIMFromSBML(sbml.model, maximize=NULL, file.name="PAH_env.sbml")
m=readSBMLmod(filename="PAH_env.sbml")
validateSBMLdocument('PAH_env.sbml')


readSBMLmod(filename, description,
            def_bnd = SYBIL_SETTINGS("MAXIMUM"),
            validateSBML = FALSE,
            extMetFlag = "b",
            bndCond = TRUE,
            ignoreNoAn = FALSE,
            mergeMet = TRUE,
            balanceReact = TRUE,
            remUnusedMetReact = TRUE,
            singletonMet = FALSE,
            deadEndMet = FALSE,
            remMet = FALSE,
            constrMet = FALSE,
            tol = SYBIL_SETTINGS("TOLERANCE"))

gene.info <- extractGeneAssociations(Recon2)
ls("package:BiGGR")
###############################################



getOrganismInfo <- function(organism){
  KEGG_INFO_BASE <- "http://rest.kegg.jp/info/"
  info_REST_url <- paste(KEGG_INFO_BASE, organism, sep="")
  info <- readLines(info_REST_url)
  info
}

getOrganismInfo("hsa")

mapPathwayToName <- function(organism) {
  KEGG_PATHWAY_LIST_BASE <- "http://rest.kegg.jp/list/pathway/"
  pathway_list_REST_url <- paste(KEGG_PATHWAY_LIST_BASE, organism, sep="")
  
  pathway_id_name <- data.frame()
  
  for (line in readLines(pathway_list_REST_url)) {
    tmp <- strsplit(line, "\t")[[1]]
    pathway_id <- strsplit(tmp[1], organism)[[1]][2]
    pathway_name <- tmp[2]
    pathway_name <- strsplit(pathway_name, "\\s+-\\s+")[[1]][1]
    pathway_id_name[pathway_id, 1] = pathway_name
    
  }
  
  names(pathway_id_name) <- "pathway_name"
  pathway_id_name
}

rval=mapPathwayToName("hsa")

mapGeneToPathway <- function(organism) {
  KEGG_PATHWAY_LINK_BASE <- "http://rest.kegg.jp/link/pathway/"
  pathway_link_REST_url <- paste(KEGG_PATHWAY_LINK_BASE, organism, sep="")
  
  gene_pathway <- data.frame()
  
  for (line in readLines(pathway_link_REST_url)) {
    tmp <- strsplit(line, "\t")[[1]]
    gene <- tmp[1]
    gene <- strsplit(gene, ":")[[1]][2]  
    pathway_id<- strsplit(tmp[2], organism)[[1]][2]
    
    if (is.null(gene_pathway[gene, 1])) {
      gene_pathway[gene,1] = pathway_id
    } else {
      if (is.na(gene_pathway[gene,1])) {
        gene_pathway[gene,1] = pathway_id
      } else {
        gene_pathway[gene,1] = paste(gene_pathway[gene, 1], pathway_id, sep=";")
      }
    }
  }
  names(gene_pathway) <- "pathway_id"
  gene_pathway
}

rval=mapGeneToPathway("hsa")

mapGeneToReaction <- function(organism) {
  organism="hsa"
  KEGG_REACTION_LINK_BASE <- "http://rest.kegg.jp/link/rn/"
  reaction_link_REST_url <- paste(KEGG_REACTION_LINK_BASE, organism, sep="")
  
  gene_reaction <- data.frame()
  
  readLines("http://rest.kegg.jp/link/pathway/hsa:48")
  readLines("http://rest.kegg.jp/list/rn/map00010")
  
  readLines("http://rest.kegg.jp/get/hsa:48")
  readLines("http://rest.kegg.jp/find/reaction")
  
  for (line in readLines(reaction_link_REST_url)) {
    tmp <- strsplit(line, "\t")[[1]]
    gene <- tmp[1]
    gene <- strsplit(gene, ":")[[1]][2]  
    pathway_id<- strsplit(tmp[2], organism)[[1]][2]
    
    if (is.null(gene_pathway[gene, 1])) {
      gene_reaction[gene,1] = reaction_id
    } else {
      if (is.na(gene_reaction[gene,1])) {
        gene_reaction[gene,1] = reaction_id
      } else {
        gene_reaction[gene,1] = paste(gene_pathway[gene, 1], reaction_id, sep=";")
      }
    }
  }
  names(gene_reaction) <- "reaction_id"
  gene_reaction
}

rval=mapGeneToReaction("hsa")

###########################################################
# RbioRXN
# biocLite("fmcsR")
library("RbioRXN")
rxn=get.kegg.byId("R01324")
names(rxn)
rxn$DEFINITION
rxn$EQUATION
rxn$ENZYME
ls("package:RbioRXN")

kegg=get.kegg.all()
save(kegg, file = "kegg.RbioRXN.RData")
dim(kegg)
names(kegg)

