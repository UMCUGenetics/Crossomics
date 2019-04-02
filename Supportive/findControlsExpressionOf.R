findControlsExpressionOf <- function(HGNC_id){
# HGNC_id=gene_map_table$hgnc_symbol[j]  
# HGNC_id=hgnc  

  # http://www.pathwaycommons.org/pc2/formats
  
  # controls-expression-of      First protein controls a conversion or a template reaction that changes expression of the second protein.
  # controls-state-change-of    First protein controls a reaction that changes the state of the second protein.
  # controls-transport-of       First protein controls a reaction that changes the cellular location of the second protein.
  # controls-phosphorylation-of First protein controls a reaction that changes the phosphorylation status of the second protein.
  # catalysis-precedes          First protein controls a reaction whose output molecule is input to another reaction controled by the second protein.
  # in-complex-with             Proteins are members of the same complex.
  # interacts-with 	            Proteins are participants of the same MolecularInteraction.
  # neighbor-of 	              Proteins are participants or controlers of the same interaction.
  # consumption-controled-by   	The small molecule is consumed by a reaction that is controled by a protein
  # controls-production-of     	The protein controls a reaction of which the small molecule is an output.
  # controls-transport-of-chemical 	The protein controls a reaction that changes cellular location of the small molecule.
  # chemical-affects           	A small molecule has an effect on the protein state.
  # reacts-with 	              Small molecules are input to a biochemical reaction.
  # used-to-produce 	          A reaction consumes a small molecule to produce another small molecule.
  
  filterInteraction <- function(results, interactionType, HGNC, directional=TRUE){
#     results=searchResults
#     interactionType="in-complex-with"
#     directional=FALSE
    rval=NULL
    tmp = results[which(results$INTERACTION_TYPE == interactionType), ]
    if (directional) {
      tmp = tmp[which(tmp$PARTICIPANT_A == HGNC), ]
      tmp = as.vector(tmp$PARTICIPANT_B)
    } else {
      tmp = c(as.vector(tmp$PARTICIPANT_A), as.vector(tmp$PARTICIPANT_B))
      tmp = tmp[-which(tmp==HGNC)]
    }
    
    if (length(tmp)>0){
      rval=tmp    
    }
    return(rval)
  }
  
  retval = NULL
  
  searchResults = try({graphPc(source = HGNC_id, kind = "neighborhood", format = "BINARY_SIF", verbose = TRUE)}, silent = TRUE) 

  if ("INTERACTION_TYPE" %in% names(searchResults)){
  # if (!is.null(searchResults)){
    retval = c(retval, filterInteraction(searchResults, "controls-expression-of", HGNC_id))
#     retval = c(retval, filterInteraction(searchResults, "controls-state-change-of", HGNC_id))
#     retval = c(retval, filterInteraction(searchResults, "controls-transport-of", HGNC_id))
#     retval = c(retval, filterInteraction(searchResults, "controls-phosphorylation-of", HGNC_id))
#     retval = c(retval, filterInteraction(searchResults, "catalysis-precedes", HGNC_id, FALSE))
#     retval = c(retval, filterInteraction(searchResults, "in-complex-with", HGNC_id, FALSE))
#     retval = c(retval, filterInteraction(searchResults, "interacts-with", HGNC_id, FALSE))
#     retval = c(retval, filterInteraction(searchResults, "neighbor-of", HGNC_id, FALSE))
  }

  return(unique(retval))  
}

#   owl = getPc("http://pathwaycommons.org/pc2/Control_8991c2a81fcb004bda83764c92184d45")
#   saveXML(t1,file="LHX3.owl")
