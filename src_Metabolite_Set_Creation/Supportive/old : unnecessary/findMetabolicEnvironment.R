findMetabolicEnvironment <- function(gene_in, model, recon2chebi, HGNC, step, reconVersion, src, rval){
# gene_in=hgnc
# HGNC=hgnc

  message(HGNC)

  if (reconVersion==2.2){
    S=model@S
  } else if (reconVersion==2.0){
    S=model$S
  }

  max_rxns=4

  rvalMets=NULL
  rvalRxns=NULL
  index_consumed=NULL
  index_produced=NULL

  # result = getMetsRecon2(59272, -1, "ser_L[c]")
  #   consumed = getMetsKEGG(gene_in,consumedMets=TRUE)
  #   produced = getMetsKEGG(gene_in,consumedMets=FALSE)

  # # Get primary reaction
  # rval = getMetsPathwayCommons(HGNC, src)
  # return(rval)
  
  if (is.null(rval)) {
    return(NULL)
  } else {

    consumed = rval$left
    produced = rval$right
    rxns = rval$rxns

    if (length(consumed)>0) {
      # Do not convert to ReconID in primary reaction!
      amount = dim(consumed)[1]
      result_mets_consumed = cbind("rxn_id"=rep(toString(gene_in),amount),
                                   "step"=rep(0,amount),
                                   "met_in"=rep(NA,amount),
                                   "left_right"=rep("left",amount),
                                   "met_short"=consumed[,"Compound"],
                                   "met_long"=consumed[,"Compound"],
                                   "hmdb"=rep(NA,amount),
                                   "kegg"=consumed[,"KEGG_id"],
                                   "chebi"=consumed[,"CheBI_id"],
                                   "pubchem"=consumed[,"PubChem_id"],
                                   "InChI_key"=consumed[,"InchI_key"],
                                   "rxn_name"=rep(NA,amount),
                                   "rxn"=rep(NA,amount),
                                   "resource"=consumed[,"Data Source"])
      tmp=NULL
      for (i in 1:dim(result_mets_consumed)[1]){
        source = rxns[grep(result_mets_consumed[i,"resource"], names(rxns), fixed=TRUE)]
        tmp=c(tmp,
        paste(source[grep(result_mets_consumed[i,"met_short"], source, fixed = TRUE)],collapse = ";"))
      }
      result_mets_consumed = cbind(result_mets_consumed,"rxn_formula"=tmp)

      if (step>0) {
        if (reconVersion==2.2){

          index_consumed = getMets2Recon2.2ID(consumed, recon2chebi)

          amount = length(convert(model@met_id[index_consumed]))

          result_mets_consumed = rbind(result_mets_consumed,
                                       cbind("rxn_id"=rep(toString(gene_in),amount),
                                             "step"=rep(0,amount),
                                             "met_in"=rep(NA,amount),
                                             "left_right"=rep("left",amount),
                                             "met_short"=convert(model@met_id[index_consumed]),
                                             "met_long"=convert(model@met_name[index_consumed]),
                                             "hmdb"=rep(NA,amount),
                                             "kegg"=rep(NA,amount),
                                             "chebi"=convert(recon2chebi[index_consumed,2]),
                                             "pubchem"=rep(NA,amount),
                                             "InChI_key"=convert(recon2chebi[index_consumed,3]),
                                             "rxn_name"=rep(NA,amount),
                                             "rxn"=rep(NA,amount)))
          # ,
          # "resource"=rep("Recon",amount),
          # "rxn_formula"=getReactionsRecon(index_consumed, model)))
        } else if (reconVersion==2.0){

          index_consumed = getMets2ReconID(consumed, model)

          # replace empty entries by NA
          for (j in 1:length(index_consumed)){
            if (length(convert(model$metHMDB[index_consumed[j]]))==0) model$metHMDB[index_consumed[j]][[1]]=list(matrix(NA))
          }

          amount = length(convert(model$mets[index_consumed]))
          rxn=NA
          if (length(index_consumed)>0) rxn = getReactionsRecon(index_consumed, model)
          
          result_mets_consumed = rbind(result_mets_consumed,
                                       cbind("rxn_id"=rep(toString(gene_in),amount),
                                             "step"=rep(0,amount),
                                             "met_in"=rep(NA,amount),
                                             "left_right"=rep("left",amount),
                                             "met_short"=convert(model$mets[index_consumed]),
                                             "met_long"=convert(model$metNames[index_consumed]),
                                             "hmdb"=convert(model$metHMDB[index_consumed]),
                                             "kegg"=convert(model$metKeggID[index_consumed]),
                                             "chebi"=convert(model$metCHEBIID[index_consumed]),
                                             "pubchem"=convert(model$metPubChemID[index_consumed]),
                                             "InChI_key"=rep(NA,amount),
                                             "rxn_name"=rep(NA,amount),
                                             "rxn"=rep(NA,amount),
                                             "resource"=rep("Recon",amount),
                                             "rxn_formula"=rep(rxn,amount)))
        }
      }
    }

    if (length(produced)>0) {
      # Do not convert to ReconID if only primary reaction used!
      amount = dim(produced)[1]
      result_mets_produced = cbind("rxn_id"=rep(toString(gene_in),amount),
                                   "step"=rep(0,amount),
                                   "met_in"=rep(NA,amount),
                                   "left_right"=rep("right",amount),
                                   "met_short"=produced[,"Compound"],
                                   "met_long"=produced[,"Compound"],
                                   "hmdb"=rep(NA,amount),
                                   "kegg"=produced[,"KEGG_id"],
                                   "chebi"=produced[,"CheBI_id"],
                                   "pubchem"=produced[,"PubChem_id"],
                                   "InChI_key"=produced[,"InchI_key"],
                                   "rxn_name"=rep(NA,amount),
                                   "rxn"=rep(NA,amount),
                                   "resource"=produced[,"Data Source"])
      tmp=NULL
      for (i in 1:dim(result_mets_produced)[1]){
        source = rxns[grep(result_mets_produced[i,"resource"], names(rxns), fixed=TRUE)]
        tmp=c(tmp,
              paste(source[grep(result_mets_produced[i,"met_short"], source, fixed = TRUE)],collapse = ";"))
      }
      result_mets_produced = cbind(result_mets_produced,"rxn_formula"=tmp)


      if (step>0) {
        if (reconVersion==2.2){

          index_produced = getMets2Recon2.2ID(produced, recon2chebi)

          amount = length(convert(model@met_id[index_produced]))

          result_mets_produced = rbind(result_mets_consumed,
                                       cbind("rxn_id"=rep(toString(gene_in),amount),
                                             "step"=rep(0,amount),
                                             "met_in"=rep(NA,amount),
                                             "left_right"=rep("right",amount),
                                             "met_short"=convert(model@met_id[index_produced]),
                                             "met_long"=convert(model@met_name[index_produced]),
                                             "hmdb"=rep(NA,amount),
                                             "kegg"=rep(NA,amount),
                                             "chebi"=convert(recon2chebi[index_produced,2]),
                                             "pubchem"=rep(NA,amount),
                                             "InChI_key"=convert(recon2chebi[index_produced,3]),
                                             "rxn_name"=rep(NA,amount),
                                             "rxn"=rep(NA,amount)))
          # ,
          # "resource"=rep("Recon",amount),
          # "rxn_formula"=getReactionsRecon(index_produced, model)))
        } else if (reconVersion==2.0){

          index_produced = getMets2ReconID(produced, model)

          # replace empty entries by NA
          for (j in 1:length(index_produced)){
            if (length(convert(model$metHMDB[index_produced[j]]))==0) model$metHMDB[index_produced[j]][[1]]=list(matrix(NA))
          }

          amount = length(convert(model$mets[index_produced]))
          rxn=NA
          if (length(index_produced)>0) rxn = getReactionsRecon(index_produced, model)
          result_mets_produced = rbind(result_mets_consumed,
                                       cbind("rxn_id"=rep(toString(gene_in),amount),
                                             "step"=rep(0,amount),
                                             "met_in"=rep(NA,amount),
                                             "left_right"=rep("right",amount),
                                             "met_short"=convert(model$mets[index_produced]),
                                             "met_long"=convert(model$metNames[index_produced]),
                                             "hmdb"= convert(model$metHMDB[index_produced]),
                                             "kegg"= convert(model$metKeggID[index_produced]),
                                             "chebi"= convert(model$metCHEBIID[index_produced]),
                                             "pubchem"= convert(model$metPubChemID[index_produced]),
                                             "InChI_key"=rep(NA,amount),
                                             "rxn_name"=rep(NA,amount),
                                             "rxn"=rep(NA,amount),
                                             "resource"=rep("Recon",amount),
                                             "rxn_formula"=rep(rxn,amount)))
        }
      }
    }

    if (length(consumed)>0 && length(produced)>0) {
      result_mets = rbind(result_mets_consumed, result_mets_produced)
    } else if (length(consumed)==0 && length(produced)>0) {
      result_mets = result_mets_produced
    } else if (length(consumed)>0 && length(produced)==0) {
      result_mets = result_mets_consumed
    } else {
      result_mets = NULL
    }

    if ((step>0) & ((length(index_consumed)==0) | (length(index_produced)==0))){
      result_mets = NULL
    }

    # rvalMets = result_mets

    result_mets_consumed=NULL
    result_mets_produced=NULL

    if (!is.null(result_mets)) {

      # rval = removeMetsFromSet(result_mets, model, index_consumed, index_produced)
      rval = list("result_mets"=result_mets, "index_produced"=index_produced, "index_consumed"=index_consumed)
      
      result_mets = rval$result_mets
      index_produced = rval$index_produced
      index_consumed = rval$index_consumed

      if (step>0) {

        index_primary = unique(c(index_produced, index_consumed))
        if (length(index_primary)>0) {
          # step -1,1 ###################################################################################
          if (reconVersion==2.2){
            A = getPreviousNextRecon2.2(model, index_primary, produced=FALSE, max_rxns, recon2chebi)
            B = getPreviousNextRecon2.2(model, index_primary, produced=TRUE, max_rxns, recon2chebi)
          } else if (reconVersion==2.0){
            A = getPreviousNext(model, index_primary, produced=FALSE, max_rxns, step=1)
            B = getPreviousNext(model, index_primary, produced=TRUE, max_rxns, step=1)
          }

          # result_mets_1 = removeMetsFromSet(rbind(A$mets,B$mets),model,NULL,NULL)$result_mets
          result_mets_1 = rbind(A$mets,B$mets)
          if (is.null(result_mets_1)) result_mets_1 = result_mets  
        }

        if (step>1) {

          tmp=result_mets_1
          colnames(tmp)[8:10]=c("KEGG_id","CheBI_id","PubChem_id")
          index_1 = getMets2ReconID(tmp, model)

          if (length(index_1)>0) {
            # step -2,2 ###################################################################################
            if (reconVersion==2.2){
              A = getPreviousNextRecon2.2(model, index_1, produced=FALSE, max_rxns, recon2chebi)
              B = getPreviousNextRecon2.2(model, index_1, produced=TRUE, max_rxns, recon2chebi)
            } else if (reconVersion==2.0){
              A = getPreviousNext(model, index_1, produced=FALSE, max_rxns, step=2)
              B = getPreviousNext(model, index_1, produced=TRUE, max_rxns, step=2)
            }

            # result_mets_2 = removeMetsFromSet(rbind(A$mets,B$mets),model,NULL,NULL)$result_mets
            result_mets_2 = rbind(A$mets,B$mets)
            if (is.null(result_mets_2)) result_mets_2 = result_mets_1
          }
        }

        if (reconVersion==2.2){
          result_mets=result_mets[,-c(which(colnames(result_mets)=="hmdb"),which(colnames(result_mets)=="kegg"),which(colnames(result_mets)=="pubchem"))]
        } else {
          which(colnames(result_mets)=="InChI_key")
          
          result_mets=result_mets[,-c(which(colnames(result_mets)=="InChI_key"))]
        }

        if (step==2){
          result_from_consumed_plus_produced = rbind(result_mets,result_mets_1,result_mets_2)
        } else if (step==1){
          result_from_consumed_plus_produced = rbind(result_mets,result_mets_1)
        } else if (step==0){
          result_from_consumed_plus_produced = result_mets
        }

        if (dim(result_from_consumed_plus_produced)[1]>0){
          # dir.create("./results", showWarnings = FALSE)
          # library(XLConnect)
          # genExcelFileShort(result_from_consumed_plus_produced, paste("./results/mets",HGNC,"xls", sep="."))


          # remove: "Exchange/demand reaction", "Transport, extracellular",
          index=grep("Exchange/demand reaction", result_from_consumed_plus_produced[,"rxn"], fixed=TRUE)
          if (length(index)>0) result_from_consumed_plus_produced = result_from_consumed_plus_produced[-index,]
          index=grep("Transport, extracellular", result_from_consumed_plus_produced[,"rxn"], fixed=TRUE)
          if (length(index)>0) result_from_consumed_plus_produced = result_from_consumed_plus_produced[-index,]

          # genExcelFileShort(as.data.frame(result_from_consumed_plus_produced), paste("./results/mets.min.", HGNC,".xls",sep=""))

          # A = result_from_consumed_plus_produced[,"met_short"]
          A = result_from_consumed_plus_produced[,"met_long"]
          unic = !duplicated(A)  ## logical vector of unique values
          index = seq_along(A)[unic]  ## indices
          result_from_consumed_plus_produced = result_from_consumed_plus_produced[index,]

          # genExcelFileShort(as.data.frame(result_from_consumed_plus_produced), paste("./results/mets.min.unic.", HGNC,".xls",sep=""))

          rvalMets=result_from_consumed_plus_produced

        }
      } else {

        # A = result_mets[,"met_short"]
        A = result_mets[,"met_long"]
        unic = !duplicated(A)  ## logical vector of unique values
        index = seq_along(A)[unic]  ## indices
        result_mets = result_mets[index,]

        rvalMets=result_mets
      }
    }

    if (reconVersion==2.2){
      getInChI <- function(InChI_key){
        # InChI_key="OVRNDRQMDRJTHS-KEWYIRBNSA-N"
        require("RCurl")
        require("XML")

        # url=sprintf("https://www.ebi.ac.uk/chembl/api/data/substructure/%s",InChI_key)
        url=sprintf("https://www.chemspider.com/InChI.asmx/InChIKeyToInChI?inchi_key=%s",InChI_key)

        myXML=getURL(url, ssl.verifyhost=FALSE, ssl.verifypeer=FALSE)
        data=xmlParse(myXML)
        xml_data=xmlToList(data)

        # message(unlist(strsplit(xml_data, "=", fixed = TRUE))[2])

        # return(unlist(strsplit(xml_data$molecules$molecule$molecule_structures$standard_inchi, "=", fixed = TRUE))[2])
        return(xml_data)
      }


      for (i in 1:dim(rvalMets)[1]){
        if (!is.na(rvalMets[i,8])){
          if (nchar(rvalMets[i,8])==27){
            tmp=NULL
            try({tmp=getInChI(as.vector(rvalMets[i,8]))}, silent = TRUE)

            if (!is.null(tmp)){
              # message(as.vector(rvalMets[i,8]))
              rvalMets[i,8]=getInChI(as.vector(rvalMets[i,8]))
            }
          }
        }
      }
    }

    return(list("mets"=rvalMets))
  }
}
