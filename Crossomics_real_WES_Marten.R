# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Session info ------------------------------------------------------------
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

cat("
    Copied from 'Crossomics_real_WES.R' in the Metab/Metabolomics/DIMS_pipeline/R_workspace_ME/Crossomics/Crossomics_SinglePatients/src folder
on 22/03/2019 to investigate and annotate the functionality of the crossomics pipeline.

package versions:
bioDist
Cairo
XLConnect
BridgeDbR 
    ")



load("./input_DIMS/BSP20160902_2,28,34-47.RData")
# The data loaded with this looks like the output of the DIMS pipeline, in which HMDB metabolite scores and Z
# scores are shown for different patients and controls

library("bioDist")
library("Cairo") 
library("XLConnect")
library("BridgeDbR")

source("./src/Crossomics/sourceDir.R")
sourceDir("./src/Crossomics")
# Load in a bunch of supportive scripts

# SPIV
patients=c(28,41:45) #patients=c(28,41:45) #patients=c(2,28:47)
controls=c(1:17)

n_patients=length(patients)
n_controls=length(controls)




# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Separate modes only selected adducts ------------------------------------
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

adductsSummed=FALSE

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Get p values of positive mode -------------------------------------------
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

p.values.assi.pos = getPvalues(peaklist = outlist.pos.stats.more, 
                               n_patients, 
                               n_controls, 
                               assi.lab = "HMDB_code",
                               patients,
                               controls
                               ) # assi.hmdb, HMDB_code For new identification ! <- no idea what this means...



# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Filter adducts ----------------------------------------------------------
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Remove adducts that do not contain 1 or 2 as extention, rename ones that do

names = lapply(rownames(p.values.assi.pos), function(x, adducts=c(1,2), adducts_long=c("[M+Na]+","[M+K]+")) {
  # x=rownames(p.values.assi.pos)[1]
  compounds=unlist(strsplit(as.vector(x),";"))
  # compounds=compounds[-1]  # since DIMS 2.1
  index=c(1:length(compounds))
  adducts_index = grep("_", compounds, fixed=TRUE)
  
  if (length(adducts_index)>0) index=index[-adducts_index]
  
  all_adducts=compounds[adducts_index]
  keep=NULL
  for (i in 1:length(all_adducts)){
    # unlist(strsplit(as.vector(all_adducts[i]),"_"))[2] %in% adducts
    keep = c(keep, unlist(strsplit(as.vector(all_adducts[i]),"_"))[2] %in% adducts)
  }
  all_adducts=all_adducts[keep]
  
  for (i in 1:length(adducts)){
    all_adducts=gsub(paste("_",toString(adducts[i]),sep=""), paste(" ", adducts_long[i], sep=""), all_adducts)
  }
  
  compounds=compounds[index]
  if (length(all_adducts)!=0){
    if (length(compounds)!=0){
      compounds=paste(compounds, all_adducts, sep=";", collapse=";")
    } else {
      compounds=paste(all_adducts, collapse=";")
    }  
  } else if (length(compounds)>1) {
    compounds=paste(compounds, collapse=";")
  }
  
  if (length(compounds)==0) compounds="" 
  
  return(compounds)
})

index=which(names=="")
if (length(index)>0) names=names[-index]
p.values.assi.pos.filt = p.values.assi.pos[-index,] 
rownames(p.values.assi.pos.filt) = names

save(p.values.assi.pos.filt, file="p.values.assi.pos.filt.RDate")



# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Get p values of negative mode -------------------------------------------
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

p.values.assi.neg = getPvalues(outlist.neg.stats.more, 
                               n_patients, 
                               n_controls, 
                               assi.lab="HMDB_code",
                               patients,
                               controls)  




# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Filter adducts ----------------------------------------------------------
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Remove adducts that do not contain 1 or 2 as extention, rename ones that do

names = lapply(rownames(p.values.assi.neg), function(x, adducts=c(1), adducts_long=c("[M+Cl]-")) {
  # x=rownames(p.values.assi.pos)[1]
  compounds=unlist(strsplit(as.vector(x),";"))
  # compounds=compounds[-1]   # since DIMS 2.1
  index=c(1:length(compounds))
  adducts_index = grep("_", compounds, fixed=TRUE)
  
  if (length(adducts_index)>0) index=index[-adducts_index]
  
  all_adducts=compounds[adducts_index]
  keep=NULL
  for (i in 1:length(all_adducts)){
    # unlist(strsplit(as.vector(all_adducts[i]),"_"))[2] %in% adducts
    keep = c(keep, unlist(strsplit(as.vector(all_adducts[i]),"_"))[2] %in% adducts)
  }
  all_adducts=all_adducts[keep]
  
  for (i in 1:length(adducts)){
    all_adducts=gsub(paste("_",toString(adducts[i]),sep=""), paste(" ", adducts_long[i], sep=""), all_adducts)
  }
  
  compounds=compounds[index]
  if (length(all_adducts)!=0){
    if (length(compounds)!=0){
      compounds=paste(compounds, all_adducts, sep=";", collapse=";")
    } else {
      compounds=paste(all_adducts, collapse=";")
    }  
  } else if (length(compounds)>1) {
    compounds=paste(compounds, collapse=";")
  }
  
  if (length(compounds)==0) compounds="" 
  
  return(compounds)
})

index=which(names=="")
if (length(index)>0) names=names[-index]
p.values.assi.neg.filt = p.values.assi.neg[-index,] 
rownames(p.values.assi.neg.filt) = names

save(p.values.assi.neg.filt, file="p.values.assi.neg.filt.RDate")

# Positive and negative mode together
rownames(p.values.assi.pos.filt) = paste(rownames(p.values.assi.pos.filt), "pos", sep="_") 
rownames(p.values.assi.neg.filt) = paste(rownames(p.values.assi.neg.filt), "neg", sep="_") 
p.values.assi.all=rbind(p.values.assi.pos.filt, p.values.assi.neg.filt)




# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Sla alles op ------------------------------------------------------------
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

save(p.values.assi.all, file = "./results/BSP20160902_2,28,34-47.RData")

dir.create("./results", showWarnings = FALSE)
path="./results/crossomics"
dir.create(path, showWarnings = FALSE)
thresh_F_pos=2
thresh_F_neg=-2
top=20
# id="InChI_key"
id="hmdb"

# ## Generate SinglePatients gene lists ###############################################################################
# # patients=c(5,28,41:45)
# for (i in 2:length(patients)){
#   getSinglePatientGeneList(patients[i])  
# }
# #####################################################################################################################

# Then performe MSEA ######################################################################################

path2 = "./results/mss_2"
overview = NULL

Sys.time()

for (i in 1:length(patients)){

  ######################## MSEA Revon2 neighbourhood ######################################################
  dir.create(paste(path,"/P",patients[i], sep=""), showWarnings = FALSE)
  dir.create(paste(path,"/P",patients[i],"/Recon2", sep=""), showWarnings = FALSE)

  # message(paste("Patient", patients[i]))
  
  # subset to 1 patient!
  p.values=p.values.assi.all[,c(grep(paste("P",patients[i],sep=""), colnames(p.values.assi.all), fixed=TRUE),
                                                  grep("C", colnames(p.values.assi.all), fixed=TRUE))] 

  gene_list = read.table(paste("./db/P",patients[i],"_HGNC.txt",sep=""), header = FALSE, sep="\t")
  mss = paste(as.vector(unlist(gene_list)),"RData",sep=".")

  metSetResult = NULL
  nMets = NULL
  
#   for (j in 1:dim(gene_map_table)[1]){
  for (j in 1:length(mss)){
    
    if (file.exists(paste(path2, mss[j], sep="/"))){
      load(paste(path2, mss[j], sep="/"))
      gene_in=strsplit(mss[j], split = "." , fixed=TRUE)[[1]][1]
  
      # A lot of missing HMDB identifiers!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      # metaboliteSet = retVal$mets
      if (path2 == "./results/mss_2"){
        metaboliteSet = result_mets_2 # Geen idee wat ie hiermee doet. result_mets_2 is niet eerder genoemd
      } else {
        metaboliteSet = result_mets_1
      }
      
      if (!is.null(metaboliteSet)){
        if (is.null(dim(metaboliteSet))){
          metaboliteSet=data.frame(t(metaboliteSet))
        }
      }  
      
      ##################################################
      # temporarily work around to be fixed in findMetabolicEnvironment
      metaboliteSet = as.matrix(metaboliteSet)
      index = which(is.na(metaboliteSet[,"hmdb"]))
      if (length(index)>0) metaboliteSet[index,"hmdb"] = "character(0)"
      ##################################################
      
      ########### Recon 2.0 ############################
      index = which(metaboliteSet[,"hmdb"] == "character(0)") 
      
      # present in BridgeDB?
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
  
      nMets=c(nMets, dim(metaboliteSet)[1])
  
      message(paste("gene: ", gene_in))
      message(paste("metaboliteSet: ", dim(metaboliteSet)[1]))
  
      retVal = performMSEAenv(metaboliteSet, p.values, patients[i], gene_in, n_patients, thresh_F_pos, thresh_F_neg, path,  top, id, adductsSummed)
      
      p_value = as.numeric(retVal$p.value)
      if (length(p_value)==0){
        p_value=NA
      }
  
      metSetResult = rbind(metSetResult, c("p.value"=p_value, "patient"=paste("P", patients[i], sep=""), "metabolite.set"=gene_in))
    }
  }
  
  tmp=data.frame("HGNC"=metSetResult[,3],"p.value"=as.numeric(metSetResult[,"p.value"]), "metabolites"=nMets)
  genExcelFileShort(tmp[order(tmp[,"p.value"]),], paste(path,"/P",patients[i],"/Recon2/MSEA_results.xls",sep=""))

  # tmp2 = rbind(rep("", 4019),rep("", 4019),rep("", 4019))
  # tmp1 = t(tmp[order(tmp[,"p.value"]),])
  # tmp2[1, 1:dim(tmp1)[2]] = tmp1[1,] 
  # tmp2[2, 1:dim(tmp1)[2]] = tmp1[2,]
  # tmp2[3, 1:dim(tmp1)[2]] = tmp1[3,]
  # 
  # overview = rbind(overview, tmp2)
  ##########################################################################################################
}
genExcelFileShort(overview, paste(path,"/MSEA_overview.xls",sep=""))
