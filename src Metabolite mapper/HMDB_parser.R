
countOccurence <- function(formula, atom) {
#   formula ="Xc2CH7NO3"
#   atom="C"
  
  nr=NULL
  
  characters = as.vector(unlist(strsplit(formula, "")))
  
  split=as.vector(unlist(strsplit(formula, split=atom,fixed=TRUE)))
  start = nchar(split[1]) + nchar(atom) + 1
  
  if (grepl("[a-z]", characters[start])) { # lower case thus second letter of abbreviation
#     if (length(split)>2) { # first letter of abbriviation occurs second time 
#       start = nchar(split[1]) + nchar(split[2]) + 2*nchar(atom) + 1
#       nr = check(start)
#     } else {
      nr = 0
#     }
  } else if ((nchar(split[1]) + nchar(atom)) == length(characters)) { # last element one time (without number)
    nr = 1
  } else if (is.na(characters[start])) {
    nr = 0
  } else if (is.na(as.numeric(characters[start]))) {
    nr = 1
  } else {
    nr = as.numeric(characters[start])

#   Ignore!!!!!!!!!!!!!!!!!!!!!!!!!     
#   ##############################################################
#   Warning message:
#   NAs introduced by coercion
#   ##############################################################
    
    if ((start + 1) <= length(characters)) {
      if (!is.na(as.numeric(characters[start + 1]))){
        nr = as.numeric((paste(characters[start], characters[start + 1], sep="")))
        if ((start + 2)<=length(characters)){
          if (!is.na(as.numeric(characters[start + 2]))){
            nr = as.numeric((paste(nr, characters[start + 2], sep="")))
          }
        }
      }
    }
  }
  
  return(nr)
  #return(check(start))
}

# formula ="Xc2CH7NO3"
# atom="H"
# countOccurence(formula, atom)

library(XML)

H = 1.00782503207 
D = 2.0141017778  
T_ = 3.0160492777 
C = 12.000000 
C13 = 13.0033548378 # C13
N = 14.0030740049 
N15 = 15.0001088983 # N15
O = 15.99491461956 
P = 30.97376163  
S = 31.97207100 
Cl = 34.96885268 
I = 126.904473 
Ca = 39.96259098  
Mg = 23.985041700 
K = 38.96370668  
Na = 22.9897692809    
Cr = 49.9460442  
Co = 58.9331951  
Cu = 62.9295976     
F_ = 18.99840323 
Fe = 53.9396106  
Al = 26.98153863 
Mo = 91.906811    
Zn = 63.9291423  
Rb = 84.911789738
Mn = 54.9380452
Se = 73.9224764 
Zr = 89.9047044
Te = 119.904020
Sn = 111.904818 
Ti = 45.9526317 
W = 179.946704
As = 74.9215965
B = 10.0129370 
Ni = 57.9353430 
Br = 78.9183371
Ge = 69.9242474
V = 49.9471585
Hg = 195.965833
Cd = 105.906460
Sr =  83.913425
Sb = 120.9038157
Au = 196.9665688
Ba = 129.9063208
Ta = 179.9474648
Pb = 203.9730436
Li = 6.015122795
Ru = 95.907599
Cs = 132.905451933
Pd = 101.905609
Pt = 189.959933
Ce = 135.907172
Nd = 141.9077233
Re = 184.9529550
Tl = 202.9723442   
Hf = 141.9077233
Th = 230.0331338
Gd = 151.9197910 
Si = 27.9769265325
Bi = 208.9803987 
Ag = 106.905097 
Ga = 68.9255736 

e = 0.00054857990946

options(digits=16)
# setwd("./src")
message(getwd())

getMNeutral <- function(formula) {
  
  MNeutral = NULL
  
  nrH = countOccurence(formula, "H")
  nrD = countOccurence(formula, "D")
  nrT = countOccurence(formula, "T")
  nrC = countOccurence(formula, "C")
  nrC13 = countOccurence(formula, "Xc") # Xc is C13
  nrN = countOccurence(formula, "N")
  nrN15 = 0
  nrO = countOccurence(formula, "O")
  nrP = countOccurence(formula, "P")
  nrS = countOccurence(formula, "S")
  nrCl = countOccurence(formula, "Cl")
  nrI = countOccurence(formula, "I") 
  nrCa = countOccurence(formula, "Ca")  
  nrMg = countOccurence(formula, "Mg") 
  nrK = countOccurence(formula, "K")  
  nrNa = countOccurence(formula, "Na")    
  nrCr = countOccurence(formula, "Cr")  
  nrCo = countOccurence(formula, "Co")  
  nrCu = countOccurence(formula, "Cu")     
  nrF = countOccurence(formula, "F") 
  nrFe = countOccurence(formula, "Fe")  
  nrAl = countOccurence(formula, "Al") 
  nrMo = countOccurence(formula, "Mo")    
  nrZn = countOccurence(formula, "Zn")  
  nrRb = countOccurence(formula, "Rb")
  nrMn = countOccurence(formula, "Mn")
  nrSe = countOccurence(formula, "Se") 
  nrZr = countOccurence(formula, "Zr")
  nrTe = countOccurence(formula, "Te")
  nrSn = countOccurence(formula, "Sn") 
  nrTi = countOccurence(formula, "Ti") 
  nrW = countOccurence(formula, "W")
  nrAs = countOccurence(formula, "As")
  nrB = countOccurence(formula, "B") 
  nrNi = countOccurence(formula, "Ni") 
  nrBr = countOccurence(formula, "Br")
  nrGe = countOccurence(formula, "Ge")
  nrV = countOccurence(formula, "V")
  nrHg = countOccurence(formula, "Hg")
  nrCd = countOccurence(formula, "Cd")
  nrSr = countOccurence(formula, "Sr")
  nrSb = countOccurence(formula, "Sb")
  nrAu = countOccurence(formula, "Au")
  nrBa = countOccurence(formula, "Ba")
  nrTa = countOccurence(formula, "Ta")
  nrPb = countOccurence(formula, "Pb")
  nrLi = countOccurence(formula, "Li")
  nrRu = countOccurence(formula, "Ru")
  nrCs = countOccurence(formula, "Cs")
  nrPd = countOccurence(formula, "Pd")
  nrPt = countOccurence(formula, "Pt")
  nrCe = countOccurence(formula, "Ce")
  nrNd = countOccurence(formula, "Nd")
  nrRe = countOccurence(formula, "Re")
  nrTl = countOccurence(formula, "Tl")   
  nrHf = countOccurence(formula, "Hf")
  nrTh = countOccurence(formula, "Th")
  nGd = countOccurence(formula, "Gd") 
  nSi = countOccurence(formula, "Si")
  nBi = countOccurence(formula, "Bi")
  nAg = countOccurence(formula, "Ag")
  nGa = countOccurence(formula, "Ga")
  
  MNeutral = nrH*H + nrD*D + nrT*T_ + nrC*C + nrC13*C13 + nrN*N + nrN15*N15 + nrO*O + nrP*P + nrS*S + nrCl*Cl +
    I*nrI + Ca*nrCa + Mg*nrMg + K*nrK + Na*nrNa + Cr*nrCr + Co*nrCo + Cu*nrCu + F_*nrF + Fe*nrFe + Al*nrAl + Mo*nrMo +
    Zn*nrZn + Rb*nrRb + Mn*nrMn + Se*nrSe + Zr*nrZr + Te*nrTe + Sn*nrSn + Ti*nrTi + W*nrW + As*nrAs + B*nrB + Ni*nrNi +
    Br*nrBr + Ge*nrGe + V*nrV + Hg*nrHg + Cd*nrCd + Sr*nrSr + Sb*nrSb + Au*nrAu + Ba*nrBa + Ta*nrTa + Pb*nrPb + Li*nrLi +
    Ru*nrRu + Cs*nrCs + Pd*nrPd + Pt*nrPt + Ce*nrCe + Nd*nrNd + Re*nrRe + Tl*nrTl + Hf*nrHf + Th*nrTh + nGd*Gd + nSi*Si +
    nBi*Bi + nAg*Ag + nGa*Ga
  
  row = c("nrH" = nrH,
          "nrD" = nrD,
          "nrT" = nrT,
          "nrC" = nrC,
          "nrC13" = nrC13,
          "nrN" = nrN,
          "nrN15" = 0,
          "nrO" = nrO,
          "nrP" = nrP,
          "nrS" = nrS,
          "nrCl" = nrCl,
          "nrI" = nrI,
          "nrCa" = nrCa,
          "nrMg" = nrMg,
          "nrK" = nrK,
          "nrNa" = nrNa,
          "nrCr" = nrCr,
          "nrCo" = nrCo,
          "nrCu" = nrCu,
          "nrF" = nrF,
          "nrFe" = nrFe,
          "nrAl" = nrAl,
          "nrMo" = nrMo,
          "nrZn" = nrZn,
          "nrRb" = nrRb,
          "nrMn" = nrMn,
          "nrSe" = nrSe,
          "nrZr" = nrZr,
          "nrTe" = nrTe,
          "nrSn" = nrSn,
          "nrTi" = nrTi,
          "nrW" = nrW,
          "nrAs" = nrAs,
          "nrB" = nrB,
          "nrNi" = nrNi,
          "nrBr" = nrBr,
          "nrGe" = nrGe,
          "nrV" = nrV,
          "nrHg" = nrHg,
          "nrCd" = nrCd,
          "nrSr" = nrSr,
          "nrSb" = nrSb,
          "nrAu" = nrAu,
          "nrBa" = nrBa,
          "nrTa" = nrTa,
          "nrPb" = nrPb,
          "nrLi" = nrLi,
          "nrRu" = nrRu,
          "nrCs" = nrCs,
          "nrPd" = nrPd,
          "nrPt" = nrPt,
          "nrCe" = nrCe,
          "nrNd" = nrNd,
          "nrRe" = nrRe,
          "nrTl" = nrTl,
          "nrHf" = nrHf,
          "nrTh" = nrTh,
          "nGd" = nGd, 
          "nSi" = nSi,
          "nBi" = nBi,
          "nAg" = nAg,
          "nGa" = nGa)
  
  
  return(list("MNeutral" = MNeutral, "row" = row))
}

getExtraMetRubenDB <- function(){
  
  gmz_list = read.csv("./db/Extra_metaboites_Ruben.txt", header = F, sep = "\t", quote = "\"")
  
  table = NULL
  
  for (i in 1:dim(gmz_list)[1]){
    
    tmp = getMNeutral(as.vector(gmz_list[i,2]))
    MNeutral = tmp$MNeutral
    MNeg = MNeutral - H + e
    
    table = rbind(table, c("CompoundName" =  as.vector(gmz_list[i,1]), 
                           "ID" = as.vector(gmz_list[i,3]),
                           "Composition" = as.vector(gmz_list[i,2]), 
                           "MNeutral" = MNeutral,
                           "MNeg" = MNeg,
                           "Mpos" = MNeutral + H - e, tmp$row))
  }
  
  return(table)
}

getGmzDB <- function(){
  
  gmz_list = read.csv("../data/DB_gmz.csv", header = TRUE, sep = "\t", quote = "\'")
  
  table = NULL
  
  for (i in 1:dim(gmz_list)[1]){
    
    tmp = getMNeutral(as.vector(gmz_list[i,"Composition"]))
    MNeutral = tmp$MNeutral
    MNeg = MNeutral - H + e
    
    table = rbind(table, c("CompoundName" =  as.vector(gmz_list[i,"CompoundName"]), 
                           "Class" = as.vector(gmz_list[i,"Class"]),
                           "Composition" = as.vector(gmz_list[i,"Composition"]), 
                           "MNeg" = MNeg,
                           "Mpos" = MNeutral + H - e, tmp$row))
    
  }
  
  return(table)
  
}

getGmzISDB <- function(){
  
  gmz_list = read.csv("../data/DB_is_gmz.csv", header = TRUE, sep = "\t", quote = "\'")
  
  table = NULL
  
  for (i in 1:dim(gmz_list)[1]){
    
    tmp = getMNeutral(as.vector(gmz_list[i,"Composition"]))
    MNeutral = tmp$MNeutral
    MNeg = MNeutral - H + e
    
    table = rbind(table, c("CompoundName" =  as.vector(gmz_list[i,"CompoundName"]), 
                           "Class" = as.vector(gmz_list[i,"Class"]),
                           "Composition" = as.vector(gmz_list[i,"Composition"]), 
                           "MNeg" = MNeg,
                           "Mpos" = MNeutral + H - e, tmp$row))
    
  }
  
  return(table)
  
}

getUrineDB <- function(){

  dir = "../data/hmdb_metabolites"
  
  #load("files.RData")
  files = list.files(dir)
  #files = c("HMDB00001.xml")
  #files = c("HMDB13711.xml")
  
  # files = c("hmdb_metabolites.xml")
  # > data = xmlParse(paste(dir, files[i], sep="/"))
  # XML declaration allowed only at the start of the document
  # Extra content at the end of the document
  # Error: 1: XML declaration allowed only at the start of the document
  # 2: Extra content at the end of the document
  
  table = NULL
  
  for (i in 1:length(files)){
    
    data = xmlParse(paste(dir, files[i], sep="/"))
  
    xml_data = xmlToList(data)
    
    row = NULL
    
    tmp = xml_data[["biofluid_locations"]]
    tmp = as.vector(unlist(tmp))

    message(i)
      
    tmp = xml_data[["predicted_properties"]]
    if (tmp=="\n  ") next
    tmp = as.vector(unlist(tmp))
      
    formula = tmp[which(tmp=="formula") + 1]
    if (length(formula)==0 | is.null(formula)) next
      
      #         formula="C1H2O0N0P0"
      
    tmp = getMNeutral(formula)      
      
    MNeg = tmp$MNeutral - H + e
    MassDenes = tmp[which(tmp=="average_mass") + 1]
    Mono_mass = tmp[which(tmp=="mono_mass") + 1]
    
    if ("Urine"%in%tmp) {
      row = c("Accession" = xml_data[["accession"]],
              "CompoundName" =  paste(xml_data[["name"]], "[Urine]", sep=" "), 
              "Average_mass" = MassDenes, 
              "Composition" = formula, 
              "MNeutral" = tmp$MNeutral,  
              "MNeg" = MNeg,
              "Mpos" = tmp$MNeutral + H - e,
              "DevMneuMavg" = as.numeric(MassDenes) - tmp$MNeutral,  
              "Mono_mass" = Mono_mass, tmp$row)
    } else {
      row = c("Accession" = xml_data[["accession"]],
              "CompoundName" =  xml_data[["name"]], 
              "Average_mass" = MassDenes, 
              "Composition" = formula, 
              "MNeutral" = tmp$MNeutral,  
              "MNeg" = MNeg,
              "Mpos" = tmp$MNeutral + H - e,
              "DevMneuMavg" = as.numeric(MassDenes) - tmp$MNeutral,  
              "Mono_mass" = Mono_mass, tmp$row)
              
    }
    
    table = rbind(table, row)      
    
  }  
  return(table)
}

getHMDB4crossomics <- function(){
  
  dir = "E:/Metabolomics/Crossomics_SinglePatients/HMDB"
  
  
  #<metabolite><chebi_id>50599</chebi_id></metabolite>
  
  #load("files.RData")
  files = list.files(dir)
  #files = c("HMDB00001.xml")
  #files = c("HMDB13711.xml")
  
  # files = c("hmdb_metabolites.xml")
  # > data = xmlParse(paste(dir, files[i], sep="/"))
  # XML declaration allowed only at the start of the document
  # Extra content at the end of the document
  # Error: 1: XML declaration allowed only at the start of the document
  # 2: Extra content at the end of the document
  
  table = NULL
  
  for (i in 1:length(files)){
    
    data = xmlParse(paste(dir, files[i], sep="/"))
    
    xml_data = xmlToList(data)
    
    row = NULL
    
#     tmp = xml_data[["biofluid_locations"]]
#     tmp = as.vector(unlist(tmp))
    
    # message(i)
    
    tmp = xml_data[["predicted_properties"]]
    if (tmp=="\n  ") next
    tmp = as.vector(unlist(tmp))
    
    formula = tmp[which(tmp=="formula") + 1]
    if (length(formula)==0 | is.null(formula)) next
    
    #         formula="C1H2O0N0P0"
    
    tmp = getMNeutral(formula)      
    
    MNeg = tmp$MNeutral - H + e
    MassDenes = tmp[which(tmp=="average_mass") + 1]
    Mono_mass = tmp[which(tmp=="mono_mass") + 1]
    
    if (is.null(xml_data[["chebi_id"]])){
      chebi=NA
    } else {
      chebi=xml_data[["chebi_id"]]
    }

    if (is.null(xml_data[["inchi"]])){
      inchi=NA
    } else {
      inchi=xml_data[["inchi"]]
    }
    
    row = c("Accession" = xml_data[["accession"]],
              "CompoundName" = xml_data[["accession"]], 
              "Average_mass" = MassDenes, 
              "Composition" = formula, 
              "MNeutral" = tmp$MNeutral,  
              "MNeg" = MNeg,
              "Mpos" = tmp$MNeutral + H - e,
              "DevMneuMavg" = as.numeric(MassDenes) - tmp$MNeutral,  
              "Mono_mass" = Mono_mass,
              "Chebi" = chebi,
              "Inchi" = inchi,tmp$row)

    table = rbind(table, row)      
    
  }  
  return(table)
}

#table=getExtraMetRubenDB()

table=getHMDB4crossomics()
rownames(table) = table[,"Accession"]
table = table[,-1]
save(table, file = "./db/HMDB4crossomics_15-11-2016.RData")
# > dim(table)
# [1] 3966   70

message("klaar")

