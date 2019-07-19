removeMetsFromSet <- function(index, model){
#index =  index_primary
  
  # primary:
  # "hydrogencarbonate" 
  # step 1:
  # "Bicarbonate","potassium","Sodium","superoxide","Chloride","sulfate","Nicotinamide"   
  
  # mets2remove = c("Water","water",
  #     "O2","dioxygen","Hydrogen","hydron","H+","proton","Proton",
  #     "ATP","ADP","AMP","dAMP","dADP","dATP",
  #     "GTP","GDP","GMP","dGMP","dGDP","dGTP",
  #     "CTP","CDP","CMP","dCMP","dCDP","dCTP",
  #     "TTP","TDP","TMP","dTMP","dTDP","dTTP",
  #     "UTP","UDP","UMP","dUMP","dUDP","dUTP",
  #     "carbon dioxide",
  #     "hydrogenphosphate","phosphate(3-)","Phosphate","Orthophosphate",
  #     "Diphosphate",
  #     "hydrosulfide",
  #     "Nicotinamide adenine dinucleotide phosphate",
  #     "Nicotinamide adenine dinucleotide phosphate - reduced",
  #     "Nicotinamide adenine dinucleotide",
  #     "Nicotinamide adenine dinucleotide - reduced",
  #     "NADH","NAD(1-)","NAD(+)","NADP(+)","NADP(3-)","NADPH(4-)","NADPH", 
  #     "Flavin adenine dinucleotide oxidized",
  #     "Flavin adenine dinucleotide reduced",
  #     "FAD","FADH2","acceptor","hydrogen donor","flavoprotein",
  #     "Coenzyme A")
  
  #candidates:
#   "Sodium",
#   "Chloride",
#   "Hydrogen peroxide",
#   "potassium",
#   "calcium(2+)"
  
  # ind=grep("Hydrogen peroxide", as.vector(unlist(model$metNames)), fixed=TRUE)
  # as.vector(unlist(model$metNames))[ind]
  # ind
  # 
  # ind=23
  # model$metNames[ind]
  # model$metCHEBIID[ind]
  # model$metKeggID[ind]
  # model$metPubChemID[ind]
  # model$metInchiString[ind]
  
  # HMDB01967 carbon dioxide => P15 SLC16A1 
  
  mets2remove = rbind(c("Water","15377","C00001","962","InChI=1S/H2O/h1H2"),
                      c("O2","15379","C00007","977","InChI=1S/O2/c1-2"),
                      c("proton","15378","C00282","1038","InChI=1S/p+1/i/hH"),
                      c("ATP","15422","C00002","5957","InChI=1S/C10H16N5O13P3/c11-8-5-9(13-2-12-8)15(3-14-5)10-7(17)6(16)4(26-10)1-25-30(21,22)28-31(23,24)27-29(18,19)20/h2-4,6-7,10,16-17H,1H2,(H,21,22)(H,23,24)(H2,11,12,13)(H2,18,19,20)/p-4/t4-,6-,7-,10-/m1/s1"),
                      c("dATP","16284","C00131","15993","InChI=1S/C10H16N5O12P3/c11-9-8-10(13-3-12-9)15(4-14-8)7-1-5(16)6(25-7)2-24-29(20,21)27-30(22,23)26-28(17,18)19/h3-7,16H,1-2H2,(H,20,21)(H,22,23)(H2,11,12,13)(H2,17,18,19)/p-4/t5-,6+,7+/m0/s1"),
                      c("ADP","16761","C00008","6022","InChI=1S/C10H15N5O10P2/c11-8-5-9(13-2-12-8)15(3-14-5)10-7(17)6(16)4(24-10)1-23-27(21,22)25-26(18,19)20/h2-4,6-7,10,16-17H,1H2,(H,21,22)(H2,11,12,13)(H2,18,19,20)/t4-,6-,7-,10-/m1/s1"),
                      c("dADP","16174","C00206","188966","InChI=1S/C10H15N5O9P2/c11-9-8-10(13-3-12-9)15(4-14-8)7-1-5(16)6(23-7)2-22-26(20,21)24-25(17,18)19/h3-7,16H,1-2H2,(H,20,21)(H2,11,12,13)(H2,17,18,19)/p-3/t5-,6+,7+/m0/s1"),
                      c("AMP","16027","C00020","6083","InChI=1S/C10H14N5O7P/c11-8-5-9(13-2-12-8)15(3-14-5)10-7(17)6(16)4(22-10)1-21-23(18,19)20/h2-4,6-7,10,16-17H,1H2,(H2,11,12,13)(H2,18,19,20)/p-2/t4-,6-,7-,10-/m1/s1"),
                      c("dAMP(2-)","17713","C00360","12599","InChI=1S/C10H14N5O6P/c11-9-8-10(13-3-12-9)15(4-14-8)7-1-5(16)6(21-7)2-20-22(17,18)19/h3-7,16H,1-2H2,(H2,11,12,13)(H2,17,18,19)/p-2/t5-,6+,7+/m0/s1"),
                      c("carbon dioxide","16526","C00011","280","InChI=1S/CO2/c2-1-3"),
                      c("hydrogenphosphate","43474","C00009","1004","InChI=1S/H3O4P/c1-5(2,3)4/h(H3,1,2,3,4)/p-2"),
                      c("Diphosphate","29888","C00013","1023","InChI=1S/H4O7P2/c1-8(2,3)7-9(4,5)6/h(H2,1,2,3)(H2,4,5,6)/p-3"),
                      c("Nicotinamide adenine dinucleotide phosphate - reduced","57783","C00005","22833512","InChI=1S/C21H30N7O17P3/c22-17-12-19(25-7-24-17)28(8-26-12)21-16(44-46(33,34)35)14(30)11(43-21)6-41-48(38,39)45-47(36,37)40-5-10-13(29)15(31)20(42-10)27-3-1-2-9(4-27)18(23)32/h1,3-4,7-8,10-11,13-16,20-21,29-31H,2,5-6H2,(H2,23,32)(H,36,37)(H,38,39)(H2,22,24,25)(H2,33,34,35)/p-4/t10-,11-,13-,14-,15-,16-,20-,21-/m1/s1"),
                      c("Nicotinamide adenine dinucleotide phosphate","18009","C00006","5886","InChI=1S/C21H28N7O17P3/c22-17-12-19(25-7-24-17)28(8-26-12)21-16(44-46(33,34)35)14(30)11(43-21)6-41-48(38,39)45-47(36,37)40-5-10-13(29)15(31)20(42-10)27-3-1-2-9(4-27)18(23)32/h1-4,7-8,10-11,13-16,20-21,29-31H,5-6H2,(H7-,22,23,24,25,32,33,34,35,36,37,38,39)/p-3/t10-,11-,13-,14-,15-,16-,20-,21-/m1/s1"),
                      c("Nicotinamide adenine dinucleotide - reduced","16908","C00004","439153","InChI=1S/C21H29N7O14P2/c22-17-12-19(25-7-24-17)28(8-26-12)21-16(32)14(30)11(41-21)6-39-44(36,37)42-43(34,35)38-5-10-13(29)15(31)20(40-10)27-3-1-2-9(4-27)18(23)33/h1,3-4,7-8,10-11,13-16,20-21,29-32H,2,5-6H2,(H2,23,33)(H,34,35)(H,36,37)(H2,22,24,25)/p-2/t10-,11-,13-,14-,15-,16-,20-,21-/m1/s1"),
                      c("Nicotinamide adenine dinucleotide","15846","C00003","5893","InChI=1S/C21H27N7O14P2/c22-17-12-19(25-7-24-17)28(8-26-12)21-16(32)14(30)11(41-21)6-39-44(36,37)42-43(34,35)38-5-10-13(29)15(31)20(40-10)27-3-1-2-9(4-27)18(23)33/h1-4,7-8,10-11,13-16,20-21,29-32H,5-6H2,(H5-,22,23,24,25,33,34,35,36,37)/p-1/t10-,11-,13-,14-,15-,16-,20-,21-/m1/s1"),
                      c("Flavin adenine dinucleotide oxidized","16238","C00016","643975","InChI=1S/C27H33N9O15P2/c1-10-3-12-13(4-11(10)2)35(24-18(32-12)25(42)34-27(43)33-24)5-14(37)19(39)15(38)6-48-52(44,45)51-53(46,47)49-7-16-20(40)21(41)26(50-16)36-9-31-17-22(28)29-8-30-23(17)36/h3-4,8-9,14-16,19-21,26,37-41H,5-7H2,1-2H3,(H,44,45)(H,46,47)(H2,28,29,30)(H,34,42,43)/t14-,15+,16+,19-,20+,21+,26+/m0/s1"),
                      c("Flavin adenine dinucleotide reduced","17877","C01352","446013","InChI=1S/C27H35N9O15P2/c1-10-3-12-13(4-11(10)2)35(24-18(32-12)25(42)34-27(43)33-24)5-14(37)19(39)15(38)6-48-52(44,45)51-53(46,47)49-7-16-20(40)21(41)26(50-16)36-9-31-17-22(28)29-8-30-23(17)36/h3-4,8-9,14-16,19-21,26,32,37-41H,5-7H2,1-2H3,(H,44,45)(H,46,47)(H2,28,29,30)(H2,33,34,42,43)/p-2/t14-,15+,16+,19-,20+,21+,26+/m0/s1"),
                      c("Coenzyme A","15346","C00010","6816","InChI=1S/C21H36N7O16P3S/c1-21(2,16(31)19(32)24-4-3-12(29)23-5-6-48)8-41-47(38,39)44-46(36,37)40-7-11-15(43-45(33,34)35)14(30)20(42-11)28-10-27-13-17(22)25-9-26-18(13)28/h9-11,14-16,20,30-31,48H,3-8H2,1-2H3,(H,23,29)(H,24,32)(H,36,37)(H,38,39)(H2,22,25,26)(H2,33,34,35)/p-4/t11-,14-,15-,16+,20-/m1/s1"),
                      c("Hydrogen peroxide","16240","C00027","784","InChI=1S/H2O2/c1-2/h1-2H"))
  
  colnames(mets2remove)=c("Name","chebi","kegg","pubchem","inchi")
  
  
  # ind = which(as.vector(unlist(model$metCHEBIID[index]))%in%mets2remove[,"chebi"])
  ind <- which(as.vector(unlist(lapply(model$metCHEBIID[index], paste) %in%mets2remove[,"chebi"])))
  # as.vector(unlist(model$metNames[index[ind]]))
  if (length(ind)>0) index = index[-ind]
  
  # ind = which(as.vector(unlist(model$metKeggID[index]))%in%mets2remove[,"kegg"])
  ind <- which(as.vector(unlist(lapply(model$metKEGGID[index], paste) %in%mets2remove[,"kegg"])))
  if (length(ind)>0) index = index[-ind]
  
  # ind = which(as.vector(unlist(model$metPubChemID[index]))%in%mets2remove[,"pubchem"])
  ind <- which(as.vector(unlist(lapply(model$metPubChemID[index], paste) %in%mets2remove[,"pubchem"])))
  if (length(ind)>0) index = index[-ind]

  # ind = which(as.vector(unlist(model$metInchiString[index]))%in%mets2remove[,"inchi"])
  ind <- which(as.vector(unlist(lapply(model$metInChIString[index], paste) %in%mets2remove[,"inchi"])))
  if (length(ind)>0) index = index[-ind]

  return(index)
  
}