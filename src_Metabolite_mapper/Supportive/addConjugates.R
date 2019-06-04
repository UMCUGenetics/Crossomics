addConjugates <- function(metaboliteSet){
  
#   metaboliteSet => HMDB id => map identified peaklist => enrichment
#   => met_long => CoA? => mass(MNeutral)
#   => Calculate mass glycine/carnitine conjugate => identification => HMDB id => add to metaboliteSet

  load("./db/HMDB4crossomics.RData")
  
  index = grep("-CoA", metaboliteSet[,"met_long"], fixed=TRUE)
  
  for (i in 1:length(index)){
    found = table[which(table[,"CompoundName"] == toString(metaboliteSet[index[i],"hmdb"])),]
    RCO_S_CoA = found$MNeutral
    RCO_S_CoA.pos = found$Mpos
    RCO_S_CoA.neg = found$MNeg
    
    # glycine conjugate
    # RCO-S-CoA + NH2CH2COOH(glycine) => RCO-NHCH2COOH(conjugate)
    
    HS_CoA = 767.115208365 # Monoisotopic Molecular Weight HMDB
    glycine = 75.032028409 # Monoisotopic Molecular Weight HMDB
    
    # HS-CoA and glycine both have a hydrogen atom too much    
    RCO_glycine = as.numeric(RCO_S_CoA) - HS.CoA + glycine
    RCO_glycine.pos = as.numeric(RCO_S_CoA.pos) - HS.CoA + glycine
    RCO_glycine.peg = as.numeric(RCO_S_CoA.neg) - HS.CoA + glycine
  
    found$CompoundName = paste(found$CompoundName, "[M+ glycine]")
    found$Composition = "Todo"
    found$MNeutral = RCO_glycine
    
    #########################################################################################  
    carnitine = 161.105193351 # Monoisotopic Molecular Weight HMDB
    
    # HS-CoA and carnitine both have a hydrogen atom too much    
    mass.carn.con = as.numeric(mass) - HS.CoA + carnitine
    
    
    # table MNeutral lijkt afgerond op 4 decimalen!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  }
}