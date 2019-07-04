###########################################################################
# SessionInfo -------------------------------------------------------------
###########################################################################

# Script to remove all trivial metabolites from the metabolite sets 


# OS  macOS Mojave 10.14.5

# R: 3.6.0


###########################################################################
# Script ------------------------------------------------------------------
###########################################################################

setwd("/Users/mkerkho7/DIMS2_repo/Crossomics/Results/")

metsToRemove <- c("water",
                  "O2",
                  "proton",
                  "H+",
                  "ATP",
                  "ADP",
                  "dATP",
                  "dADP",
                  "AMP",
                  "dAMP(2-)",
                  "carbon dioxide",
                  "hydrogenphosphate",
                  "Diphosphate",
                  "Nicotinamide adenine dinucleotide phosphate - reduced",
                  "Nicotinamide adenine dinucleotide phosphate",
                  "Nicotinamide adenine dinucleotide - reduced",
                  "Nicotinamide adenine dinucleotide",
                  "Flavin adenine dinucleotide oxidized",
                  "Flavin adenine dinucleotide reduced",
                  "Coenzyme A",
                  "Hydrogen peroxide",
                  "GTP",
                  "GDP",
                  "Na+",
                  "hydron")
metsToRemove <- c("Water",
                  "Nicotinamide Adenine Dinucleotide Phosphate - Reduced",
                  "Nicotinamide Adenine Dinucleotide - Reduced",
                  "Nicotinamide Adenine Dinucleotide Phosphate",
                  "Nicotinamide Adenine Dinucleotide",
                  "Flavin Adenine Dinucleotide Oxidized",
                  "Flavin Adenine Dinucleotide Reduced",
                  "Oxygen",
                  "Proton",
                  "Hydrogen Peroxide",
                  "Guanosine-5'-Diphosphate",
                  "Guanosine-5'-Triphosphate",
                  "Guanosine-5'-Monophosphate",
                  # "Guanosine",
                  # "chloride",
                  # "chloride salt",
                  "Adenosine Triphosphate",
                  "Adenosine Diphosphate",
                  "Adenosine Monophosphate",
                  "Carbon Dioxide",
                  "Coenzyme A",
                  # "Diphosphate",
                  # "Glycine"
                  )


files <- NULL
# All file names inside mss_0, mss_1 etc should be exactly the same
files <- list.files("mss_0")
for(directory in c("mss_0","mss_1","mss_2","mss_3","mss_4")){
  
}


# Function to load R data and be able to assign it to a variable
loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

for(j in names(files)){
  for(i in files[j]){
    name <- paste(j,i, sep = "/")
    # load(name)
    print(name)
    metaboliteset <- loadRData(name)
    # metaboliteset <- load(paste(getwd(), name, sep = "/"))
    # load("/Users/mkerkho7/DIMS2_repo/Crossomics/Results/mss_4/AACS.RData")
    print(head(metaboliteset,2))
  }
}




metabolites <- as.data.frame(metaboliteset, stringsAsFactors = FALSE)[,c("met_short","met_long")]


index <- NULL
for (i in 1:nrow(metabolites)){
  if (tolower(metabolites[i,1]) %in% tolower(metsToRemove) | tolower(metabolites[i,2]) %in% tolower(metsToRemove)){
    index <- c(index, i)
  }
}

