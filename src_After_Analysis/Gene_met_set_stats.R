# Little script to save the status of all gene-metabolite sets

library(data.table)
library(rstudioapi)


setwd("/Users/mkerkho7/DIMS2_repo/Crossomics/Data/2019-08-12/")
code_dir <- paste0(dirname(rstudioapi::getActiveDocumentContext()$path),"/../Results/")

maxrxns <- c(8, 10, 12, 15, 17, 19)
steps <- c(0, 1, 2, 3, 4, 5)

date <- "2019-12-10"

gn_fi_list <- list.files("../mss/")
gn_fi_list <- sub(".RData",".RDS", gn_fi_list)



nr_dt_rows <- length(maxrxns)*length(steps)*length(gn_fi_list)

char_vec <- rep("NA", nr_dt_rows)
int_vec <- rep(0L, nr_dt_rows)

DT_met_sets <- data.table(Gene = char_vec,
                 Nr_mets = int_vec,
                 Step = int_vec,
                 Maxrxn = int_vec)

i <- 0

for(mx in maxrxns){
  for(st in steps){
    subdir <- paste0("maxrxn",mx,"/mss_",st,"_HMDBtranslated/")
    for(gn in gn_fi_list){
      i <- i + 1
      met_set <- tryCatch(readRDS(paste0(subdir,gn)), error = function(e) paste("file",gn,"not found for Step",st,"Maxrxn",mx,"\n"))
      hmdb_codes <- met_set[grep("HMDB", met_set[,"hmdb"]),"hmdb"]
      nr_unq_mets <- uniqueN(hmdb_codes)
      DT_met_sets[i, c("Gene","Nr_mets","Step","Maxrxn") := list(gn, nr_unq_mets, st, mx)]
    }
  }
}
DT_met_sets[, Step := as.factor(Step)]
DT_met_sets[, Maxrxn := as.factor(Maxrxn)]
DT_met_sets[, Nr_mets_log := log2(Nr_mets + 1)]

DT_met_sets[, Empty := sum(Nr_mets) == 0, by = Gene]
empty_genes <- unique(DT_met_sets[Empty == TRUE, Gene])

DT_met_sets <- DT_met_sets[ Empty == FALSE, ]
DT_met_sets[, Gene := str_replace(Gene, pattern = ".RDS", replacement = "")]

saveRDS(DT_met_sets, file = paste0(code_dir,"/../Results/",date,"/Metabolite_set_sizes.RDS"))

