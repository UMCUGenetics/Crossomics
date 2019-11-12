# Little script to check the status of all gene-metabolite sets

setwd("/Users/mkerkho7/DIMS2_repo/Crossomics/Data/2019-08-12/")

maxrxns <- c(8, 10, 12, 15, 17, 19)
steps <- c(0, 1, 2, 3, 4, 5)

gn_fi_list <- list.files("../mss/")
gn_fi_list <- sub(".RData",".RDS", gn_fi_list)


library(data.table)

nr_dt_rows <- length(maxrxns)*length(steps)*length(gn_fi_list)

char_vec <- rep("NA", nr_dt_rows)
int_vec <- rep(0L, nr_dt_rows)

DT <- data.table(Gene = char_vec,
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
      nr_unq_mets <- length(unique(met_set[,"hmdb"]))
      DT[i, c("Gene","Nr_mets","Step","Maxrxn") := list(gn, nr_unq_mets, st, mx)]
    }
  }
}

DT[, median(Nr_mets), by = .(Step, Maxrxn)]
DT[, range(Nr_mets), by = .(Step, Maxrxn)]
