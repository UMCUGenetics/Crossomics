# Convert xlsx bioinformatics file to RDS

# Go to some directory where the xlsx file is (only 1 xlsx file should be present) and run the whole thing
setwd("/Users/mkerkho7/DIMS2_repo/Crossomics/Data/Project 2017_008 MetabolomicsDiagnosis_DBS/RES_DBS_20170420_MetabolomicsDiagnosis_RUN2/Bioinformatics_metIS/")

name_xlsx <- list.files(pattern = "*.xlsx")
tmp_xlsx <- read.xlsx(name_xlsx)
tmp_name <- unlist(strsplit(name_xlsx, split = ".xlsx")[1])
output_name <- paste0(tmp_name, ".RDS")
saveRDS(object = tmp_xlsx, file = output_name)
