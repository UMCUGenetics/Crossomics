library(rstudioapi)

code_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
path <- paste0(code_dir,"/../Results/")

for (seed in (c(2341, 6734892, 83, 698, 991))){
  set.seed(seed = seed)

  # Load dataset containing patients with their disease genes (dataset is called xls_data)
  load(paste0(code_dir,"/../Data/Crossomics_DBS_Marten_Training.RData"))
  xls_data[,"PatientID"] <- sapply(strsplit(paste(xls_data$Dataset, xls_data$Patient.number, sep = "^"), split = "\\."), `[`, 1)
  uni_dat_pat <- xls_data[!duplicated(xls_data$PatientID),c("Gene","PatientID")]
  uni_genes <- unique(xls_data$Gene)
  patients <- NULL
  
  for(i in uni_genes){
    patient <- sample(uni_dat_pat[uni_dat_pat$Gene == i,"PatientID"], 1)
    patients <- c(patients,patient)
  }
  
  write.table(patients, file = paste0(path, "/patient_subset_seed",seed,".txt"), row.names = FALSE, col.names = FALSE)
}