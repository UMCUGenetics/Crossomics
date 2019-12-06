library(data.table)
library(stringr)
library(rstudioapi)

code_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)

date <- "2019-12-04"
seed_date <- "2019-10-22"


                         
# seeds <- sub(x = list.files(path = paste0(code_dir, "/../Results/Mock_genes/",seed_date,"/")),
#              pattern = ".*seed([^seed]*?)\\.txt",
#              replacement = "\\1")

seed_files <- list.files(path = paste0(code_dir, "/../Results/Mock_genes/",seed_date,"/"))
seed_files <- seed_files[grep("mock_genes", seed_files)]
seeds <- sub(x = seed_files, pattern = ".*seed([^seed]*?)\\.txt", replacement = "\\1")


patient_info_file <- "Crossomics_DBS_Marten_trimmed20191205.RData"

load(paste0(code_dir,"/../Data/", patient_info_file))

# if("Patient.number.in.set" %in% colnames(xls_data)){
#   xls_data$Patient <- unlist(lapply(xls_data$`Patient.number.in.set`, function(x) unlist(strsplit(x, split = "\\."))[1]))
# } else if ("Patient.number" %in% colnames(xls_data)){
#   xls_data$Patient <- unlist(lapply(xls_data$Patient.number, function(x) unlist(strsplit(x, split = "\\."))[1]))
# } else {
#   stop("column names of xls_data do not follow a known pattern")
# }

xls_data_backup <- xls_data
xls_data <- as.data.table(xls_data_backup)

# rbenchmark::benchmark()
# xls_data[, "Patient" := tstrsplit(Patient, "[P.]")[2]]
# xls_data[, "Patient" := paste0("P",str_pad(Patient, width = 3, pad = "0", side = "left"))]
xls_data[, Gene := str_replace_all(Gene, " ", "")]
# xls_data[, PatientID := do.call(paste, c(.SD, sep = ":")), .SDcols= c("Patient", "Dataset")]
# xls_data[, "DBS" := .N, by = PatientID]


# create a new column `x` with the three columns collapsed together
xls_unique <- xls_data[!duplicated(apply( xls_data[ , c( 'Dataset' , 'Patient' ) ] , 1 , paste , collapse = "_" )),]

# Remove patients/disease genes who are artifacts
# xls_unique <- xls_unique[Gene != "NA" & (is.na(xls_unique[,`Judith.`]) | `Judith.` == "Inclusion"),]
xls_unique <- as.data.table(xls_unique)


seed_list <- list()
dir.create("missing_seeds")

# for(i in 1:nrow(xls_unique)){
for(i in 1:nrow(xls_unique)){
  incomplete_seeds <- NULL
  # patient <- xls_unique$Patient[i]
  # dataset <- xls_unique$Dataset[i]
  # patient_folder <- paste0(date,"/", patient, "_",dataset)
  patient_folder <- paste0(date, "/", xls_unique$file_patientIDs[i])
  path <- paste0(code_dir,"/../Results/",patient_folder)
  for(seed in seeds){
    path2 <- paste0(path,"/seed",seed,"/MSEA_results.RData")
    if(file.exists(path2)){
      next
    } else {
      incomplete_seeds <- c(incomplete_seeds, seed)
      seed_list[[i]] <- incomplete_seeds
      # write.table(as.numeric(incomplete_seeds), file = paste0("missed_seeds/patient_",i,"_seeds.txt"), row.names = FALSE, col.names = FALSE)
    }
  }
}
