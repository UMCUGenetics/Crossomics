# R_location <- "/hpc/local/CentOS7/dbg_mz/R_libs/3.6.0/lib64" 
# Supply the maxrxn and the directory where this script is placed
cmd_args <- commandArgs(trailingOnly = TRUE)
patient_number <- as.numeric(cmd_args[1])
thresholds <- cmd_args[2]
max_rxns <- cmd_args[3]
max_rxns <- as.numeric(unlist(strsplit(max_rxns, split = ",")))
steps <- cmd_args[4]
steps <- as.numeric(unlist(strsplit(steps, split = ",")))
code_dir <- cmd_args[5]
seed_file <- cmd_args[6]
seed <- as.integer(sub(".txt", "", sub(".*seed", "", seed_file)))
R_location <- cmd_args[7]
R_location <- sub("/bin", "/lib64", R_location)
outdir <- "Results"

thresh_df <- sapply(strsplit(unlist(strsplit(thresholds, split = ",")), ";"), `[`)

thresh_neg_list <- as.numeric(thresh_df[1,])
thresh_pos_list <- as.numeric(thresh_df[2,])

suppressMessages(library("stringr",lib.loc = R_location)) # string manipulation, add leading 0's
suppressMessages(library("Cairo",lib.loc = R_location))
suppressMessages(library("backports",lib.loc = R_location))
suppressMessages(library("crayon",lib.loc = R_location))
suppressMessages(library("vctrs",lib.loc = R_location))
suppressMessages(library("tidyselect",lib.loc = R_location))
suppressMessages(library("tidyr",lib.loc = R_location))
suppressMessages(library("data.table",lib.loc = R_location))
suppressMessages(library("dplyr",lib.loc = R_location))

message(patient_number, "\n", 
        paste(thresh_neg_list, collapse = " "), "\n", 
        paste(thresh_pos_list, collapse = " "), "\n", 
        paste(max_rxns, collapse = " "), "\n", 
        paste(steps, collapse = " "), "\n", 
        code_dir, "\n", 
        seed, "\n", 
        R_location)

test1 <- c("test","testtie")

save(test1, file = paste0(code_dir,"/../", outdir, "test"))