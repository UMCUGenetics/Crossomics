genExcelFileShort <- function(list, wbfile) { 
#   list=result_from_consumed_plus_produced$mets
#   wbfile=paste("./results/mets",gene_in,"xls", sep=".")
  
  filelist <- "metabolites"
  wb <- loadWorkbook(wbfile, create = TRUE)
  createSheet(wb, name = filelist)
  writeWorksheet(wb, list, sheet = filelist)
  saveWorkbook(wb)
  rm(wb)
}
