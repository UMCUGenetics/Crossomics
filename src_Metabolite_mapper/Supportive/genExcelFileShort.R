genExcelFileShort <- function(list, wbfile) { 
#   list=result_from_consumed_plus_produced$mets
#   wbfile=paste("./results/mets",gene_in,"xls", sep=".")
  
  filelist <- "metabolites"
  wb <- loadWorkbook(wbfile, create = TRUE)
  createSheet(wb, name = filelist)
  writeWorksheet(wb, list, sheet = filelist)
  setColumnWidth(wb, sheet = filelist, column = c(1:ncol(list)), width = -1 )
  saveWorkbook(wb)
  rm(wb)
}
