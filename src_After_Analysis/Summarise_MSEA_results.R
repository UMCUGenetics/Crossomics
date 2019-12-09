#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Info --------------------------------------------------------------------
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# This script uses the output of Summarise_MSEA_results.R to generate graphic representations of
# the Crossomics results

# R version:  3.6.1 (2019-07-05)
# platform:   x86_64-apple-darwin15.6.0 (64-bit)
# OS:         macOS Mojave 10.14.6
# 
# libraries:
# rstudioapi    0.10
# data.table    1.12.6



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Libraries ---------------------------------------------------------------
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

library(rstudioapi)
library(data.table)
# library(RColorBrewer)
# library(ggplot2)
# library(scales)
# library(ggforce) # zoom in specific parts of plot
# library(grid) # manual grob adjustment
# library(gtable) # manual grob adjustment



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Adjustable settings -----------------------------------------------------
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Date of data
date <- "2019-12-04"



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Pre-processing ----------------------------------------------------------
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

##### Read data -----------------------------------------------------------
code_dir <- paste0(dirname(rstudioapi::getActiveDocumentContext()$path),"/../Results/")

DT <- data.table::as.data.table(readRDS(paste0(code_dir,date,"/MSEA_DT_compiled_all.RDS")))
