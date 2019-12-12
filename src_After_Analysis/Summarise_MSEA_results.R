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
library(stringr)
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
date <- "2019-12-10"



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Pre-processing ----------------------------------------------------------
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

##### Read data -----------------------------------------------------------
code_dir <- paste0(dirname(rstudioapi::getActiveDocumentContext()$path),"/../Results/")

DT <- data.table::as.data.table(readRDS(paste0(code_dir,date,"/MSEA_DT_compiled_all.RDS")))
nr_seeds <- length(unique(DT$Seed))


##### Modify initial datatable --------------------------------------------
DT[,c("Prioritised15","Prioritised10","Prioritised05","Prioritised02") := 
     list(Position <= 15, Position <= 10, Position <= 5, Position <= 2)]



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# data table per parameter ------------------------------------------------
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

DT_per_parameter <- NULL
DT_tmp1 <- data.table()
DT_tmp2 <- data.table()
DT_tmp1[, c("Step", "Z_threshold", "Max_rxn", "Prior.frac15","Prior.frac10","Prior.frac05","Prior.frac02","Rank.frac.av.rev","Rank.frac.av","Missed","Missed.frac",
            "Max_Tot.Genes","Min_Tot.Genes","Prior.sd15", "Prior.sd10", "Prior.sd05", "Prior.sd02","Rank.frac.sd.rev","Rank.frac.sd") :=  
          DT[Include == TRUE, list(
            mean(Prioritised15), # Prior.frac15
            mean(Prioritised10), # Prior.frac10
            mean(Prioritised05), # Prior.frac05
            mean(Prioritised02), # Prior.frac02
            mean(Rev.Rank.frac), # Prior.pos.frac.av.rev
            mean(Rank.frac), # Prior.pos.frac.av
            sum(P.value==1)/nr_seeds, # Missed; if a gene is missed in 1 seed, it is missed in all --> count once per parametercombination
            mean(P.value==1), # Missed.frac
            max(Last_position), # Max_Tot.Genes
            min(Last_position), # Min_Tot.Genes
            sd(Prioritised15), # Prior.sd15
            sd(Prioritised10), # Prior.sd10
            sd(Prioritised05), # Prior.sd05
            sd(Prioritised02), # Prior.sd02
            sd(Rev.Rank.frac), # Prior.pos.frac.sd.rev
            sd(Rank.frac) # Prior.pos.frac.sd
          ), by = .(Step, Z_threshold, Max_rxn)]]
DT_tmp2[, c("Step", "Z_threshold", "Max_rxn","Av_non_missed","Sd_non_missed") := DT[P.value != 1, list(mean(Position), sd(Position)), by = .(Step, Z_threshold, Max_rxn)]]
DT_per_parameter <- merge(DT_tmp1, DT_tmp2, by = c("Step","Z_threshold", "Max_rxn"))
rm(DT_tmp1,DT_tmp2)

# DT_per_parameter[, Max_rxn:=factor(Max_rxn, levels = max_rxns)]

# return 1 or 0 depending on whether the parameter combination is the best scoring or not
DT_per_parameter[,c("best15","best10","best05","best02") := list(
  ifelse(Prior.frac15 == max(Prior.frac15), 1, 0),
  ifelse(Prior.frac10 == max(Prior.frac10), 1, 0),
  ifelse(Prior.frac05 == max(Prior.frac05), 1, 0),
  ifelse(Prior.frac02 == max(Prior.frac02), 1, 0))]

# Determine best parameters according to non-missed disease genes (frac. missed < 0.5)
DT_tmp <- DT_per_parameter[,Av_non_missed, Missed.frac]
min_values <- unique(head(sort(DT_tmp[Missed.frac < 0.5, Av_non_missed]),5))

DT_per_parameter[,"best_av_NM" := Missed.frac < 0.5 & Av_non_missed %in% min_values]
# DT_per_parameter[,"best_av_NM"] <- DT_per_parameter$Missed.frac < 0.5 & DT_per_parameter$Av_non_missed %in% min_values
DT_per_parameter[best_av_NM == TRUE, best_order_NM := rank(DT_per_parameter[best_av_NM == TRUE, Av_non_missed])]

best_ratio_in_top10 <- head(sort(DT_per_parameter[,Prior.frac10], decreasing = TRUE),5)
DT_per_parameter[,"best_top10" := Prior.frac10 %in% best_ratio_in_top10]
DT_per_parameter[best_top10 == TRUE, best_order_top10 := rank(-DT_per_parameter[best_top10 == TRUE, Prior.frac10])]

best_ratio_in_top05 <- head(sort(DT_per_parameter[,Prior.frac05], decreasing = TRUE),5)
DT_per_parameter[,"best_top05" := Prior.frac05 %in% best_ratio_in_top05]
DT_per_parameter[best_top05 == TRUE, best_order_top05 := rank(-DT_per_parameter[best_top05 == TRUE, Prior.frac05])]

DT_per_parameter <- data.table(DT_per_parameter)

saveRDS(DT_per_parameter, paste0(code_dir,date,"/MSEA_DT_compiled_per_parameter.RDS"))



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# data table per patient --------------------------------------------------
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# patients <- as.vector(unique(DT$PatientID))
# patients <- as.vector(unique(DT$PatientID))


DT_per_patient <- data.table()
DT_per_patient[, c("Step", "Z_threshold", "Max_rxn", "PatientID", "Gene", "DBS", "Mets_in_set", "Av_rank", "Sd_rank", 
                   "Max_rank", "Min_rank", "Missed", "Av_Rev_Rel_rank", "Sd_Rev_Rel_rank",
                   "Av_Rel_rank", "Sd_Rel_rank","P.value_check", "Max_Last_position", "Min_Last_position") := 
                 DT[Include == TRUE, list(
                   unique(Gene),
                   unique(DBS),
                   unique(Mets_in_set),
                   mean(Position), # Av_rank
                   sd(Position), # Sd_rank
                   max(Position), # Max_rank
                   min(Position), # Min_rank
                   sum(P.value == 1), # Missed
                   mean(Rev.Rank.frac), # Av_Rev_Rel_rank
                   sd(Rev.Rank.frac), # Sd_Rev_Rel_rank
                   mean(Rank.frac), # Av_Rel_rank
                   sd(Rank.frac), # Sd_Rel_rank
                   length(unique(P.value)) == 1, # P.value_check
                   max(Last_position), # Max_tot_gen
                   min(Last_position) # Min_tot_gen
                 ), by = .(Step, Z_threshold, Max_rxn, PatientID)]]
# DT_per_patient[, Max_rxn:=factor(Max_rxn, levels = max_rxns)]

DT_per_patient[, Av_rank_excl_miss := ifelse(Missed > 0, NA, Av_rank)]
DT_per_patient[, Sd_rank_excl_miss := ifelse(Missed > 0, NA, Sd_rank)]
DT_per_patient[, Dataset := str_replace_all(PatientID, "P[0-9]+\\^","")]

# DT_per_patient[DT_per_patient$Missed > 0, Sd_rank := -1]
# DT_per_patient[, Colour := ifelse(DT_per_patient[,Transporter], "Red","Black")]

# DT_per_patient[, c("Sd_rank", "Sd_Rev_Rel_rank", "Sd_Rel_rank") := ifelse(Missed, -1, Sd_rank)]

# For some reason, probably something with the list() creation of the DT, the DT doesn't correctly register as one
# So this 'rectify' is necessary
DT_per_patient <- data.table(DT_per_patient)



saveRDS(DT_per_patient, paste0(code_dir,date,"/MSEA_DT_compiled_per_patient.RDS"))



