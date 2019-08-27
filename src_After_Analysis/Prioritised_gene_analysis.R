###########################################################################
# Info --------------------------------------------------------------------
###########################################################################

# This script uses the output of Summarise_MSEA_results.R to generate graphic representations of
# the Crossomics results

# R version:  3.6.0 (2019-04-26)
# platform:   x86_64-apple-darwin15.6.0 (64-bit)
# OS:         macOS Mojave 10.14.6
# 
# libraries:
# rstudioapi    0.10
# data.table    1.12.2
# ggplot2       3.1.1
# RColorBrewer  1.1-2
# scales        1.0.0


###########################################################################
# Libraries ---------------------------------------------------------------
###########################################################################

library(rstudioapi)
library(data.table)
library(RColorBrewer)
library(ggplot2)
library(scales)



###########################################################################
# Read data ---------------------------------------------------------------
###########################################################################

code_dir <- paste0(dirname(rstudioapi::getActiveDocumentContext()$path),"/../Results/")
date <- "2019-08-15"
DT <- data.table::as.data.table(readRDS(paste0(code_dir,date,"/MSEA_DT_compiled.RDS")))


# temp$size_f = factor(temp$size, levels=c('50%','100%','150%','200%'))
Z_thresholds <- levels(DT$Z_threshold)
max_rxns <- levels(DT$Max_rxn)
steps <- levels(DT$Step)
subset_patients <- TRUE

if(subset_patients){
  patients_excluded <- c("P56","P57","P58","P59","P68")
  DT <- DT[!Patient %in% patients_excluded]
  sub_name <- "sub_P"
} else {
  sub_name <- "all_P"
}


# seeds <- c(2341, 6734892, 83, 698, 991)


###########################################################################
# Create all datatables ---------------------------------------------------
###########################################################################

###########################################################################
# Correct / incorrect prioritisation datatable ----------------------------
###########################################################################
# gene correctly prioritised?
DT$Prioritised40 <- DT$Position <= 40
DT$Prioritised20 <- DT$Position <= 20
DT$Prioritised10 <- DT$Position <= 10
DT$Prioritised5 <- DT$Position <= 5
DT$Prioritised50 <- DT$Position <= 50
DT[, PatientID:=do.call(paste0,.SD), .SDcols=c("Patient","Gene")]

DT_noTrans <- DT
DT_noTrans <- DT_noTrans[Transporter==FALSE,,]

outdir_name <- paste0("/Users/mkerkho7/DIMS2_repo/Crossomics/TestResults/",Sys.Date())
if (!file.exists(outdir_name))dir.create(outdir_name, recursive = TRUE)

# for(z in c(TRUE, FALSE)){
z <- TRUE
  NoTrans <- z
  
  # tmpDT <- DT_noTrans
  # tmpDT <- DT
  if(NoTrans) tmpDT <- DT_noTrans else tmpDT <- DT
  var_name <- ifelse(NoTrans, "NoTrans", "AllGenes")
  
  DT_per_parameter <- NULL
  # for(maxrxn in max_rxns){
  #   for(threshs in Z_thresholds){
  #     thresh <- unlist(strsplit(threshs, ", "))
  #     for(step in steps){
  #       # There are 255 rows per combination (5seeds * 51 patients/disease genes)
  #       # # fraction of disease genes in the parameter combination to be in top x (prior = prioritised)
  #       # Prior.frac40 <- mean(tmpDT[Step==step & Z_threshold==threshs & Max_rxn==maxrxn, Prioritised40])
  #       # Prior.frac20 <- mean(tmpDT[Step==step & Z_threshold==threshs & Max_rxn==maxrxn, Prioritised20])
  #       # Prior.frac10 <- mean(tmpDT[Step==step & Z_threshold==threshs & Max_rxn==maxrxn, Prioritised10])
  #       # Prior.frac5 <- mean(tmpDT[Step==step & Z_threshold==threshs & Max_rxn==maxrxn, Prioritised5])
  #       # prior_av <- mean(tmpDT[Step==step & Z_threshold==threshs & Max_rxn==maxrxn, Rev_Pos_frac])
  #       # # fraction of disease genes missed (p=1)
  #       # missed_dis <- mean(tmpDT[Step==step & Z_threshold==threshs & Max_rxn==maxrxn, Position==Total_genes])
  #       # Prior.sd40 <- sd(tmpDT[Step==step & Z_threshold==threshs & Max_rxn==maxrxn, Prioritised40])
  #       # Prior.sd20 <- sd(tmpDT[Step==step & Z_threshold==threshs & Max_rxn==maxrxn, Prioritised20])
  #       # Prior.sd10 <- sd(tmpDT[Step==step & Z_threshold==threshs & Max_rxn==maxrxn, Prioritised10])
  #       # Prior.sd5 <- sd(tmpDT[Step==step & Z_threshold==threshs & Max_rxn==maxrxn, Prioritised5])
  #       # sd_av <- sd(tmpDT[Step==step & Z_threshold==threshs & Max_rxn==maxrxn, Rev_Pos_frac])
  #       # sd_missed <- sd(tmpDT[Step==step & Z_threshold==threshs & Max_rxn==maxrxn, Position==Total_genes])
  #       # 
  #       # # Out/inside top x genes:
  #       # out_top50 <- sum(!tmpDT[Step==step & Z_threshold==threshs & Max_rxn==maxrxn, Prioritised50])
  #       # in_top50 <- sum(tmpDT[Step==step & Z_threshold==threshs & Max_rxn==maxrxn, Prioritised50])
  #       # av_top50 <- mean(tmpDT[Step==step & Z_threshold==threshs & Max_rxn==maxrxn & Prioritised50, Position])
  #       # sd_top50 <- sd(tmpDT[Step==step & Z_threshold==threshs & Max_rxn==maxrxn & Prioritised50, Position])
  #       
  #       DT_var_specific <- tmpDT[Step==step & Z_threshold==threshs & Max_rxn==maxrxn]
  #       # fraction of disease genes in the parameter combination to be in top x (prior = prioritised)
  #       Prior.frac40 <- mean(DT_var_specific[, Prioritised40])
  #       Prior.frac20 <- mean(DT_var_specific[, Prioritised20])
  #       Prior.frac10 <- mean(DT_var_specific[, Prioritised10])
  #       Prior.frac5 <- mean(DT_var_specific[, Prioritised5])
  #       prior_av <- mean(DT_var_specific[, Rev_Pos_frac])
  #       # fraction of disease genes missed (p=1)
  #       missed_dis <- mean(DT_var_specific[, Position==Total_genes])
  #       Prior.sd40 <- sd(DT_var_specific[, Prioritised40])
  #       Prior.sd20 <- sd(DT_var_specific[, Prioritised20])
  #       Prior.sd10 <- sd(DT_var_specific[, Prioritised10])
  #       Prior.sd5 <- sd(DT_var_specific[, Prioritised5])
  #       sd_av <- sd(DT_var_specific[, Rev_Pos_frac])
  #       sd_missed <- sd(DT_var_specific[, Position==Total_genes])
  #       
  #       # Out/inside top x genes:
  #       out_top50 <- sum(!DT_var_specific[, Prioritised50])
  #       in_top50 <- sum(DT_var_specific[, Prioritised50])
  #       av_top50 <- mean(DT_var_specific[Prioritised50 == TRUE, Position])
  #       sd_top50 <- sd(DT_var_specific[Prioritised50 == TRUE, Position])
  #       
  #       varDT <- data.table(Z_threshold = threshs,
  #                           Step = step,
  #                           Max_rxn = maxrxn,
  #                           
  #                           Prior.frac40 = Prior.frac40,
  #                           Prior.frac20 = Prior.frac20,
  #                           Prior.frac10 = Prior.frac10,
  #                           Prior.frac5 = Prior.frac5,
  #                           Prior.sd40 = Prior.sd40,
  #                           Prior.sd20 = Prior.sd20,
  #                           Prior.sd10 = Prior.sd10,
  #                           Prior.sd5 = Prior.sd5,
  #                           
  #                           Prior_av = prior_av,
  #                           Sd_av = sd_av,
  #                           
  #                           Missed.frac = missed_dis,
  #                           Sd_missed = sd_missed,
  #                           
  #                           Out_top50 = out_top50,
  #                           In_top50 = in_top50,
  #                           Av_top50 = av_top50,
  #                           Sd_top50 = sd_top50
  #       )
  #       DT_per_parameter <- rbind(DT_per_parameter, varDT)
  #     }
  #   }
  # }
  # DT_per_parameter[, Max_rxn:=factor(Max_rxn, levels = max_rxns)]
  
#   # return 1 or 0 depending on whether the parameter combination is the best scoring or not
#   DT_per_parameter[,"best40" := ifelse(Prior.frac40 == max(Prior.frac40), 1, 0)]
#   DT_per_parameter[,"best20" := ifelse(Prior.frac20 == max(Prior.frac20), 1, 0)]
#   DT_per_parameter[,"best10" := ifelse(Prior.frac10 == max(Prior.frac10), 1, 0)]
#   DT_per_parameter[,"best5" := ifelse(Prior.frac5 == max(Prior.frac5), 1, 0)]
#   
#   # DT_per_parameter[,"signif40"] <- ifelse(DT_per_parameter[,"Prior.frac40"] == max(DT_per_parameter[,"Prior.frac40"]),1,0)
#   # DT_per_parameter[,"signif20"] <- ifelse(DT_per_parameter[,"Prior.frac20"] == max(DT_per_parameter[,"Prior.frac20"]),1,0)
#   # DT_per_parameter[,"signif10"] <- ifelse(DT_per_parameter[,"Prior.frac10"] == max(DT_per_parameter[,"Prior.frac10"]),1,0)
#   # DT_per_parameter[,"signif5"] <- ifelse(DT_per_parameter[,"Prior.frac5"] == max(DT_per_parameter[,"Prior.frac5"]),1,0)
#   DT_per_parameter[,"Frac.inTop50" := In_top50/(In_top50+Out_top50)]
#   
#   DT_tmp <- DT_per_parameter[,Av_top50, Frac.inTop50]
#   min_values <- unique(head(sort(DT_tmp[Frac.inTop50 >= 0.5, Av_top50]),5))
#   
#   DT_per_parameter[,"best_av50" := Frac.inTop50 >= 0.5 & Av_top50 %in% min_values]
#   DT_per_parameter[,"best_av50"] <- DT_per_parameter$Frac.inTop50 >= 0.5 & DT_per_parameter$Av_top50 %in% min_values
#   # DT_per_parameter[best_av50 == TRUE, best_order := order(DT_per_parameter[DT_per_parameter$best_av50 == TRUE,"Av_top50"])]
#   DT_per_parameter[best_av50 == TRUE, best_order := rank(DT_per_parameter[best_av50 == TRUE, Av_top50])]
#   
# }
  
# New version of creating the complete datatable without for loops, using the power of data.table:
  
  if(NoTrans) tmpDT <- DT_noTrans else tmpDT <- DT
  
DT_tmp1 <- data.table()
DT_tmp2 <- data.table()
DT_tmp1[, c("Step", "Z_threshold", "Max_rxn", "Prior.frac40","Prior.frac20","Prior.frac10","Prior.frac5","Prior.pos.frac.av","Missed","Missed.frac",
            "Max_Tot.Genes","Min_Tot.Genes","Prior.sd40", "Prior.sd20", "Prior.sd10", "Prior.sd5","Prior.pos.frac.sd","Missed.sd", "Out_top50", "In_top50") := 
          tmpDT[, list(
            mean(Prioritised40), # Prior.frac40
            mean(Prioritised20), # Prior.frac20
            mean(Prioritised10), # Prior.frac10
            mean(Prioritised5), # Prior.frac5
            mean(Rev_Pos_frac), # Prior.pos.frac.av
            sum(P.value==1)/5, # Missed; if a gene is missed in 1 seed, it is missed in all (5) --> count once per parametercombination
            mean(P.value==1), # Missed.frac
            max(Total_genes), # Max_Tot.Genes
            min(Total_genes), # Min_Tot.Genes
            sd(Prioritised40), # Prior.sd40
            sd(Prioritised20), # Prior.sd20
            sd(Prioritised10), # Prior.sd10
            sd(Prioritised5), # Prior.sd5
            sd(Rev_Pos_frac), # Prior.pos.frac.sd
            sd(Position==Total_genes), # Missed.sd
            sum(!Prioritised50), # Out_top50
            sum(Prioritised50) # In_top50
            ), by = .(Step, Z_threshold, Max_rxn)]]
DT_tmp2[, c("Step", "Z_threshold", "Max_rxn","Av_top50","Sd_top50") := tmpDT[Prioritised50 == TRUE, list(mean(Position), sd(Position)), by = .(Step, Z_threshold, Max_rxn)]]
DT_per_parameter <- merge(DT_tmp1, DT_tmp2, by = c("Step","Z_threshold", "Max_rxn"))
rm(DT_tmp1,DT_tmp2)

DT_per_parameter[, Max_rxn:=factor(Max_rxn, levels = max_rxns)]

# return 1 or 0 depending on whether the parameter combination is the best scoring or not
DT_per_parameter[,c("best40","best20","best10","best5") := list(
  ifelse(Prior.frac40 == max(Prior.frac40), 1, 0),
  ifelse(Prior.frac20 == max(Prior.frac20), 1, 0),
  ifelse(Prior.frac10 == max(Prior.frac10), 1, 0),
  ifelse(Prior.frac5 == max(Prior.frac5), 1, 0))]
DT_per_parameter[,"Frac.inTop50" := In_top50/(In_top50+Out_top50)]

DT_tmp <- DT_per_parameter[,Av_top50, Frac.inTop50]
min_values <- unique(head(sort(DT_tmp[Frac.inTop50 >= 0.5, Av_top50]),5))

DT_per_parameter[,"best_av50" := Frac.inTop50 >= 0.5 & Av_top50 %in% min_values]
DT_per_parameter[,"best_av50"] <- DT_per_parameter$Frac.inTop50 >= 0.5 & DT_per_parameter$Av_top50 %in% min_values
# DT_per_parameter[best_av50 == TRUE, best_order := order(DT_per_parameter[DT_per_parameter$best_av50 == TRUE,"Av_top50"])]
DT_per_parameter[best_av50 == TRUE, best_order := rank(DT_per_parameter[best_av50 == TRUE, Av_top50])]
  
# For some reason, probably something with the list() creation of the DT, the DT doesn't correctly register as one
# So this 'rectify' is necessary
DT_per_parameter <- data.table(DT_per_parameter)




# patients <- as.vector(unique(DT$PatientID))
patients <- as.vector(unique(tmpDT$PatientID))


DT_per_patient <- data.table()
DT_per_patient[, c("Step", "Z_threshold", "Max_rxn", "PatientID", "Av_rank", "Sd_rank", "Max_rank", "Min_rank","Missed", "Out_top50", "Av_Rel_rank", "Sd_Rel_rank",
            "Transporter", "P.value_check", "Max_tot_gen", "Min_tot_gen") := 
          tmpDT[, list(
            mean(Position), # Av_rank
            sd(Position), # Sd_rank
            max(Position), # Max_rank
            min(Position), # Min_rank
            sum(P.value == 1), # Missed
            sum(Position > 50), # Out_top50
            mean(Rev_Pos_frac), # Av_Rel_rank
            sd(Rev_Pos_frac), # Sd_Rel_rank
            unique(Transporter), # Transporter
            length(unique(P.value)) == 1, # P.value_check
            max(Total_genes), # Max_tot_gen
            min(Total_genes) # Min_tot_gen
          ), by = .(Step, Z_threshold, Max_rxn, PatientID)]]
DT_per_patient[, Max_rxn:=factor(Max_rxn, levels = max_rxns)]
DT_per_patient[DT_per_patient$Missed == 5, Sd_rank := -1]
DT_per_patient[, Colour := ifelse(DT_per_patient[,Transporter], "Red","Black")]

# For some reason, probably something with the list() creation of the DT, the DT doesn't correctly register as one
# So this 'rectify' is necessary
DT_per_patient <- data.table(DT_per_patient)



###########################################################################
# Make plot of effectiveness of variables ---------------------------------
###########################################################################
# 
# 
# 
# # Violin plots of relative rank -------------------------------------------
# 
# 
# # p <- ggplot(DT, aes(x=Z_threshold, y=Pos_frac, color=as.factor(Max_rxn))) +
# p <- ggplot(DT, aes(x=Z_threshold, y=Pos_frac, color=Protein_function)) +
#   geom_boxplot() + 
#   geom_point(position = "jitter") +
#   # geom_point(position = "jitter", aes(shape=as.factor(Step))) +
#   scale_x_discrete(name ="Z-value threshold")
# # p + facet_grid(. ~ Step)
# p <- ggplot(DT, aes(x = factor(Step), y = Rev_Pos_frac, fill = Transporter)) +
#   geom_boxplot(width=0.4, fatten = 3, outlier.size = 0.1) +
#   geom_violin(fill = "white", alpha = 0.1) +
#   geom_dotplot(binaxis = "y", stackdir = "centerwhole", binwidth = 0.015) +
#   scale_x_discrete(name ="Z-value threshold") +
#   theme_light()
# # p <- ggplot(DT, aes(x = factor(Step), y = Position, fill = Transporter)) +
# #   geom_boxplot(width=0.1) +
# #   geom_violin(fill = "white", alpha = 0.1) +
# #   geom_dotplot(binaxis = "y", stackdir = "centerwhole", binwidth = 1.5) +
# #   scale_x_discrete(name ="Z-value threshold")
#   
#   
# # png(file = "/Users/mkerkho7/DIMS2_repo/Crossomics/TestResults/AllPars.png", width = 1200, height = 800, res=1200)
# p + facet_grid(Max_rxn ~ Z_threshold, labeller = label_both)
# try(dev.off(), silent = TRUE)
# ggsave("/Users/mkerkho7/DIMS2_repo/Crossomics/TestResults/AllPars.png", 
#        width = 300, height = 200, dpi=600, units = "mm")
# 
# 
# 
# 
# # Subset of rows ----------------------------------------------------------
# 
# # DT_small <- DT[Z_threshold == "-1, 1.5" & (Step == 3) & (Max_rxn == 10 | Max_rxn == 12)]
# # No transporter disease genes
# DT_small <- DT[Transporter != TRUE]
# p <- ggplot(DT_small, aes(x = factor(Step), y = Rev_Pos_frac, fill = Gene)) +
#   stat_summary(fun.data="mean_sdl",  fun.args = list(mult=1), 
#                geom="pointrange", color = "red") +
#   geom_boxplot(fill = "red", alpha = 0.1, width=0.2, fatten = 5) +
#   geom_violin(fill = "white", alpha = 0.1) +
#   geom_dotplot(binaxis = "y", stackdir = "centerwhole", binwidth = 0.008) +
#   scale_x_discrete(name ="Z-value threshold") +
#   theme_linedraw()
# png(file = "/Users/mkerkho7/DIMS2_repo/Crossomics/TestResults/AllParsNoTransporters.png", width = 2000, height = 2000)
# p + facet_grid(Max_rxn ~ Z_threshold)
# try(dev.off(), silent = TRUE)
# 
# png(file = "/Users/mkerkho7/DIMS2_repo/Crossomics/TestResults/ByProtFunc.png", width = 800, height = 400)
# ggplot(DT, aes(x = factor(Step), y = Rev_Pos_frac, fill = Protein_function)) +
#   # geom_dotplot(binaxis = "y", stackdir = "center", binwidth = 0.008) +
#   geom_boxplot() +
#   theme_dark()
# try(dev.off(), silent = TRUE)
# 
# DT_small$Z_threshold_order <- factor(DT_small$Z_threshold, levels = Z_thresholds)
# 
# q <- ggplot(DT_small, aes(x = factor(Step), y = Pos_frac)) +
#   geom_dotplot(binaxis = "y", stackdir = "center", binwidth = 0.005) +
#   scale_x_discrete(name ="Z-value threshold") +
#   geom_violin(alpha = 0.1)
# 
# 
# # library(ggrepel)
# 
# bw <- 0.005
# 
# built <- ggplot_build(q)
# point.pos <- built$data[[1]]
# 
# # Order rows 
# idx <- order(DT_small[,Pos_frac])
# DT_small <- DT_small[idx]
# 
# # Get the dimensions of the target device in pixels
# size <- dev.size(units = 'px')
# # Get the range of x and y domain values
# extent <- with(built$layout$panel_params[[1]], abs(c(diff(x.range), diff(y.range))))
# 
# DT_small$ytext <- point.pos$y
# DT_small$xtext <- point.pos$x
# 
# try(dev.off(), silent = TRUE)
# png(file = "/Users/mkerkho7/DIMS2_repo/Crossomics/TestResults/testRplotsmallstep3.png", width = 500, height = 1200)
# q <- ggplot(DT_small, aes(x = factor(Step), y = Pos_frac)) +
#   geom_dotplot(binaxis = "y", binwidth = bw, stackdir = 'center') +
#   geom_text_repel(
#     aes(xtext, ytext, label = paste(DT_small$Patient,DT_small$Gene)),
#     box.padding = unit(2 * size[1] * bw / extent[1], 'points'),
#     color = 'red'
#   ) +
#   scale_x_discrete(name ="Step") +
#   geom_violin(alpha = 0.1)
# q + facet_grid(Max_rxn ~ .)
# try(dev.off(), silent = TRUE)
# 


  ###########################################################################
  # Combination plot of correctly prioritised and missed genes --------------
  ###########################################################################
  
  # Colour scheme to use in comparison plots
  my_palette = rev(brewer.pal(5, "Greens"))[c(2:5)]
  
  
  # png(file = "/Users/mkerkho7/DIMS2_repo/Crossomics/TestResults/20190815PriorMissed_noTrans.png", width = 1600, height = 1400)
  
p <- ggplot(DT_per_parameter, aes(x = Step)) +
  geom_line(aes(y = Prior.frac40, colour = my_palette[4], group = 1)) +
  
  # Add shaded area for what is expected if the results were random:
  geom_ribbon(aes(x= as.numeric(Step), ymin=1*40/Max_Tot.Genes, ymax=1*40/Min_Tot.Genes), colour = adjustcolor(my_palette[4],alpha.f=0.5), fill = my_palette[4], alpha="0.5") +
  geom_ribbon(aes(x= as.numeric(Step), ymin=1*20/Max_Tot.Genes, ymax=1*20/Min_Tot.Genes), colour = adjustcolor(my_palette[3],alpha.f=0.5), fill = my_palette[3], alpha="0.5") +
  geom_ribbon(aes(x= as.numeric(Step), ymin=1*10/Max_Tot.Genes, ymax=1*10/Min_Tot.Genes), colour = adjustcolor(my_palette[2],alpha.f=0.5), fill = my_palette[2], alpha="0.5") +
  geom_ribbon(aes(x= as.numeric(Step), ymin=1*5/Max_Tot.Genes, ymax=1*5/Min_Tot.Genes), colour = adjustcolor(my_palette[1],alpha.f=0.6), fill = my_palette[1], alpha="0.6") +
  
  geom_line(aes(y = Prior.frac40, colour = my_palette[4], group = 1)) +
  geom_point(aes(y = Prior.frac40, colour = my_palette[4]), size=0.5) +
  geom_line(aes(y = Prior.frac20, colour = my_palette[3], group = 1)) +
  geom_point(aes(y = Prior.frac20, colour = my_palette[3]), size=0.5) +
  geom_line(aes(y = Prior.frac10, colour = my_palette[2], group = 1)) +
  geom_point(aes(y = Prior.frac10, colour = my_palette[2]), size=0.5) +
  geom_line(aes(y = Prior.frac5, colour = my_palette[1], group = 1)) +
  geom_point(aes(y = Prior.frac5, colour = my_palette[1]), size=0.5) +
  geom_line(aes(y = Missed.frac, colour = "Red", group = 1)) +
  geom_point(aes(y = Missed.frac, colour = "Red"), size=0.5) 
  # geom_point(aes(x= Step, y = 100*best_order, shape = "Best"), size=8) + scale_shape_manual(values = 42)

  
  # theme_dark() + 
  # ylab("Disease genes / Total dis. genes") +
  # xlab("Max distance to primary reaction") +
  # scale_color_manual(name = "Prioritised Genes", 
  #                    labels = c("In Top 5", "In Top 10", "In Top 20", "In Top 40", "Missed"),
  #                    values=c(my_palette[1], my_palette[2], my_palette[3], my_palette[4],'RED',"black"),
  #                    guide = guide_legend(override.aes = list(
  #                      shape = c(rep(NA, 5), "*")
  #                    )))
# scale_colour_manual(values = c("purple", "green", "blue", "yellow", "magenta","orange", "cyan", "red", "black"),
#                     guide = guide_legend(override.aes = list(
#                       linetype = c(rep("blank", 7), "solid", "dashed"),
#                       shape = c(rep(16, 7), NA, NA))))
p <- p + facet_grid(Max_rxn ~ Z_threshold, labeller = label_both) +
  # geom_point(aes(x= Step, y = 100*best_order, shape = "Best"), size=8) + scale_shape_manual(values = 42)
  geom_point(data = DT_per_parameter[DT_per_parameter$best40 ==1, ],aes(x=Step, y=Prior.frac40, shape = "Best"), size=8) + 
  scale_shape_manual(name = "",
                     labels = "Best param.\nperformance",
                     values = 42) +
  # geom_point(data = DT_per_parameter[DT_per_parameter$best40 ==1, ],aes(x=Step, y=Prior.frac40),shape = "*", size=8, show.legend = FALSE, colour = "black", fill = "black") +
  geom_point(data = DT_per_parameter[DT_per_parameter$best20 ==1, ],aes(x=Step, y=Prior.frac20),shape = "*", size=8, show.legend = FALSE, colour = "black", fill = "black") +
  geom_point(data = DT_per_parameter[DT_per_parameter$best10 ==1, ],aes(x=Step, y=Prior.frac10),shape = "*", size=8, show.legend = FALSE, colour = "black", fill = "black") +
  geom_point(data = DT_per_parameter[DT_per_parameter$best5 ==1, ],aes(x=Step, y=Prior.frac5),shape = "*", size=8, show.legend = FALSE, colour = "black", fill = "black")
p <- p + theme_dark() + 
  ylab("Disease genes / Total dis. genes") +
  xlab("Max distance to primary reaction") +
  scale_color_manual(name = "Prioritised Genes",
                     labels = c("In Top 5", "In Top 10", "In Top 20", "In Top 40", "Missed"),
                     values=c(my_palette[1], my_palette[2], my_palette[3], my_palette[4],'RED'))
try(dev.off(), silent = TRUE)
# ggsave("/Users/mkerkho7/DIMS2_repo/Crossomics/TestResults/20190823PriorMissed_NoTrans.png", plot = p,
# ggsave("/Users/mkerkho7/DIMS2_repo/Crossomics/TestResults/20190823PriorMissed_AllGenes.png", plot = p,
ggsave(paste0(outdir_name,"/PriorMissed_",var_name,"_",sub_name,".png"), plot = p,
       width = 300, height = 200, dpi=600, units = "mm")



  ###########################################################################
  # Plot of average of top 50 genes -----------------------------------------
  ###########################################################################
  # DT_per_parameter[,"signif40"] <- ifelse(DT_per_parameter[,"Prior.frac40"] == max(DT_per_parameter[,"Prior.frac40"]),1,0)
  
  my_sig_palette = rev(brewer.pal(6, "RdYlGn"))[2:6]
  
  p <- ggplot(DT_per_parameter, aes(label=Av_top50)) +
    geom_line(aes(x = Step, y = Av_top50, colour = "Average", group = 1)) +
    geom_point(aes(x = Step, y = Av_top50, colour = "Average"), size=0.5) +
    # geom_errorbar(aes(x = Step, ymax = Av_top50 + Sd_top50, ymin = Av_top50 - Sd_top50), position = "dodge") +
    # geom_line(aes(x = Step, y = Out_top50/4, colour = "Red", group = 1)) +
    # geom_point(aes(x = Step, y = Out_top50/4, colour = "Red"), size=0.5) +
    # geom_line(aes(x = Step, y = In_top50/4, colour = my_palette[1], group = 1)) +
    # geom_point(aes(x = Step, y = In_top50/4, colour = my_palette[1]), size=0.5) +
    geom_line(aes(x = Step, y = In_top50/(In_top50+Out_top50)*40, colour = "Frac.InTop", group = 1)) +
    geom_point(aes(x = Step, y = In_top50/(In_top50+Out_top50)*40, colour = "Frac.InTop"), size=0.5) +
    scale_y_continuous(limits = c(0,40), sec.axis = sec_axis(~./40, name = "Frac. genes in top 50")) +
    # geom_line(aes(x = Step, y = In_top50, colour = my_palette[3], group = 1)) +
    # geom_point(aes(x = Step, y = In_top50, colour = my_palette[3]), size=0.5) +
    theme_dark() +
    scale_color_manual(name = "Legend",
                       labels = c("Av. rank of top 50", "Frac. in top 50"),
                       values=c('Black', 'RED'))
  p <- p + facet_grid(Max_rxn ~ Z_threshold, labeller = label_both) +
    labs(y = "Average disease gene rank")
    # geom_point(data = DT_per_parameter[best_av50==TRUE, ],aes(x=Step, y=Av_top50),shape = "*", size=8, show.legend = FALSE, colour = "black")+
  for(i in c(1:5)){
    p <- p + geom_point(data = DT_per_parameter[as.vector(DT_per_parameter[,"best_order"]==i),],aes(x=Step, y=Av_top50),shape = "*", size=8, show.legend = FALSE, colour = my_sig_palette[i])
    p <- p + geom_text(data = DT_per_parameter[as.vector(DT_per_parameter[,"best_order"]==i),],
                       aes(x=Step, y=Av_top50, label = signif(Av_top50, digits = 5), group = best_order), 
                       size=3, 
                       # show.legend = FALSE, 
                       colour = my_sig_palette[i],
                       position = position_dodge(width = 2),
                       vjust = -0.5)
  }
   
  p <- p + geom_hline(yintercept=20, linetype="dashed", color = "darksalmon")
  
  try(dev.off(), silent = TRUE)
  # ggsave("/Users/mkerkho7/DIMS2_repo/Crossomics/TestResults/20190823Top50_NoTrans.png", plot = p,
  # ggsave("/Users/mkerkho7/DIMS2_repo/Crossomics/TestResults/20190823Top50_AllGenes.png", plot = p,
  #        width = 300, height = 200, dpi=600, units = "mm")
  ggsave(paste0(outdir_name,"/Top50_",var_name,"_",sub_name,".png"), plot = p,
         width = 300, height = 200, dpi=600, units = "mm")
  
  
  
  
  ###########################################################################
  # Standard deviations of disease gene rank per patient --------------------
  ###########################################################################
  
  # DT[, PatientID:=do.call(paste0,.SD), .SDcols=c("Patient","Gene")]
  
  # if(NoTrans) tmpDT <- DT_noTrans else tmpDT <- DT
  # 
  # patients <- as.vector(unique(DT$PatientID))
  # patients <- as.vector(unique(tmpDT$PatientID))
  # 
  # 
  # DT_per_patient <- NULL
  # for(maxrxn in max_rxns){
  #   for(threshs in Z_thresholds){
  #     thresh <- unlist(strsplit(threshs, ", "))
  #     for(step in steps){
  #       for(patient in patients){
  #         varDT <- NULL
  #         DT_var_specific <- tmpDT[Step==step & Z_threshold==threshs & Max_rxn==maxrxn & PatientID == patient]
  #         # There should be 9180 rows (6 steps * 5 thresholds * 6 maxrxns * 51 patients/disease genes)
  #         # av_rank <- mean(tmpDT[Step==step & Z_threshold==threshs & Max_rxn==maxrxn & PatientID == patient, Position])
  #         # min_rank <- min(tmpDT[Step==step & Z_threshold==threshs & Max_rxn==maxrxn & PatientID == patient, Position])
  #         # max_rank <- max(tmpDT[Step==step & Z_threshold==threshs & Max_rxn==maxrxn & PatientID == patient, Position])
  #         # missed <- sum(tmpDT[Step==step & Z_threshold==threshs & Max_rxn==maxrxn & PatientID == patient, P.value] == 1)
  #         # out_top50 <- sum(tmpDT[Step==step & Z_threshold==threshs & Max_rxn==maxrxn & PatientID == patient, Position] >= 50)
  #         # sd_rank <- sd(tmpDT[Step==step & Z_threshold==threshs & Max_rxn==maxrxn & PatientID == patient, Position])
  #         # av_rel_rank <- mean(tmpDT[Step==step & Z_threshold==threshs & Max_rxn==maxrxn & PatientID == patient, Rev_Pos_frac])
  #         # sd_rel_rank <- sd(tmpDT[Step==step & Z_threshold==threshs & Max_rxn==maxrxn & PatientID == patient, Rev_Pos_frac])
  #         # transporter <- unique(tmpDT[Step==step & Z_threshold==threshs & Max_rxn==maxrxn & PatientID == patient, Transporter])
  #         # P.value_check <- length(unique(tmpDT[Step==step & Z_threshold==threshs & Max_rxn==maxrxn & PatientID == patient, P.value]))==1
  # 
  #         av_rank <- mean(DT_var_specific[, Position])
  #         min_rank <- min(DT_var_specific[, Position])
  #         max_rank <- max(DT_var_specific[, Position])
  #         missed <- sum(DT_var_specific[, P.value] == 1)
  #         out_top50 <- sum(DT_var_specific[, Position] >= 50)
  #         sd_rank <- sd(DT_var_specific[, Position])
  #         av_rel_rank <- mean(DT_var_specific[, Rev_Pos_frac])
  #         sd_rel_rank <- sd(DT_var_specific[, Rev_Pos_frac])
  #         transporter <- unique(DT_var_specific[, Transporter])
  #         P.value_check <- length(unique(DT_var_specific[, P.value]))==1
  # 
  #         max_total_genes <- max(DT_var_specific[, Total_genes])
  #         min_total_genes <- min(DT_var_specific[, Total_genes])
  # 
  #         varDT <- data.table(Av_rank = av_rank,
  #                             Min_rank = min_rank,
  #                             Max_rank = max_rank,
  #                             Missed = missed,
  #                             Out_top50 = out_top50,
  # 
  #                             Max_tot_gen = max_total_genes,
  #                             Min_tot_gen = min_total_genes,
  # 
  #                             Sd_rank = sd_rank,
  #                             Av_Rel_rank = av_rel_rank,
  #                             Sd_Rel_rank = sd_rel_rank,
  # 
  #                             PatientID = patient,
  #                             Z_threshold = threshs,
  #                             Step = step,
  #                             Max_rxn = maxrxn,
  #                             Transporter = transporter,
  #                             P.value_check = P.value_check)
  # 
  #         DT_per_patient <- rbind(DT_per_patient, varDT)
  # 
  #       }
  #     }
  #   }
  # }
  # DT_per_patient[, Max_rxn:=factor(Max_rxn, levels = max_rxns)]
  # DT_per_patient[DT_per_patient$Missed == 5, Sd_rank := -1]

  # DT_per_patient[, Colour := ifelse(DT_per_patient[,Transporter], "Red","Black")]

  # p <- ggplot(data = DT_per_patient, aes(x = Step, y = Sd_rank, fill = PatientID)) +
  #   geom_bar(position = "dodge", stat="identity") +
  #   facet_grid(Max_rxn ~ Z_threshold, labeller = label_both)
  # ggsave(paste0(outdir_name,"/Sd_rank_PerPatient_",var_name,".png"), plot = p,
  #        width = 300, height = 200, dpi=600, units = "mm")
  # 
  # p <- ggplot(data = DT_per_patient, aes(x = Step, y = Sd_Rel_rank, fill = PatientID)) +
  #   geom_bar(position = "dodge", stat="identity") +
  #   facet_grid(Max_rxn ~ Z_threshold, labeller = label_both)
  # ggsave(paste0(outdir_name,"/Sd_Rel_rank_PerPatient_",var_name,".png"), plot = p,
  #        width = 300, height = 200, dpi=600, units = "mm")
  
  
  # DT_per_parameter <- DT_per_patient[,sum(Missed > 0), by = c("Step","Max_rxn", "Z_threshold")]
  # names(DT_per_parameter)[4] <- "Tot.Missed"
  
  p <- ggplot(data = DT_per_patient, aes(x = Step, y = Sd_rank)) +
    geom_bar(position = "dodge", stat="identity", aes(fill = PatientID)) +
    geom_line(data = DT_per_parameter, aes(x = Step, y = Missed/4, group = 1, colour = "Genes Missed")) +
    scale_y_continuous(sec.axis = sec_axis(~.*4, name = "Tot.Genes Missed")) +
    facet_grid(Max_rxn ~ Z_threshold, labeller = label_both)
  ggsave(paste0(outdir_name,"/Sd_rank_PerPatient_",var_name,"_",sub_name,".png"), plot = p,
         width = 300, height = 200, dpi=600, units = "mm")
  
  p <- ggplot(data = DT_per_patient, aes(x = Step, y = Sd_Rel_rank)) +
    geom_bar(position = "dodge", stat="identity", aes(fill = PatientID)) +
    geom_line(data = DT_per_parameter, aes(x = Step, y = Missed/800, group = 1, colour = "Genes Missed")) +
    scale_y_continuous(sec.axis = sec_axis(~.*800, name = "Tot.Genes Missed")) +
    facet_grid(Max_rxn ~ Z_threshold, labeller = label_both) 
  ggsave(paste0(outdir_name,"/Sd_Rel_rank_PerPatient_",var_name,"_",sub_name,".png"), plot = p,
         width = 300, height = 200, dpi=600, units = "mm")
  
}

# p <- ggplot(data = DT_per_patient, aes(x = Step, y = Av_Rel_rank)) +
#   geom_bar(position = "dodge", stat="identity", aes(fill = PatientID)) +
#   # theme(axis.text.x = element_text(angle = 90))
#   geom_line(data = DT_per_patient, aes(x = Step, y = Av_Rel_rank, group = 1, colour = "Genes Missed")) +
#   # geom_line(data = DT_per_parameter, aes(x = Step, y = Tot.Missed/4, group = 1, colour = "Genes Missed")) +
#   # scale_y_continuous(sec.axis = sec_axis(~.*4, name = "Tot.Genes Missed")) +
#   facet_grid(Max_rxn ~ Z_threshold, labeller = label_both)
# ggsave(paste0(outdir_name,"/Av_Rel_rank_PerPatient_",var_name,"_",sub_name,".png"), plot = p,
#        width = 300, height = 200, dpi=600, units = "mm")

