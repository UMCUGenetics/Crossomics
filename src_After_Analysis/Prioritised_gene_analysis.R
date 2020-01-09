#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Info --------------------------------------------------------------------
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

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


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Libraries ---------------------------------------------------------------
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

library(rstudioapi)
library(data.table)
library(RColorBrewer) # Colour schemes
library(ggplot2)
library(scales)
# library(ggforce) # zoom in specific parts of plot
library(grid) # manual grob adjustment
library(gtable) # manual grob adjustment
library(gridExtra)
library(stringr)

library(svglite)




#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Adjustable settings -----------------------------------------------------
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

##### Create general plot aesthetics --------------------------------------
Thresh_labs <- c("Min.-1; Max.1.5","Min.-1.5; Max.2","Min.-3; Max.3","Min.-5; Max.5")
names(Thresh_labs) <- c("-1, 1.5", "-1.5, 2", "-3, 3", "-5, 5")
Rxn_labs <- c("<= 8","<= 10","<= 12","<= 15","<= 17","<= 19")
names(Rxn_labs) <- c(8, 10, 12, 15, 17, 19)

# Colour scheme 
my_greens <- rev(brewer.pal(5, "Greens"))[c(2:5)]
my_blues <- rev(brewer.pal(5, "Blues"))[c(1:5)]
my_greens <- my_blues
my_reds <- rev(brewer.pal(5, "Reds"))[c(2:5)]
my_sig_palette <- rev(brewer.pal(6, "RdYlGn"))[2:6]
my_sig_palette <- colorRamps::green2red(5)

# Image resolution
high_res <- 600
low_res <- 100
resolution <- low_res
digit_significance <- 3

##### Other ---------------------------------------------------------------
# Exclude patients / diseases
subset_patients <- FALSE
# trainingset <- "combined" # possible: TRUE, FALSE, NULL (for all data), "combined" (for all, but separated)

# Date of data
date <- "2019-12-10"



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Pre-processing ----------------------------------------------------------
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

##### Read data -----------------------------------------------------------
code_dir <- paste0(dirname(rstudioapi::getActiveDocumentContext()$path),"/../Results/")

DT <- data.table::as.data.table(readRDS(paste0(code_dir,date,"/MSEA_DT_compiled_all.RDS")))
# DT_per_parameter_val <- data.table::as.data.table(readRDS(paste0(code_dir,date,"/MSEA_DT_validation_per_parameter.RDS")))
# DT_per_patient_val <- data.table::as.data.table(readRDS(paste0(code_dir,date,"/MSEA_DT_validation_per_patient.RDS")))
# DT_per_parameter_tra <- data.table::as.data.table(readRDS(paste0(code_dir,date,"/MSEA_DT_training_per_parameter.RDS")))
# DT_per_patient_tra <- data.table::as.data.table(readRDS(paste0(code_dir,date,"/MSEA_DT_training_per_patient.RDS")))

DT_per_parameter <- readRDS(paste0(code_dir,date,"/MSEA_DT_per_parameter.RDS"))
DT_per_parameter_tra_val_sum <- readRDS(paste0(code_dir,date,"/MSEA_DT_per_parameter_tra_val_summary.RDS"))

# DT_validation_per_parameter <- readRDS(paste0(code_dir,date,"/MSEA_DT_by_validation_per_parameter.RDS"))

# test <- readRDS(paste0(code_dir,date,"/Metabolite_set_sizes.RDS"))
# DT_met_sets <- readRDS(paste0(code_dir,date,"/Metabolite_set_sizes.RDS"))

# DT_per_parameter_tra[, rank_10 := frank(-Prior.frac10)]

# #####
# # See how many disease genes have a metabolite set
# genes_w_mets <- list.files(paste0(code_dir, "../Data/2019-08-12/maxrxn19/mss_5_HMDBtranslated"))
# genes_w_mets <- str_remove(genes_w_mets, ".RDS")
# 
# sum(!unique(DT[Include == TRUE, Gene]) %in% genes_w_mets)

# if(trainingset == "combined"){
#   train_val <- "trainingset"
#   DT1 <- data.table::as.data.table(readRDS(paste0(code_dir,date,"/MSEA_DT_compiled_",train_val,".RDS")))
#   DT1[,Set := "training"]
#   train_val <- "validationset"
#   DT2 <- data.table::as.data.table(readRDS(paste0(code_dir,date,"/MSEA_DT_compiled_",train_val,".RDS")))
#   DT2[,Set := "validation"]
#   DT <- rbind(DT1, DT2)
#   rm(DT1, DT2)
# } else {
#   if(is.null(trainingset)) {
#     train_val <- "all"
#     DT <- data.table::as.data.table(readRDS(paste0(code_dir,date,"/MSEA_DT_compiled_",train_val,".RDS")))
#     DT[,Set := "all"]
#   } else {
#     train_val <- ifelse(trainingset, "trainingset", "validationset")
#     DT <- data.table::as.data.table(readRDS(paste0(code_dir,date,"/MSEA_DT_compiled_",train_val,".RDS")))
#     ifelse(trainingset, DT[,Set := "training"], DT[,Set := "validation"])
#   }
# }

# DT <- data.table::as.data.table(readRDS(paste0(code_dir,date,"/MSEA_DT_compiled_",train_val,".RDS")))

##### Determine parameters ------------------------------------------------
# temp$size_f = factor(temp$size, levels=c('50%','100%','150%','200%'))
# Z_thresholds <- levels(DT$Z_threshold)
# max_rxns <- levels(DT$Max_rxn)
# steps <- levels(DT$Step)
# seeds <- unique(DT$Seed)

##### Exclude predetermined patients and determine output names -----------
if(subset_patients){
  patients_excluded <- c("P56","P57","P58","P59","P68")
  DT <- DT[!Patient %in% patients_excluded]
  sub_name <- "Sub_P"
} else {
  sub_name <- "All_P"
}

outdir_name <- paste0(dirname(rstudioapi::getActiveDocumentContext()$path),"/../Plots/",date)
if (!file.exists(outdir_name))dir.create(outdir_name, recursive = TRUE)




#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Functions ---------------------------------------------------------------
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

##### Manual standard grob addition ---------------------------------------
pretty_plot <- function(p, theme = "light", secondary_y_axis = TRUE){
  z <- ggplotGrob(p)
  right_lab_loc <- 5 + 2*length(Thresh_labs) - secondary_y_axis
  right_lab_bottom <- 6 + 2*length(Rxn_labs)
  top_lab_width <- 5 + 2*length(Thresh_labs) - 2
  
  # label right strip
  z <- gtable_add_cols(z, unit(1, 'cm'), right_lab_loc)
  z <- gtable_add_grob(z, 
                       list(rectGrob(gp = gpar(col = NA, fill = ifelse(theme=="dark", "lightgray", gray(0.5)))),
                            textGrob("Extension stringency", rot = -90, gp = gpar(col = ifelse(theme=="dark", "black", "white")))),
                       # 8, 15, 18, 15, name = paste(runif(2)))
                       8, right_lab_loc+1, right_lab_bottom, right_lab_loc+1, name = paste(runif(2)))
  
  # label top strip
  z <- gtable_add_rows(z, unit(1, 'cm'), 6)
  z <- gtable_add_grob(z, 
                       list(rectGrob(gp = gpar(col = NA, fill = ifelse(theme=="dark", "lightgray", gray(0.5)))),
                            textGrob("Biochemical stringency", gp = gpar(col = ifelse(theme=="dark", "black", "white")))),
                       # 4, 5, 4, 13, name = paste(runif(2)))
                       7, 5, 7, top_lab_width, name = paste(runif(2)))
  
  # margins
  # z <- gtable_add_cols(z, unit(1/8, "line"), 7)
  # z <- gtable_add_rows(z, unit(1/8, "line"), 3)
  
  return(z)
}




#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Simple plots, meant to use as explaning for the bigger, later plots -----
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Empty facet plot
p <- ggplot(DT_per_parameter, aes(x = Step)) +
  facet_grid(Max_rxn ~ Z_threshold, labeller = labeller(Max_rxn = Rxn_labs, Z_threshold = Thresh_labs)) +
  theme_dark() + 
  ylab("Disease genes / Total dis. genes") +
  xlab("Max distance to primary reaction") +
  ggtitle("Correct disease gene prioritisation") +
  scale_color_manual(name = "Prioritised Genes",
                     labels = c("In Top 15", "Missed"),
                     values=c(my_greens[1],'RED'))
z <- pretty_plot(p, theme = "dark")
ggsave(paste0(outdir_name,"/Empty_Facet_",sub_name,".png"), plot = z,
       width = 250, height = 200, dpi=resolution, units = "mm")

# Simple plot, just top 10, all parameter combinations
p <- ggplot(DT_per_parameter, aes(x = Step)) +
  geom_line(aes(y = Prior.frac10, colour = "line10", group = 1)) +
  geom_point(aes(y = Prior.frac10, colour = "line10"), size=0.5) +
  geom_line(aes(y = Missed.frac, colour = "Missed", group = 1)) +
  geom_point(aes(y = Missed.frac, colour = "Missed"), size=0.5) +
  facet_grid(Max_rxn ~ Z_threshold, labeller = labeller(Max_rxn = Rxn_labs, Z_threshold = Thresh_labs)) +
  theme_dark() + 
  ylab("Disease genes / Total dis. genes") +
  xlab("Max distance to primary reaction") +
  ggtitle("Correct disease gene prioritisation") +
  scale_color_manual(name = "Prioritised Genes",
                     labels = c("In Top 10", "Missed"),
                     values=c(my_greens[4],'RED'))
z <- pretty_plot(p, theme = "dark")
ggsave(paste0(outdir_name,"/",train_val,"_Ranks_And_Missed_Top15_",sub_name,".png"), plot = z,
       width = 300, height = 200, dpi=resolution, units = "mm")

# Single panel from big facet plot (top 2, 5, 10 and 15)
p <- ggplot(DT_per_parameter[Z_threshold == "-3, 3" & Max_rxn == 15], aes(x = Step)) +
  geom_line(aes(y = Prior.frac15, colour = "line15", group = 1)) +
  
  # geom_ribbon(aes(x= as.numeric(Step), ymin=1*15/Max_Tot.Genes, ymax=1*15/Min_Tot.Genes, fill = "band15")) +
  # geom_ribbon(aes(x= as.numeric(Step), ymin=1*10/Max_Tot.Genes, ymax=1*10/Min_Tot.Genes, fill = "band10")) +
  # geom_ribbon(aes(x= as.numeric(Step), ymin=1*5/Max_Tot.Genes, ymax=1*5/Min_Tot.Genes, fill = "band05")) +
  # geom_ribbon(aes(x= as.numeric(Step), ymin=1*2/Max_Tot.Genes, ymax=1*2/Min_Tot.Genes, fill = "band02")) +
  
  geom_line(aes(y = Prior.frac15, colour = "line15", group = 1), size = 1.3) +
  geom_point(aes(y = Prior.frac15, colour = "line15"), size=1) +
  geom_line(aes(y = Prior.frac10, colour = "line10", group = 1), size = 1.3) +
  geom_point(aes(y = Prior.frac10, colour = "line10"), size=1) +
  geom_line(aes(y = Prior.frac05, colour = "line05", group = 1), size = 1.3) +
  geom_point(aes(y = Prior.frac05, colour = "line05"), size=1) +
  geom_line(aes(y = Prior.frac02, colour = "line02", group = 1), size = 1.3) +
  geom_point(aes(y = Prior.frac02, colour = "line02"), size=1) +
  geom_line(aes(y = Missed.frac, colour = "missed", group = 1), size = 1.3) +
  geom_point(aes(y = Missed.frac, colour = "missed"), size=1) 

p <- p + theme_dark() + 
  ylab("Disease genes / Total dis. genes") +
  xlab("Max distance to primary reaction") +
  ylim(c(0,1)) +
  ggtitle("Correct disease gene prioritisation") +
  scale_color_manual(name = "Prioritised Genes",
                     labels = c("In Top 2", "In Top 5", "In Top 10", "In Top 15", "Missed"),
                     values=c(my_greens,'RED')) +
  scale_fill_manual(name = "Prioritised if random Gene",
                    labels = c("In Top 2", "In Top 5", "In Top 10", "In Top 15", "Missed"),
                    values=c(my_greens)) +
  guides(colour = guide_legend(order = 1, reverse = TRUE),
         fill = guide_legend(order = 2, reverse = TRUE))
# colour = adjustcolor(my_greens[4],alpha.f=0.5), fill = my_greens[4], alpha="0.5"
ggsave(paste0(outdir_name,"/",train_val,"_Ranks_And_Missed_Single_Par_Comb_",sub_name,".png"), plot = p,
       width = 300, height = 200, dpi=resolution, units = "mm")




#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Combination plot of correctly prioritised and missed genes --------------
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#####
# Disease gene prioritization per parameter, separated by training/validation runs/patients
p <- ggplot(DT_per_parameter_tra_val_sum, aes(x = Step)) +
  geom_line(data = DT_per_parameter_tra_val_sum[Validation == TRUE, ], aes(y = Av.prior.frac10, colour = "Validation", group = "Validation")) +
  geom_errorbar(data = DT_per_parameter_tra_val_sum[Validation == TRUE, ], 
                aes(ymax = Av.prior.frac10 + Sd.prior.frac10, 
                    ymin = Av.prior.frac10 - Sd.prior.frac10, 
                    colour = "Validation"), 
                width = 0.4) +
  geom_line(data = DT_per_parameter_tra_val_sum[Validation == FALSE, ], aes(y = Av.prior.frac10, colour = "Training", group = "Training")) +
  geom_errorbar(data = DT_per_parameter_tra_val_sum[Validation == FALSE, ], 
                aes(ymax = Av.prior.frac10 + Sd.prior.frac10, 
                    ymin = Av.prior.frac10 - Sd.prior.frac10, 
                    colour = "Training"), 
                width = 0.4) +
  geom_line(data = DT_per_parameter, aes(y = Missed.frac, colour = "Missed", group = 1)) +
  geom_point(data = DT_per_parameter,aes(y = Missed.frac, colour = "Missed"), size=0.5) +
  facet_grid(Max_rxn ~ Z_threshold, labeller = labeller(Max_rxn = Rxn_labs, Z_threshold = Thresh_labs)) +
  theme_light() + 
  ylab("Disease genes / Total dis. genes") +
  xlab("Max distance to primary reaction") +
  ggtitle("Correct disease gene prioritisation") +
  ylim(0,1)
for(i in c(1:5)){
  p <- p + geom_point(data = DT_per_parameter_tra_val_sum[Validation == FALSE & as.vector(DT_per_parameter_tra_val_sum[Validation == FALSE, "best_order_top10"]==i),],
                      aes(x=Step, y=Av.prior.frac10, shape = "Best"), 
                      shape = "*", 
                      size=8)
  p <- p + geom_text(data = DT_per_parameter_tra_val_sum[Validation == FALSE & as.vector(DT_per_parameter_tra_val_sum[Validation == FALSE, "best_order_top10"]==i),],
                     aes(x=Step, y=Av.prior.frac10, label = signif(Av.prior.frac10, digits = digit_significance), group = best_order_top10),
                     size=3,
                     colour = my_sig_palette[i],
                     position = position_dodge(width = 2),
                     vjust = -0.8)
}

z <- pretty_plot(p, theme = "light")
ggsave(paste0(outdir_name,"/CorrectPrior_top10_tra_val_",sub_name,".svg"), plot = z,
       width = 300, height = 200, units = "mm")

#####
# Disease gene prioritization per parameter, all runs/patients averaged to one value
p <- ggplot(DT_per_parameter, aes(x = Step)) +
  geom_line(aes(y = Missed.frac, colour = "Missed", group = 1)) +
  geom_point(aes(y = Missed.frac, colour = "Missed"), size=0.5) +
  geom_line(aes(y = Prior.frac10, colour = "line10", group = 1)) +
  geom_point(aes(y = Prior.frac10, colour = "line10"), size=0.5) +
  facet_grid(Max_rxn ~ Z_threshold, labeller = labeller(Max_rxn = Rxn_labs, Z_threshold = Thresh_labs)) +
  theme_light() + 
  ylab("Disease genes / Total dis. genes") +
  xlab("Max distance to primary reaction") +
  ggtitle("Correct disease gene prioritisation") +
  scale_color_manual(name = "Prioritised Genes",
                     labels = c("In Top 10", "Missed"),
                     values=c(my_greens[1],'RED')) +
  ylim(0,1)
for(i in c(1:5)){
  p <- p + geom_point(data = DT_per_parameter[as.vector(DT_per_parameter[, "best_order_top10"]==i),],
                      aes(x=Step, y=Prior.frac10, shape = "Best"), 
                      shape = "*", 
                      size=8)
  p <- p + geom_text(data = DT_per_parameter[as.vector(DT_per_parameter[, "best_order_top10"]==i),],
                     aes(x=Step, y=Prior.frac10, label = signif(Prior.frac10, digits = digit_significance), group = best_order_top10),
                     size=3,
                     colour = my_sig_palette[i],
                     position = position_dodge(width = 2),
                     vjust = -0.8)
}

z <- pretty_plot(p, theme = "light")
ggsave(paste0(outdir_name,"/CorrectPrior_top10_all_",sub_name,".svg"), plot = z,
       width = 300, height = 200, units = "mm")



# p <- ggplot(DT_per_parameter, aes(x = Step)) +
#   geom_line(aes(y = Prior.frac05, colour = "line05", group = 1)) +
#   geom_point(aes(y = Prior.frac05, colour = "line05"), size=0.5) +
#   geom_line(aes(y = Prior.frac10, colour = "line10", group = 1)) +
#   geom_point(aes(y = Prior.frac10, colour = "line10"), size=0.5) +
#   geom_line(aes(y = Missed.frac, group = 1, linetype="Training"), colour = my_reds[1]) +
#   geom_point(aes(y = Missed.frac), size=0.5, colour = my_reds[1]) 
# # geom_hline(data = DT_per_parameter[Set == "training"], aes(yintercept = max(Prior.frac05), size = "Training", colour="line05t"), alpha = 0.8, linetype = "dashed") +
# # geom_hline(data = DT_per_parameter[Set == "validation"], aes(yintercept = max(Prior.frac05), size = "Validation", colour="line05v"), alpha = 0.8, linetype = "dashed")
# p <- p + facet_grid(Max_rxn ~ Z_threshold, labeller = labeller(Max_rxn = Rxn_labs, Z_threshold = Thresh_labs)) 
# p <- p + theme_light() + 
#   ylab("Ratio [rank disease gene] / [total genes]") +
#   # ylab(expression(paste("Ratio of: ", frac(`disease genes`,`total disease genes`)))) +
#   xlab("Maximum distance to primary reaction") +
#   ylim(0, 1) +
#   # ggtitle("Performance of disease gene prioritization") +
#   scale_color_manual(name = NULL,
#                      labels = c("Correctly prioritized genes (top 5)","Correctly prioritized genes (top 10)"),
#                      values= c(my_greens[3], my_greens[2]),
#                      guide = guide_legend(override.aes = list(linetype = "solid",
#                                                               shape = 16), 
#                                           order = 1)
#   ) +
#   scale_linetype_manual(name = NULL,
#                         labels = "Missed genes",
#                         values = "solid", 
#                         guide = guide_legend(override.aes = list(colour = my_reds[3]), 
#                                              order = 2)
#   )
# pp <- pretty_plot(p, theme = "light")
# ggsave(paste0(outdir_name,"/",train_val,"_Ranks_And_Missed_",sub_name,".png"), plot = pp,
#        width = 300, height = 200, dpi=resolution, units = "mm")


# Some example plot for making multiple legends out of 1 aesthetic
# ggplot(data=dfr, mapping=aes(x=id, y=value)) +
#   geom_line(mapping=aes(colour=group), show_guide=TRUE) +
#   geom_hline(
#     mapping=aes(yintercept=c(-1,1)*qnorm(0.95), fill="95% CI"),
#     color="orange"
#   ) +
#   geom_hline(
#     mapping=aes(yintercept=c(-1,1)*qnorm(0.99), fill="99% CI"),
#     color="darkred"
#   ) +
#   scale_color_hue("Group", guide=guide_legend(order=1)) +
#   scale_fill_manual("CI horizontal line", values=rep(1,4),
#                     guide=guide_legend(
#                       override.aes = list(colour=c("orange", "darkred")),
#                       order=2
#                     ),
#                     labels=c("CI of 95%", "CI of 99%")
#   )



# ##### Average rank non-missed genes + total missed genes ------------------
# p <- ggplot(DT_per_parameter, aes(label=Av_non_missed)) +
#   geom_line(aes(x = Step, y = Av_non_missed, colour = "Average", group = 1), size = 1.3) +
#   geom_point(aes(x = Step, y = Av_non_missed, colour = "Average"), size=0.5) +
#   # geom_line(aes(x = Step, y = Missed.frac*15, colour = "Frac.Missed", group = 1), size = 1.3) +
#   # geom_point(aes(x = Step, y = Missed.frac*15, colour = "Frac.Missed"), size=0.5) +
#   # scale_y_continuous(sec.axis = sec_axis(~./15, name = "Frac. genes Missed", breaks = c(0, 0.5, 1))) +
#   theme_dark() +
#   scale_color_manual(name = "",
#                      labels = c("Av. rank of non-missed genes", "Frac. Missed"),
#                      values=c('Black', my_greens[3]))
# p <- p + facet_grid(Max_rxn ~ Z_threshold, labeller = labeller(Max_rxn = Rxn_labs, Z_threshold = Thresh_labs))
# # geom_point(data = DT_per_parameter[best_av50==TRUE, ],aes(x=Step, y=Av_top50),shape = "*", size=8, show.legend = FALSE, colour = "black")+
# p <- p + geom_point(data = DT_per_parameter[as.vector(DT_per_parameter[,"best_order_NM"]==1),],aes(x=Step, y=Av_non_missed, shape = "Best"), size=8, ) +
#   scale_shape_manual(name = "",
#                      labels = "Best param.\nranks (missed \nfrac. <0.5)",
#                      values = 42)
# # Annotate the 5 best scoring parameter combinations (determined ~l.280)
# for(i in c(1:5)){
#   p <- p + geom_point(data = DT_per_parameter[as.vector(DT_per_parameter[,"best_order_NM"]==i),],aes(x=Step, y=Av_non_missed, shape = "Best"), shape = "*", size=8, colour = my_sig_palette[i])
#   p <- p + geom_text(data = DT_per_parameter[as.vector(DT_per_parameter[,"best_order_NM"]==i),],
#                      aes(x=Step, y=Av_non_missed, label = signif(Av_non_missed, digits = digit_significance), group = best_order_NM),
#                      size=3,
#                      # show.legend = FALSE,
#                      colour = my_sig_palette[i],
#                      position = position_dodge(width = 2),
#                      vjust = -0.5)
# }
# y_val <- ggplot_build(p)$layout$panel_scales_y[[1]]$range$range[2]
# p <- p +
#   ylab("Average non-missed disease gene rank") +
#   xlab("Max. distance to primary reaction") +
#   ggtitle("Av. gene rank for non-missed genes") +
#   # geom_hline(yintercept=y_val/2, linetype="dashed", color = "salmon") +
#   guides(shape = guide_legend(order = 1),
#          colour = guide_legend(order = 2))
# pp <- pretty_plot(p, theme = "dark")
# ggsave(paste0(outdir_name,"/Best_Ranks_With_>0.5_Non_Missed_",sub_name,".png"), plot = pp,
#        width = 300, height = 200, dpi=resolution, units = "mm")


##### all ranks, average rank and average standard deviation --------------
# DT_no_missed <- tmpDT[P.value < 1,]
# DT_missed <- tmpDT[P.value == 1,]
# p <- ggplot(DT_no_missed, aes(x = Step)) +
#   # geom_point(aes(y = Position), position=position_dodge(width = 5.5)) +
#   geom_jitter(aes(y = Position, colour = "Per-patient rank"), size = 0.005) +
#   geom_jitter(data =DT_missed, aes(x = Step, y = Position, fill = "Missed genes"), colour = "blue", size = 0.005, alpha = 0.2) +
#   geom_line(data = DT_per_parameter, aes(x = Step, y = Av_non_missed, group = 1, colour = "Average rank")) +
#   geom_errorbar(data = DT_per_parameter, aes(x = Step, 
#                                              ymax = Av_non_missed + Av_Sd_excl_miss, 
#                                              ymin = Av_non_missed - Av_Sd_excl_miss,
#                                              colour = "Average rank"), 
#                 # colour = "cornflowerblue",
#                 width = 0.5) +
#   ylim(c(0,40)) +
#   geom_hline(aes(yintercept = 5), colour = "darksalmon", linetype="dashed") +
#   facet_grid(Max_rxn ~ Z_threshold, labeller = labeller(Max_rxn = Rxn_labs, Z_threshold = Thresh_labs)) +
#   theme_dark() + 
#   ylab("Disease gene rank") +
#   xlab("Max distance to primary reaction") +
#   ggtitle("Stability of method") +
# geom_point(data = DT_per_parameter[DT_per_parameter$best_order_top05 ==1, ], aes(x=Step, y=35, shape = "Best"), size=8) + 
#   scale_shape_manual(name = "",
#                      labels = "Ratio dis.genes\nin top 5",
#                      values = 42) 
# for(i in c(1:5)){
#   p <- p + geom_point(data = DT_per_parameter[as.vector(DT_per_parameter[,"best_order_top05"]==i),],aes(x=Step, y=35, shape = "Best"), shape = "*", size=8)
#   p <- p + geom_text(data = DT_per_parameter[as.vector(DT_per_parameter[,"best_order_top05"]==i),],
#                      aes(x=Step, y=35, label = signif(Prior.frac05, digits = digit_significance), group = best_order_top05),
#                      size=3,
#                      # show.legend = FALSE,
#                      colour = my_sig_palette[6-i],
#                      position = position_dodge(width = 2),
#                      vjust = -0.5)
# }
# p <- p + scale_color_manual(name = "Non-missed genes",
#                             labels = c("Average rank", "Per-patient rank"),
#                             values=c("red","black")) + 
#   scale_fill_manual(name = "Missed genes",
#                             labels = c("Per-patient rank"),
#                             values=c("blue"))
# p <- pretty_plot(p, secondary_y_axis = FALSE, theme = "dark")
# ggsave(paste0(outdir_name,"/",train_val,"_Average_Patient_And_Ranks_And_Missed_",sub_name,".png"), plot = p,
#        width = 300, height = 200, dpi=resolution, units = "mm")


# ##### all ranks, average rank and average standard deviation --------------
# # p <- ggplot(tmpDT[Missed == FALSE], aes(x = Step)) +
# p <- ggplot(DT_per_parameter, aes(x = Step)) +
#   # geom_point(aes(y = Position), position=position_dodge(width = 5.5)) +
#   # geom_jitter(data = tmpDT, aes(x = Step, y = Position, colour = Missed), size = 0.005, alpha = 0.2) +
#   # geom_jitter(data = DT, aes(y = Position, colour = Missed), size = 0.005, alpha = 0.2) +
#   # geom_jitter(aes(y = Position, colour = "Per-patient rank"), size = 0.005) +
#   # geom_jitter(data =DT_test[Missed == TRUE], aes(x = Step, y = Position, fill = "Missed genes"), colour = "blue", size = 0.005, alpha = 0.2) +
#   geom_line(aes(x = Step, y = Av_non_missed, group = 1, colour = "Average rank")) +
#   geom_errorbar(aes(x = Step, 
#                                              ymax = Av_non_missed + Sd_non_missed, 
#                                              ymin = Av_non_missed - Sd_non_missed,
#                                              colour = "Average rank"), 
#                 # colour = "cornflowerblue",
#                 width = 0.5) +
#   ylim(c(-5,35)) +
#   geom_hline(aes(yintercept = 5), colour = "darksalmon", linetype="dashed") +
#   facet_grid(Max_rxn ~ Z_threshold, labeller = labeller(Max_rxn = Rxn_labs, Z_threshold = Thresh_labs)) +
#   theme_dark() + 
#   ylab("Disease gene rank") +
#   xlab("Max distance to primary reaction") +
#   ggtitle("Gene prioritization performance") +
#   geom_point(data = DT_per_parameter[DT_per_parameter$best_order_top05 ==1, ], aes(x=Step, y=35, shape = "Best"), size=8) + 
#   scale_shape_manual(name = "",
#                      labels = "Best ratios dis.genes\nin top 5",
#                      values = 42) 
# for(i in c(1:5)){
#   # p <- p + geom_point(data = DT_per_parameter[as.vector(DT_per_parameter[,"best_order_top05_training"]==i),],aes(x=Step, y=35, shape = "Best"), shape = "*", size=8)
#   # p <- p + geom_text(data = DT_per_parameter[as.vector(DT_per_parameter[,"best_order_top05_training"]==i),],
#   #                    aes(x=Step, y=35, label = signif(Prior.frac05_training, digits = digit_significance), group = best_order_top05_training),
#   #                    size=3,
#   #                    # show.legend = FALSE,
#   #                    colour = "black",
#   #                    position = position_dodge(width = 2),
#   #                    vjust = -0.5)
#   # p <- p + geom_point(data = DT_per_parameter[as.vector(DT_per_parameter[,"best_order_top05_validation"]==i),],aes(x=Step, y=35, shape = "Best"), shape = "*", size=8)
#   # p <- p + geom_text(data = DT_per_parameter[as.vector(DT_per_parameter[,"best_order_top05_validation"]==i),],
#   #                    aes(x=Step, y=35, label = signif(Prior.frac05_validation, digits = digit_significance), group = best_order_top05_validation),
#   #                    size=3,
#   #                    # show.legend = FALSE,
#   #                    colour = "black",
#   #                    position = position_dodge(width = 2),
#   #                    vjust = -0.5)
#   p <- p + geom_point(data = DT_per_parameter[as.vector(DT_per_parameter[,"best_order_top05"]==i),],aes(x=Step, y=25, shape = "Best"), shape = "*", size=8)
#   p <- p + geom_text(data = DT_per_parameter[as.vector(DT_per_parameter[,"best_order_top05"]==i),],
#                      aes(x=Step, y=25, label = signif(Prior.frac05, digits = digit_significance), group = best_order_top05),
#                      size=3,
#                      # show.legend = FALSE,
#                      colour = my_sig_palette[6-i],
#                      position = position_dodge(width = 2),
#                      vjust = -0.5)
# }
# p <- p + scale_color_manual(name = "Gene ranks",
#                             labels = c("Average non-missed", "Per patient", "Missed per patient"),
#                             values=c("red","black","blue")) 
# p <- p + guides(shape = guide_legend(override.aes = list(size = 5)),
#                 color = guide_legend(override.aes = list(linetype=c(1,NA,NA), 
#                                                          shape=c(NA,16,16), 
#                                                          size = c(0.5,2,2),
#                                                          alpha = c(NA,1,1))))
# pp <- pretty_plot(p, secondary_y_axis = FALSE, theme = "dark")
# # ggsave(paste0(outdir_name,"/",train_val,"_Average_Patient_And_Ranks_And_Missed_",sub_name,".png"), plot = pp,
# #        width = 300, height = 200, dpi=resolution, units = "mm")
# ggsave(paste0(outdir_name,"/test_",sub_name,".png"), plot = p,
#        width = 300, height = 200, dpi=resolution, units = "mm")


##### Boxplot of a single panel ---------------------------------------
# nm <- 0
# DT[,Patient_follow_number := 0]
# DT[, c("Min_Rank_Patient","Max_Rank_Patient") := list(min(Position),max(Position)), by = .(Step, Z_threshold, Max_rxn, PatientID)]
# for(i in unique(DT$PatientID)){
#   nm <- nm + 1
#   DT[PatientID == i, Patient_follow_number := nm]
# }

# # Relative ranks
# # p <- ggplot(DT[Z_threshold == "-3, 3" & Max_rxn == 10 & P.value != 1 & Step == 5], aes(x = Patient_follow_number, y = Rank.frac)) +
# p <- ggplot(DT[Z_threshold == "-3, 3" & Max_rxn == 10 & P.value != 1 & Step == 5], aes(x = PatientID, y = Position/201)) +
#   geom_boxplot(aes(fill = Patient_follow_number, group = PatientID)) + theme(legend.position = "none") +
#   geom_hline(yintercept = 11/201, linetype = "dashed")
#   # geom_ribbon(aes(x = Patient_follow_number, ymin=11/201,ymax=11/Min_Rank_Patient),alpha=0.2, fill = "yellow") +
#   # geom_ribbon(aes(x = Patient_follow_number, ymin=11/201,ymax=0.6),alpha=0.2, fill = "red") +
#   # geom_ribbon(aes(x = Patient_follow_number, ymin = 0, ymax=6/201),alpha=0.2, fill = "green") +
#   # theme(axis.text.x = element_text(size = 6, angle = 60, hjust = 1)) +
#   # ylim(0,1) +
#   # scale_x_continuous(labels = unique(DT[Z_threshold == "-3, 3" & Max_rxn == 10 & P.value != 1 & Step == 5, PatientID]),
#   #                    breaks = unique(DT[Z_threshold == "-3, 3" & Max_rxn == 10 & P.value != 1 & Step == 5,Patient_follow_number]),
#   #                    name = "Patient ID")
# ggsave(paste0(outdir_name,"/",train_val,"Relative_Rank_",sub_name,".png"), plot = p,
#        width = 400, height = 260, dpi=resolution, units = "mm")

# Absolute ranks
tmpDT <- DT[Z_threshold == "-3, 3" & Max_rxn == 15 & P.value != 1 & Step == 4 & Include == TRUE]
tmpDT[ , P.value := as.double(P.value)]
tmpDT <- tmpDT[order(tmpDT[, Gene]),]
tmpDT[, PatientID := paste(Gene, PatientID, sep = "^")]
disease_genes <- NULL
for(i in unique(tmpDT[,PatientID])){
  tmp_gene <- unique(tmpDT[PatientID == i, Gene])
  cat(tmp_gene, i, "\n")
  disease_genes <- c(disease_genes, tmp_gene)
}


p <- ggplot(tmpDT, aes(x = PatientID, y = Position)) +
  # geom_boxplot(aes(fill = P.colour, group = PatientID)) + 
  # geom_boxplot(aes(fill = DBS, group = PatientID)) + 
  # geom_boxplot(aes(fill = log(P.value), group = PatientID)) + 
  # geom_boxplot(aes(fill = P.value, group = PatientID)) + 
  # scale_fill_gradient(high = "#c92222", low = "#22c939") +
  geom_boxplot(aes(group = PatientID)) +
  # theme(legend.position = "none") +
  # geom_boxplot(aes(fill = Patient_follow_number, group = PatientID)) + theme(legend.position = "none") +
  geom_hline(yintercept = 10, linetype = "dashed") +
  # add DBS to plot
  # geom_point(aes(x = Patient_follow_number, y = DBS*5), shape = 8, colour = 'red', alpha = 0.8) +
  theme(axis.text.x = element_text(size = 6, angle = 60, hjust = 1))+
  ylab("Disease gene rank") +
  scale_x_discrete(labels = disease_genes,
                   name = "Disease genes of individual patients")
ggsave(paste0(outdir_name,"/1Par_Absolute_Rank_",sub_name,".svg"), plot = p,
       width = 400, height = 150, units = "mm")


# # visualise Sd_rank & Av_rank ~ DBS per patient
# DT_per_patient[, DBS := factor(DBS)]
# p <- ggplot(DT_per_patient, aes(DBS, Sd_rank_excl_miss)) + 
#   geom_jitter(colour = "darkgreen", alpha = 0.8) +
#   geom_boxplot(alpha = 0.8)
# ggsave(paste0(outdir_name,"/Sd_Abs_Rank_(NoMiss)_vs_DBS.png"), plot = p,
#        width = 400, height = 260, dpi=resolution, units = "mm")
# 
# p <- ggplot(DT_per_patient, aes(DBS, Av_rank_excl_miss)) + 
#   geom_jitter(colour = "darkred", alpha = 0.8) +
#   geom_boxplot(alpha = 0.8)
# ggsave(paste0(outdir_name,"/Av_Abs_Rank_(NoMiss)_vs_DBS.png"), plot = p,
#        width = 400, height = 260, dpi=resolution, units = "mm")
# 
# # visualise Sd_rank & Av_rank ~ Dataset/Run
# DT_per_patient[, Dataset := factor(Dataset)]
# ggplot(DT_per_patient, aes(Dataset, Sd_rank_excl_miss)) + 
#   geom_jitter(colour = "darkgreen", alpha = 0.8) +
#   geom_boxplot(alpha = 0.8) +
#   theme(axis.text.x = element_text(size = 6, angle = 45, hjust = 1))
# ggplot(DT_per_patient, aes(Dataset, Av_rank_excl_miss)) + 
#   geom_jitter(colour = "darkred", alpha = 0.8) +
#   geom_boxplot(aes(fill = Set), alpha = 0.8) +
#   theme(axis.text.x = element_text(size = 6, angle = 45, hjust = 1))
# 
# # visualise Sd_rank & Av_rank ~ metabolite_set_size per patient (unbinned)
# ggplot(DT_per_patient,aes(Mets_in_set, Sd_rank_excl_miss)) + 
#   geom_point() 
# ggplot(DT_per_patient,aes(Mets_in_set, Av_rank_excl_miss)) + 
#   geom_point() 
# 
# # visualise Sd_rank & Av_rank ~ metabolite_set_size per patient (binned)
# DT_per_patient[, Mets_in_set := as.numeric(Mets_in_set)]
# test <- copy(DT_per_patient)
# test[,Bin := "0"]
# test[Mets_in_set > 0, Bin := "0 < x ≤ 5"]
# test[Mets_in_set > 5, Bin := "5 < x ≤ 10"]
# test[Mets_in_set > 10, Bin := "10 < x ≤ 20"]
# test[Mets_in_set > 20, Bin := "20 < x ≤ 40"]
# test[Mets_in_set > 40, Bin := "40 < x ≤ 80"]
# test[Mets_in_set > 80, Bin := "80 < x ≤ 140"]
# test[Mets_in_set > 140, Bin := "140 < x ≤ 200"]
# test[Mets_in_set > 200, Bin := "200 < x ≤ 300"]
# test[Mets_in_set > 300, Bin := "300 < x"]
# test[,Bin := factor(Bin)]
# levels(test$Bin) <- c("0", "0 < x ≤ 5", "5 < x ≤ 10", "10 < x ≤ 20", "20 < x ≤ 40", "40 < x ≤ 80", "80 < x ≤ 140",
#                       "140 < x ≤ 200", "200 < x ≤ 300", "300 < x")
# 
# p <- ggplot(test, aes(x = Bin, y = Sd_rank_excl_miss)) +
#   geom_jitter(colour = "darkgreen", alpha = 0.8) +
#   geom_boxplot(alpha = 0.8)
# ggsave(paste0(outdir_name,"/Sd_Abs_Rank_(NoMiss)_vs_Met_Set_Size.png"), plot = p,
#        width = 400, height = 260, dpi=resolution, units = "mm")
# p <- ggplot(test, aes(x = Bin, y = Av_rank_excl_miss)) +
#   geom_jitter(colour = "darkred", alpha = 0.8) +
#   geom_boxplot(alpha = 0.8)
# ggsave(paste0(outdir_name,"/Av_Abs_Rank_(NoMiss)_vs_Met_Set_Size.png"), plot = p,
#        width = 400, height = 260, dpi=resolution, units = "mm")


#####
# difference of #DBS on the rankings

# Determine which genes are present 1. in multiple patients and 2. with different number of DBS
# geneDBS <- unique(DT[, paste(Gene, DBS, sep = ";")])
# geneDBS <- names(which(table(unlist(lapply(geneDBS, function(x) strsplit(x, split = ";")[[1]][1]))) > 1))
# 
# 
# p <-ggplot(DT_per_patient[Gene %in% geneDBS & Missed == 0, ], aes(x = as.factor(DBS), y = Av_rank)) +
#   geom_jitter(aes(colour = Dataset), alpha = 0.8) +
#   geom_boxplot(aes(colour = Dataset), alpha = 0.8) + 
#   # theme(legend.position="bottom") +
#   facet_wrap(. ~ Gene, scales='free_x')
# ggsave(paste0(outdir_name,"/Genes_With_Multiple_DBS.png"), plot = p,
#        width = 400, height = 260, dpi=resolution, units = "mm")
# 
# 
# p2 <- ggplot(DT_per_patient[Gene %in% geneDBS, ], aes(x = as.factor(DBS), y = Sd_rank)) +
#   geom_jitter(aes(colour = PatientID), alpha = 0.8) +
#   geom_boxplot(aes(colour = PatientID), alpha = 0.8) +
#   facet_wrap(. ~ Gene, scales='free_x')
# 
# 
# p1 <-ggplot(DT_per_patient[Gene=="CPT1A",], aes(x = as.factor(DBS), y = Av_rank)) +
#   geom_jitter(aes(colour = PatientID), alpha = 0.8) +
#   geom_boxplot(aes(colour = PatientID), alpha = 0.8) + 
#   theme(legend.position="bottom")
# p2 <- ggplot(DT_per_patient[Gene=="CPT1A",], aes(x = as.factor(DBS), y = Sd_rank)) +
#   geom_jitter(aes(colour = PatientID), alpha = 0.8) +
#   geom_boxplot(aes(colour = PatientID), alpha = 0.8)
# 
# g_legend<-function(a.gplot){
#   tmp <- ggplot_gtable(ggplot_build(a.gplot))
#   leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
#   legend <- tmp$grobs[[leg]]
#   return(legend)}
# mylegend<-g_legend(p1)
# p3 <- gridExtra::grid.arrange(arrangeGrob(p1 + theme(legend.position="none"),
#                          p2 + theme(legend.position="none"),
#                          nrow=1),
#              mylegend, nrow=2,heights=c(10, 1))
# # gridExtra::grid.arrange(p1, p2, ncol=2)
# ggsave(paste0(outdir_name,"/CPT1A_Av_Abs_Rank_And_Sd_vs_DBS.png"), plot = p3,
#        width = 400, height = 260, dpi=resolution, units = "mm")



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Per patient plots -------------------------------------------------------
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# ##### Average rank per patient --------------------------------------------
# p <- ggplot(data = DT_per_patient, aes(x = Step, y = Av_rank)) +
#   geom_bar(position = "dodge", stat="identity", aes(fill = PatientID)) +
#   geom_line(data = DT_per_parameter, aes(x = Step, y = Missed/0.25, group = 1, colour = "Genes Missed")) +
#   scale_y_continuous(sec.axis = sec_axis(~.*0.25, name = "Tot.Genes Missed")) +
#   ylab("Average Rank") +
#   xlab("Max. distance to primary reaction") +
#   ggtitle("Disease gene rank per patient & parameter combination") +
#   labs(colour = "", fill = "Patients") +
#   guides(colour = guide_legend(order = 1),
#          fill = guide_legend(order = 2)) +
#   facet_grid(Max_rxn ~ Z_threshold, labeller = labeller(Max_rxn = Rxn_labs, Z_threshold = Thresh_labs))
# p <- pretty_plot(p)
# ggsave(paste0(outdir_name,"/_Av_Rank_Per_Patient_",sub_name,".png"), plot = p,
#        width = 1000, height = 200, dpi=resolution, units = "mm")
# 
# 
# ##### Standard deviation of average rank per patient ----------------------
# p <- ggplot(data = DT_per_patient, aes(x = Step, y = Sd_rank)) +
#   geom_bar(position = "dodge", stat="identity", aes(fill = PatientID)) +
#   # geom_line(data = DT_per_parameter, aes(x = Step, y = Missed/4, group = 1, colour = "Genes Missed")) +
#   # scale_y_continuous(sec.axis = sec_axis(~.*4, name = "Tot.Genes Missed")) +
#   # scale_y_continuous(sec.axis = sec_axis(~.*4, name = "Tot.Genes Missed", breaks = seq(0, 50, 10))) +
#   ylab("Sd disease gene rank") +
#   xlab("Max. distance to primary reaction") +
#   ggtitle("Disease gene St.dev. per patient & parameter combination") +
#   labs(colour = "", fill = "Patients") +
#   guides(colour = guide_legend(order = 1),
#          fill = guide_legend(order = 2)) +
#   facet_grid(Max_rxn ~ Z_threshold, labeller = labeller(Max_rxn = Rxn_labs, Z_threshold = Thresh_labs))
# p <- pretty_plot(p)
# ggsave(paste0(outdir_name,"/",train_val,"_Sd_Rank_Per_Patient_",sub_name,".png"), plot = p,
#        width = 300, height = 200, dpi=resolution, units = "mm")
# 
# 
# ##### Average relative rank per patient -----------------------------------
# p <- ggplot(data = DT_per_patient, aes(x = Step, y = Av_Rel_rank)) +
#   geom_bar(position = "dodge", stat="identity", aes(fill = PatientID)) +
#   # geom_line(data = DT_per_parameter, aes(x = Step, y = Missed/50, group = 1, colour = "Genes Missed")) +
#   # scale_y_continuous(sec.axis = sec_axis(~.*50, name = "Tot.Genes Missed")) +
#   ylab("Average Relative Rank (rank/Tot. genes") +
#   xlab("Max. distance to primary reaction") +
#   ggtitle("Disease gene rel.rank per patient & parameter combination") +
#   labs(colour = "", fill = "Patients") +
#   guides(colour = guide_legend(order = 1),
#          fill = guide_legend(order = 2)) +
#   facet_grid(Max_rxn ~ Z_threshold, labeller = labeller(Max_rxn = Rxn_labs, Z_threshold = Thresh_labs))
# p <- pretty_plot(p)
# ggsave(paste0(outdir_name,"/",train_val,"_Av_Rel_Rank_Per_Patient_",sub_name,".png"), plot = p,
#        width = 300, height = 200, dpi=resolution, units = "mm")
# 
# 
# ##### Standard deviation of average relative rank per patient -------------
# p <- ggplot(data = DT_per_patient, aes(x = Step, y = Sd_Rel_rank)) +
#   geom_bar(position = "dodge", stat="identity", aes(fill = PatientID)) +
#   # geom_line(data = DT_per_parameter, aes(x = Step, y = Missed/800, group = 1, colour = "Genes Missed")) +
#   # scale_y_continuous(sec.axis = sec_axis(~.*800, name = "Tot.Genes Missed")) +
#   ylab("Sd disease gene relative Rank") +
#   xlab("Max. distance to primary reaction") +
#   ggtitle("Disease gene St.dev. per patient & parameter combination") +
#   labs(colour = "", fill = "Patients") +
#   guides(colour = guide_legend(order = 1),
#          fill = guide_legend(order = 2)) +
#   facet_grid(Max_rxn ~ Z_threshold, labeller = labeller(Max_rxn = Rxn_labs, Z_threshold = Thresh_labs))
# p <- pretty_plot(p)
# ggsave(paste0(outdir_name,"/",train_val,"_Sd_Rel_Rank_Per_Patient_",sub_name,".png"), plot = p,
#        width = 300, height = 200, dpi=resolution, units = "mm")
# 
# 
# ##### Average relative rank per patient, reversed (1 = rank 1) ------------
# p <- ggplot(data = DT_per_patient, aes(x = Step, y = Av_Rev_Rel_rank)) +
#   geom_bar(position = "dodge", stat="identity", aes(fill = PatientID)) +
#   # geom_line(data = DT_per_parameter, aes(x = Step, y = Missed/50, group = 1, colour = "Genes Missed")) +
#   # scale_y_continuous(sec.axis = sec_axis(~.*50, name = "Tot.Genes Missed")) +
#   ylab("Average reverse Relative Rank (1-((rank-1)/(Tot. genes in set-1))") +
#   xlab("Max. distance to primary reaction") +
#   ggtitle("Disease gene rel.rank per patient & parameter combination") +
#   labs(colour = "", fill = "Patients") +
#   guides(colour = guide_legend(order = 1),
#          fill = guide_legend(order = 2)) +
#   facet_grid(Max_rxn ~ Z_threshold, labeller = labeller(Max_rxn = Rxn_labs, Z_threshold = Thresh_labs))
# p <- pretty_plot(p)
# ggsave(paste0(outdir_name,"/",train_val,"_Av_Rev_Rel_Rank_Per_Patient_",sub_name,".png"), plot = p,
#        width = 300, height = 200, dpi=resolution, units = "mm")
# 
# 
# ##### St. dev. of average relative (reversed) rank per patient ------------
# p <- ggplot(data = DT_per_patient, aes(x = Step, y = Sd_Rev_Rel_rank)) +
#   geom_bar(position = "dodge", stat="identity", aes(fill = PatientID)) +
#   # geom_line(data = DT_per_parameter, aes(x = Step, y = Missed/800, group = 1, colour = "Genes Missed")) +
#   # scale_y_continuous(sec.axis = sec_axis(~.*800, name = "Tot.Genes Missed")) +
#   ylab("Sd disease gene reverse relative Rank") +
#   xlab("Max. distance to primary reaction") +
#   ggtitle("Disease gene reverse St.dev. per patient & parameter combination") +
#   labs(colour = "", fill = "Patients") +
#   guides(colour = guide_legend(order = 1),
#          fill = guide_legend(order = 2)) +
#   facet_grid(Max_rxn ~ Z_threshold, labeller = labeller(Max_rxn = Rxn_labs, Z_threshold = Thresh_labs))
# p <- pretty_plot(p)
# ggsave(paste0(outdir_name,"/",train_val,"_Sd_Rev_Rel_Rank_Per_Patient_",sub_name,".png"), plot = p,
#        width = 300, height = 200, dpi=resolution, units = "mm")

