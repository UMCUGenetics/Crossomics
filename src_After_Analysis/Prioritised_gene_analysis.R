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
library(RColorBrewer)
library(ggplot2)
library(scales)
# library(ggforce) # zoom in specific parts of plot
library(grid) # manual grob adjustment
library(gtable) # manual grob adjustment



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Adjustable settings -----------------------------------------------------
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

##### Create general plot aesthetics --------------------------------------
Thresh_labs <- c("Min.-0.5; Max.1","Min.-1; Max.1.5","Min.-1.5; Max.2","Min.-3; Max.3","Min.-5; Max.5")
names(Thresh_labs) <- c("-0.5, 1", "-1, 1.5", "-1.5, 2", "-3, 3", "-5, 5")
Rxn_labs <- c("<= 8","<= 10","<= 12","<= 15","<= 17","<= 19")
names(Rxn_labs) <- c(8, 10, 12, 15, 17, 19)

# Colour scheme 
my_greens <- rev(brewer.pal(5, "Greens"))[c(2:5)]
my_sig_palette <- rev(brewer.pal(6, "RdYlGn"))[2:6]

# Image resolution
high_res <- 600
low_res <- 100
resolution <- low_res
digit_significance <- 3

##### Other ---------------------------------------------------------------
# Exclude patients / diseases
NoTrans <- FALSE
subset_patients <- FALSE

# Date of data
date <- "2019-09-03"



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Pre-processing ----------------------------------------------------------
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

##### Read data -----------------------------------------------------------
code_dir <- paste0(dirname(rstudioapi::getActiveDocumentContext()$path),"/../Results/")
DT <- data.table::as.data.table(readRDS(paste0(code_dir,date,"/MSEA_DT_compiled.RDS")))

##### Determine parameters ------------------------------------------------
# temp$size_f = factor(temp$size, levels=c('50%','100%','150%','200%'))
Z_thresholds <- levels(DT$Z_threshold)
max_rxns <- levels(DT$Max_rxn)
steps <- levels(DT$Step)
seeds <- unique(DT$Seed)

##### Exclude predetermined patients and determine output names -----------
if(subset_patients){
  patients_excluded <- c("P56","P57","P58","P59","P68")
  DT <- DT[!Patient %in% patients_excluded]
  sub_name <- "Sub_P"
} else {
  sub_name <- "All_P"
}

outdir_name <- paste0(dirname(rstudioapi::getActiveDocumentContext()$path),"/../Plots/",Sys.Date())
var_name <- ifelse(NoTrans, "No_Trans", "All_Genes")
if (!file.exists(outdir_name))dir.create(outdir_name, recursive = TRUE)




#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Functions ---------------------------------------------------------------
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

##### Manual standard grob addition ---------------------------------------
pretty_plot <- function(p, theme = "light", secondary_y_axis = TRUE){
  z <- ggplotGrob(p)
  right_lab_loc <- 5 + 2*length(Z_thresholds) - secondary_y_axis
  right_lab_bottom <- 6 + 2*length(max_rxns)
  top_lab_width <- 5 + 2*length(Z_thresholds) - 2
  
  # label right strip
  z <- gtable_add_cols(z, unit(z$widths[[7]], 'cm'), right_lab_loc)
  z <- gtable_add_grob(z, 
                       list(rectGrob(gp = gpar(col = NA, fill = ifelse(theme=="dark", "lightgray", gray(0.5)))),
                            textGrob("Extension stringency", rot = -90, gp = gpar(col = ifelse(theme=="dark", "black", "white")))),
                       # 8, 15, 18, 15, name = paste(runif(2)))
                       8, right_lab_loc+1, right_lab_bottom, right_lab_loc+1, name = paste(runif(2)))
  
  # label top strip
  z <- gtable_add_rows(z, unit(z$heights[[3]], 'cm'), 3)
  z <- gtable_add_grob(z, 
                       list(rectGrob(gp = gpar(col = NA, fill = ifelse(theme=="dark", "lightgray", gray(0.5)))),
                            textGrob("Biochemical stringency", gp = gpar(col = ifelse(theme=="dark", "black", "white")))),
                       # 4, 5, 4, 13, name = paste(runif(2)))
                       4, 5, 4, top_lab_width, name = paste(runif(2)))
  
  # margins
  z <- gtable_add_cols(z, unit(1/8, "line"), 7)
  z <- gtable_add_rows(z, unit(1/8, "line"), 3)
  
  return(z)
}




#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Create patient and parameter datatables ---------------------------------
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

##### Modify initial datatable --------------------------------------------
DT[,c("Prioritised15","Prioritised10","Prioritised05","Prioritised02") := 
     list(Position <= 15, Position <= 10, Position <= 5, Position <= 2)]
# DT$Prioritised15 <- DT$Position <= 15
# DT$Prioritised10 <- DT$Position <= 10
# DT$Prioritised05 <- DT$Position <= 5
# DT$Prioritised02 <- DT$Position <= 2
# DT$Prioritised50 <- DT$Position <= 50
DT[, PatientID := do.call(paste0,.SD), .SDcols=c("Patient","Gene")]

# DT_noTrans <- DT
# DT_noTrans <- DT_noTrans[Transporter==FALSE,,]

if(NoTrans) tmpDT <- DT[Transporter==FALSE,,] else tmpDT <- DT


##### data table per parameter
DT_per_parameter <- NULL
DT_tmp1 <- data.table()
DT_tmp2 <- data.table()
DT_tmp3 <- data.table()
DT_tmp1[, c("Step", "Z_threshold", "Max_rxn", "Prior.frac15","Prior.frac10","Prior.frac05","Prior.frac02","Prior.pos.frac.av.rev","Prior.pos.frac.av","Missed","Missed.frac",
            "Max_Tot.Genes","Min_Tot.Genes","Prior.sd15", "Prior.sd10", "Prior.sd05", "Prior.sd02","Prior.pos.frac.sd.rev","Prior.pos.frac.sd","Missed.sd") :=  
            # "Out_top50", "In_top50"
          tmpDT[, list(
            mean(Prioritised15), # Prior.frac15
            mean(Prioritised10), # Prior.frac10
            mean(Prioritised05), # Prior.frac05
            mean(Prioritised02), # Prior.frac02
            mean(Rev_Pos_frac), # Prior.pos.frac.av.rev
            mean(Pos_frac), # Prior.pos.frac.av
            sum(P.value==1)/length(seeds), # Missed; if a gene is missed in 1 seed, it is missed in all --> count once per parametercombination
            mean(P.value==1), # Missed.frac
            max(Total_genes), # Max_Tot.Genes
            min(Total_genes), # Min_Tot.Genes
            sd(Prioritised15), # Prior.sd15
            sd(Prioritised10), # Prior.sd10
            sd(Prioritised05), # Prior.sd05
            sd(Prioritised02), # Prior.sd02
            sd(Rev_Pos_frac), # Prior.pos.frac.sd.rev
            sd(Pos_frac), # Prior.pos.frac.sd
            sd(Position==Total_genes) # Missed.sd
            # sum(!Prioritised50), # Out_top50
            # sum(Prioritised50) # In_top50
          ), by = .(Step, Z_threshold, Max_rxn)]]
# DT_tmp2[, c("Step", "Z_threshold", "Max_rxn","Av_top50","Sd_top50") := tmpDT[Prioritised50 == TRUE, list(mean(Position), sd(Position)), by = .(Step, Z_threshold, Max_rxn)]]
DT_tmp3[, c("Step", "Z_threshold", "Max_rxn","Av_non_missed","Sd_non_missed") := tmpDT[P.value != 1, list(mean(Position), sd(Position)), by = .(Step, Z_threshold, Max_rxn)]]
DT_per_parameter <- merge(DT_tmp1, DT_tmp3, by = c("Step","Z_threshold", "Max_rxn"))
# DT_per_parameter <- merge(DT_per_parameter, DT_tmp2, by = c("Step","Z_threshold", "Max_rxn"))
rm(DT_tmp1,DT_tmp2,DT_tmp3)

DT_per_parameter[, Max_rxn:=factor(Max_rxn, levels = max_rxns)]

# return 1 or 0 depending on whether the parameter combination is the best scoring or not
DT_per_parameter[,c("best15","best10","best05","best02") := list(
  ifelse(Prior.frac15 == max(Prior.frac15), 1, 0),
  ifelse(Prior.frac10 == max(Prior.frac10), 1, 0),
  ifelse(Prior.frac05 == max(Prior.frac05), 1, 0),
  ifelse(Prior.frac02 == max(Prior.frac02), 1, 0))]
# DT_per_parameter[,"Frac.inTop50" := In_top50/(In_top50+Out_top50)]

# Determine best parameters according to disease genes within top 50 (frac. in top 50 >= 0.5)
# DT_tmp <- DT_per_parameter[,Av_top50, Frac.inTop50]
# min_values <- unique(head(sort(DT_tmp[Frac.inTop50 >= 0.5, Av_top50]),5))

# DT_per_parameter[,"best_av50" := Frac.inTop50 >= 0.5 & Av_top50 %in% min_values]
# DT_per_parameter[,"best_av50"] <- DT_per_parameter$Frac.inTop50 >= 0.5 & DT_per_parameter$Av_top50 %in% min_values
# DT_per_parameter[best_av50 == TRUE, best_order := rank(DT_per_parameter[best_av50 == TRUE, Av_top50])]

# Determine best parameters according to non-missed disease genes (frac. missed < 0.5)
DT_tmp <- DT_per_parameter[,Av_non_missed, Missed.frac]
min_values <- unique(head(sort(DT_tmp[Missed.frac < 0.5, Av_non_missed]),5))

DT_per_parameter[,"best_av_NM" := Missed.frac < 0.5 & Av_non_missed %in% min_values]
# DT_per_parameter[,"best_av_NM"] <- DT_per_parameter$Missed.frac < 0.5 & DT_per_parameter$Av_non_missed %in% min_values
DT_per_parameter[best_av_NM == TRUE, best_order_NM := rank(DT_per_parameter[best_av_NM == TRUE, Av_non_missed])]

best_ratio_in_top05 <- head(sort(DT_per_parameter[,Prior.frac05], decreasing = TRUE),5)
DT_per_parameter[,"best_top05" := Prior.frac05 %in% best_ratio_in_top05]
DT_per_parameter[best_top05 == TRUE, best_order_top05 := rank(DT_per_parameter[best_top05 == TRUE, Prior.frac05])]




# For some reason, probably something with the list() creation of the DT, the DT doesn't correctly register as one
# So this 'rectify' is necessary
DT_per_parameter <- data.table(DT_per_parameter)



##### Data table per patient
# patients <- as.vector(unique(DT$PatientID))
patients <- as.vector(unique(tmpDT$PatientID))


DT_per_patient <- data.table()
DT_per_patient[, c("Step", "Z_threshold", "Max_rxn", "PatientID", "Av_rank", "Sd_rank", "Max_rank", "Min_rank","Missed", "Out_top50", "Av_Rev_Rel_rank", "Sd_Rev_Rel_rank",
                   "Av_Rel_rank", "Sd_Rel_rank","Transporter", "P.value_check", "Max_tot_gen", "Min_tot_gen") := 
                 tmpDT[, list(
                   mean(Position), # Av_rank
                   sd(Position), # Sd_rank
                   max(Position), # Max_rank
                   min(Position), # Min_rank
                   sum(P.value == 1), # Missed
                   sum(Position > 50), # Out_top50
                   mean(Rev_Pos_frac), # Av_Rev_Rel_rank
                   sd(Rev_Pos_frac), # Sd_Rev_Rel_rank
                   mean(Pos_frac), # Av_Rel_rank
                   sd(Pos_frac), # Sd_Rel_rank
                   unique(Transporter), # Transporter
                   length(unique(P.value)) == 1, # P.value_check
                   max(Total_genes), # Max_tot_gen
                   min(Total_genes) # Min_tot_gen
                 ), by = .(Step, Z_threshold, Max_rxn, PatientID)]]
DT_per_patient[, Max_rxn:=factor(Max_rxn, levels = max_rxns)]

DT_per_patient[, Av_rank_excl_miss := ifelse(Missed > 0, NA, Av_rank)]
DT_per_patient[, Sd_rank_excl_miss := ifelse(Missed > 0, NA, Sd_rank)]

# DT_per_patient[DT_per_patient$Missed > 0, Sd_rank := -1]
# DT_per_patient[, Colour := ifelse(DT_per_patient[,Transporter], "Red","Black")]

# DT_per_patient[, c("Sd_rank", "Sd_Rev_Rel_rank", "Sd_Rel_rank") := ifelse(Missed, -1, Sd_rank)]

# For some reason, probably something with the list() creation of the DT, the DT doesn't correctly register as one
# So this 'rectify' is necessary
DT_per_patient <- data.table(DT_per_patient)


##### Data table per parameter, missed genes 'removed'
# DT_tmp1 <- copy(DT_per_patient)
DT_tmp2 <- data.table()
DT_tmp2[, c("Step", "Z_threshold", "Max_rxn", "Av_Rank_excl_miss", "Av_Sd_excl_miss") :=
          DT_per_patient[, list(
            mean(Av_rank_excl_miss, na.rm=TRUE), # "Av_Rank"
            mean(Sd_rank_excl_miss, na.rm=TRUE) # "Av_Sd"
          ), by = .(Step, Z_threshold, Max_rxn)]]
DT_per_parameter <- merge(DT_per_parameter, DT_tmp2, by = c("Step","Z_threshold", "Max_rxn"))
DT_per_parameter <- data.table(DT_per_parameter)
rm(DT_tmp2)




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
ggsave(paste0(outdir_name,"/Empty_Facet_",var_name,"_",sub_name,".png"), plot = z,
       width = 250, height = 200, dpi=resolution, units = "mm")

# Simple plot, just top 15, all parameter combinations
p <- ggplot(DT_per_parameter, aes(x = Step)) +
  geom_line(aes(y = Prior.frac15, colour = "line15", group = 1)) +
  geom_point(aes(y = Prior.frac15, colour = "line15"), size=0.5) +
  geom_line(aes(y = Missed.frac, colour = "Missed", group = 1)) +
  geom_point(aes(y = Missed.frac, colour = "Missed"), size=0.5) +
  facet_grid(Max_rxn ~ Z_threshold, labeller = labeller(Max_rxn = Rxn_labs, Z_threshold = Thresh_labs)) +
  theme_dark() + 
  ylab("Disease genes / Total dis. genes") +
  xlab("Max distance to primary reaction") +
  ggtitle("Correct disease gene prioritisation") +
  scale_color_manual(name = "Prioritised Genes",
                     labels = c("In Top 15", "Missed"),
                     values=c(my_greens[4],'RED'))
z <- pretty_plot(p, theme = "dark")
ggsave(paste0(outdir_name,"/Ranks_And_Missed_Top15_",var_name,"_",sub_name,".png"), plot = z,
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
  geom_line(aes(y = Prior.frac5, colour = "line02", group = 1), size = 1.3) +
  geom_point(aes(y = Prior.frac5, colour = "line02"), size=1) +
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
ggsave(paste0(outdir_name,"/Ranks_And_Missed_Single_Par_Comb",var_name,"_",sub_name,".png"), plot = p,
       width = 300, height = 200, dpi=resolution, units = "mm")




#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Combination plot of correctly prioritised and missed genes --------------
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

p <- ggplot(DT_per_parameter, aes(x = Step)) +
  # geom_line(aes(y = Prior.frac15, colour = my_greens[4], group = 1)) +
  
  # Add shaded area for what is expected if the results were random:
  # geom_ribbon(aes(x= as.numeric(Step), ymin=1*15/Max_Tot.Genes, ymax=1*15/Min_Tot.Genes), colour = adjustcolor(my_greens[4],alpha.f=0.5), fill = my_greens[4], alpha="0.5") +
  # geom_ribbon(aes(x= as.numeric(Step), ymin=1*10/Max_Tot.Genes, ymax=1*10/Min_Tot.Genes), colour = adjustcolor(my_greens[3],alpha.f=0.5), fill = my_greens[3], alpha="0.5") +
  # geom_ribbon(aes(x= as.numeric(Step), ymin=1*5/Max_Tot.Genes, ymax=1*5/Min_Tot.Genes), colour = adjustcolor(my_greens[2],alpha.f=0.5), fill = my_greens[2], alpha="0.5") +
  # geom_ribbon(aes(x= as.numeric(Step), ymin=1*2/Max_Tot.Genes, ymax=1*2/Min_Tot.Genes), colour = adjustcolor(my_greens[1],alpha.f=0.6), fill = my_greens[1], alpha="0.6") +
  
  geom_line(aes(y = Prior.frac15, colour = "line15", group = 1)) +
  geom_point(aes(y = Prior.frac15, colour = "line15"), size=0.5) +
  geom_line(aes(y = Prior.frac10, colour = "line10", group = 1)) +
  geom_point(aes(y = Prior.frac10, colour = "line10"), size=0.5) +
  geom_line(aes(y = Prior.frac05, colour = "line05", group = 1)) +
  geom_point(aes(y = Prior.frac05, colour = "line05"), size=0.5) +
  geom_line(aes(y = Prior.frac02, colour = "line02", group = 1)) +
  geom_point(aes(y = Prior.frac02, colour = "line02"), size=0.5) +
  geom_line(aes(y = Missed.frac, colour = "missed", group = 1)) +
  geom_point(aes(y = Missed.frac, colour = "missed"), size=0.5) 

p <- p + facet_grid(Max_rxn ~ Z_threshold, labeller = labeller(Max_rxn = Rxn_labs, Z_threshold = Thresh_labs)) +
  geom_point(data = DT_per_parameter[DT_per_parameter$best15 ==1, ], aes(x=Step, y=Prior.frac15, shape = "Best"), size=8) + 
  scale_shape_manual(name = "",
                     labels = "Best param.\nperformance",
                     values = 42) +
  # geom_point(data = DT_per_parameter[DT_per_parameter$best15 ==1, ],aes(x=Step, y=Prior.frac15),shape = "*", size=8, show.legend = FALSE, colour = "black", fill = "black") +
  geom_point(data = DT_per_parameter[DT_per_parameter$best10 ==1, ],aes(x=Step, y=Prior.frac10),shape = "*", size=8, show.legend = FALSE, colour = "black", fill = "black") +
  geom_point(data = DT_per_parameter[DT_per_parameter$best05 ==1, ],aes(x=Step, y=Prior.frac05),shape = "*", size=8, show.legend = FALSE, colour = "black", fill = "black") +
  geom_point(data = DT_per_parameter[DT_per_parameter$best02 ==1, ],aes(x=Step, y=Prior.frac02),shape = "*", size=8, show.legend = FALSE, colour = "black", fill = "black")
p <- p + theme_dark() + 
  ylab("Disease genes / Total dis. genes") +
  xlab("Max distance to primary reaction") +
  ggtitle("Correct disease gene prioritisation") +
  scale_color_manual(name = "Prioritised Genes",
                     labels = c("In Top 2", "In Top 5", "In Top 10", "In Top 15", "Missed"),
                     values=c(my_greens,'RED'))
p <- pretty_plot(p, theme = "dark")
ggsave(paste0(outdir_name,"/Ranks_And_Missed_",var_name,"_",sub_name,".png"), plot = p,
       width = 300, height = 200, dpi=resolution, units = "mm")


# Average of top 50 genes -------------------------------------------------
p <- ggplot(DT_per_parameter, aes(label=Av_top50)) +
  geom_line(aes(x = Step, y = Av_top50, colour = "Average", group = 1), size = 1.3) +
  geom_point(aes(x = Step, y = Av_top50, colour = "Average"), size=0.5) +
  geom_line(aes(x = Step, y = Out_top50/(In_top50+Out_top50)*15, colour = "Frac.OutTop", group = 1), size = 1.3) +
  geom_point(aes(x = Step, y = Out_top50/(In_top50+Out_top50)*15, colour = "Frac.OutTop"), size=0.5) +
  scale_y_continuous(limits = c(0,50), sec.axis = sec_axis(~./15, name = "Frac. genes outside top 50", breaks = c(0, 0.5, 1))) +
  theme_dark() +
  scale_color_manual(name = "",
                     labels = c("Av. rank of top 50", "Frac. outside top 50"),
                     values=c('Black', my_greens[3]))
p <- p + facet_grid(Max_rxn ~ Z_threshold, labeller = labeller(Max_rxn = Rxn_labs, Z_threshold = Thresh_labs)) 
# geom_point(data = DT_per_parameter[best_av50==TRUE, ],aes(x=Step, y=Av_top50),shape = "*", size=8, show.legend = FALSE, colour = "black")+
p <- p + geom_point(data = DT_per_parameter[as.vector(DT_per_parameter[,"best_order"]==1),],aes(x=Step, y=Av_top50, shape = "Best"), size=8, ) + 
  scale_shape_manual(name = "",
                     labels = "Best param.\nperformance",
                     values = 42)
for(i in c(1:5)){
  p <- p + geom_point(data = DT_per_parameter[as.vector(DT_per_parameter[,"best_order"]==i),],aes(x=Step, y=Av_top50, shape = "Best"), shape = "*", size=8, colour = my_sig_palette[i])
  p <- p + geom_text(data = DT_per_parameter[as.vector(DT_per_parameter[,"best_order"]==i),],
                     aes(x=Step, y=Av_top50, label = signif(Av_top50, digits = digit_significance), group = best_order), 
                     size=3, 
                     # show.legend = FALSE, 
                     colour = my_sig_palette[i],
                     position = position_dodge(width = 2),
                     vjust = -0.5)
}

p <- p + 
  ylab("Average disease gene rank") +
  xlab("Max. distance to primary reaction") +
  ggtitle("Av. gene rank for top 50 genes") + 
  geom_hline(yintercept=10, linetype="dashed", color = "salmon") +
  guides(shape = guide_legend(order = 1),
         colour = guide_legend(order = 2)) 
z <- pretty_plot(p, theme = "dark")

ggsave(paste0(outdir_name,"/Top50_",var_name,"_",sub_name,".png"), plot = z,
       width = 300, height = 200, dpi=resolution, units = "mm")


##### Average rank non-missed genes + total missed genes ------------------
p <- ggplot(DT_per_parameter, aes(label=Av_non_missed)) +
  geom_line(aes(x = Step, y = Av_non_missed, colour = "Average", group = 1), size = 1.3) +
  geom_point(aes(x = Step, y = Av_non_missed, colour = "Average"), size=0.5) +
  geom_line(aes(x = Step, y = Missed.frac*80, colour = "Frac.Missed", group = 1), size = 1.3) +
  geom_point(aes(x = Step, y = Missed.frac*80, colour = "Frac.Missed"), size=0.5) +
  scale_y_continuous(sec.axis = sec_axis(~./80, name = "Frac. genes Missed", breaks = c(0, 0.5, 1))) +
  theme_dark() +
  scale_color_manual(name = "",
                     labels = c("Av. rank of non-missed genes", "Frac. Missed"),
                     values=c('Black', my_greens[3]))
p <- p + facet_grid(Max_rxn ~ Z_threshold, labeller = labeller(Max_rxn = Rxn_labs, Z_threshold = Thresh_labs)) 
# geom_point(data = DT_per_parameter[best_av50==TRUE, ],aes(x=Step, y=Av_top50),shape = "*", size=8, show.legend = FALSE, colour = "black")+
p <- p + geom_point(data = DT_per_parameter[as.vector(DT_per_parameter[,"best_order_NM"]==1),],aes(x=Step, y=Av_non_missed, shape = "Best"), size=8, ) + 
  scale_shape_manual(name = "",
                     labels = "Best param.\nperformance \n(missed frac. <0.5)",
                     values = 42)
# Annotate the 5 best scoring parameter combinations (determined ~l.280)
for(i in c(1:5)){
  p <- p + geom_point(data = DT_per_parameter[as.vector(DT_per_parameter[,"best_order_NM"]==i),],aes(x=Step, y=Av_non_missed, shape = "Best"), shape = "*", size=8, colour = my_sig_palette[i])
  p <- p + geom_text(data = DT_per_parameter[as.vector(DT_per_parameter[,"best_order_NM"]==i),],
                     aes(x=Step, y=Av_non_missed, label = signif(Av_non_missed, digits = digit_significance), group = best_order_NM), 
                     size=3, 
                     # show.legend = FALSE, 
                     colour = my_sig_palette[i],
                     position = position_dodge(width = 2),
                     vjust = -0.5)
}
p <- p + 
  ylab("Average non-missed disease gene rank") +
  xlab("Max. distance to primary reaction") +
  ggtitle("Av. gene rank for non-missed genes") + 
  geom_hline(yintercept=15, linetype="dashed", color = "salmon") +
  guides(shape = guide_legend(order = 1),
         colour = guide_legend(order = 2)) 
p <- pretty_plot(p, theme = "dark")
ggsave(paste0(outdir_name,"/Ranks_Non_Missed_And_Missed",var_name,"_",sub_name,".png"), plot = p,
       width = 300, height = 200, dpi=resolution, units = "mm")


##### all ranks, average rank and average standard deviation --------------
DT_no_missed <- tmpDT[P.value < 1,]
DT_missed <- tmpDT[P.value == 1,]
p <- ggplot(DT_no_missed, aes(x = Step)) +
  # geom_point(aes(y = Position), position=position_dodge(width = 5.5)) +
  geom_jitter(aes(y = Position, colour = "Per-patient rank"), size = 0.005) +
  geom_jitter(data =DT_missed, aes(x = Step, y = Position, fill = "Missed genes"), colour = "blue", size = 0.005, alpha = 0.2) +
  geom_line(data = DT_per_parameter, aes(x = Step, y = Av_non_missed, group = 1, colour = "Average rank")) +
  geom_errorbar(data = DT_per_parameter, aes(x = Step, 
                                             ymax = Av_non_missed + Av_Sd_excl_miss, 
                                             ymin = Av_non_missed - Av_Sd_excl_miss,
                                             colour = "Average rank"), 
                # colour = "cornflowerblue",
                width = 0.5) +
  ylim(c(0,40)) +
  geom_hline(aes(yintercept = 5), colour = "darksalmon", linetype="dashed") +
  facet_grid(Max_rxn ~ Z_threshold, labeller = labeller(Max_rxn = Rxn_labs, Z_threshold = Thresh_labs)) +
  theme_dark() + 
  ylab("Disease gene rank") +
  xlab("Max distance to primary reaction") +
  ggtitle("Stability of method") +
geom_point(data = DT_per_parameter[DT_per_parameter$best_order_top05 ==1, ], aes(x=Step, y=35, shape = "Best"), size=8) + 
  scale_shape_manual(name = "",
                     labels = "Ratio dis.genes\nin top 5",
                     values = 42) 
for(i in c(1:5)){
  p <- p + geom_point(data = DT_per_parameter[as.vector(DT_per_parameter[,"best_order_top05"]==i),],aes(x=Step, y=35, shape = "Best"), shape = "*", size=8)
  p <- p + geom_text(data = DT_per_parameter[as.vector(DT_per_parameter[,"best_order_top05"]==i),],
                     aes(x=Step, y=35, label = signif(Prior.frac05, digits = digit_significance), group = best_order_top05),
                     size=3,
                     # show.legend = FALSE,
                     colour = my_sig_palette[6-i],
                     position = position_dodge(width = 2),
                     vjust = -0.5)
}
p <- p + scale_color_manual(name = "Non-missed genes",
                            labels = c("Average rank", "Per-patient rank"),
                            values=c("red","black")) + 
  scale_fill_manual(name = "Missed genes",
                            labels = c("Per-patient rank"),
                            values=c("blue"))
p <- pretty_plot(p, secondary_y_axis = FALSE, theme = "dark")
ggsave(paste0(outdir_name,"/Average_Patient_And_Ranks_And_Missed_",var_name,"_",sub_name,".png"), plot = p,
       width = 300, height = 200, dpi=resolution, units = "mm")




#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Per patient plots -------------------------------------------------------
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

##### Average rank per patient --------------------------------------------
p <- ggplot(data = DT_per_patient, aes(x = Step, y = Av_rank)) +
  geom_bar(position = "dodge", stat="identity", aes(fill = PatientID)) +
  geom_line(data = DT_per_parameter, aes(x = Step, y = Missed/0.25, group = 1, colour = "Genes Missed")) +
  scale_y_continuous(sec.axis = sec_axis(~.*0.25, name = "Tot.Genes Missed")) +
  ylab("Average Rank") +
  xlab("Max. distance to primary reaction") +
  ggtitle("Disease gene rank per patient & parameter combination") +
  labs(colour = "", fill = "Patients") +
  guides(colour = guide_legend(order = 1),
         fill = guide_legend(order = 2)) +
  facet_grid(Max_rxn ~ Z_threshold, labeller = labeller(Max_rxn = Rxn_labs, Z_threshold = Thresh_labs))
p <- pretty_plot(p)
ggsave(paste0(outdir_name,"/Av_Rank_Per_Patient_",var_name,"_",sub_name,".png"), plot = p,
       width = 300, height = 200, dpi=resolution, units = "mm")


##### Standard deviation of average rank per patient ----------------------
p <- ggplot(data = DT_per_patient, aes(x = Step, y = Sd_rank)) +
  geom_bar(position = "dodge", stat="identity", aes(fill = PatientID)) +
  # geom_line(data = DT_per_parameter, aes(x = Step, y = Missed/4, group = 1, colour = "Genes Missed")) +
  # scale_y_continuous(sec.axis = sec_axis(~.*4, name = "Tot.Genes Missed")) +
  # scale_y_continuous(sec.axis = sec_axis(~.*4, name = "Tot.Genes Missed", breaks = seq(0, 50, 10))) +
  ylab("Sd disease gene rank") +
  xlab("Max. distance to primary reaction") +
  ggtitle("Disease gene St.dev. per patient & parameter combination") +
  labs(colour = "", fill = "Patients") +
  guides(colour = guide_legend(order = 1),
         fill = guide_legend(order = 2)) +
  facet_grid(Max_rxn ~ Z_threshold, labeller = labeller(Max_rxn = Rxn_labs, Z_threshold = Thresh_labs))
p <- pretty_plot(p)
ggsave(paste0(outdir_name,"/Sd_Rank_Per_Patient_",var_name,"_",sub_name,".png"), plot = p,
       width = 300, height = 200, dpi=resolution, units = "mm")


##### Average relative rank per patient -----------------------------------
p <- ggplot(data = DT_per_patient, aes(x = Step, y = Av_Rel_rank)) +
  geom_bar(position = "dodge", stat="identity", aes(fill = PatientID)) +
  # geom_line(data = DT_per_parameter, aes(x = Step, y = Missed/50, group = 1, colour = "Genes Missed")) +
  # scale_y_continuous(sec.axis = sec_axis(~.*50, name = "Tot.Genes Missed")) +
  ylab("Average Relative Rank (rank/Tot. genes") +
  xlab("Max. distance to primary reaction") +
  ggtitle("Disease gene rel.rank per patient & parameter combination") +
  labs(colour = "", fill = "Patients") +
  guides(colour = guide_legend(order = 1),
         fill = guide_legend(order = 2)) +
  facet_grid(Max_rxn ~ Z_threshold, labeller = labeller(Max_rxn = Rxn_labs, Z_threshold = Thresh_labs))
p <- pretty_plot(p)
ggsave(paste0(outdir_name,"/Av_Rel_Rank_Per_Patient_",var_name,"_",sub_name,".png"), plot = p,
       width = 300, height = 200, dpi=resolution, units = "mm")


##### Standard deviation of average relative rank per patient -------------
p <- ggplot(data = DT_per_patient, aes(x = Step, y = Sd_Rel_rank)) +
  geom_bar(position = "dodge", stat="identity", aes(fill = PatientID)) +
  # geom_line(data = DT_per_parameter, aes(x = Step, y = Missed/800, group = 1, colour = "Genes Missed")) +
  # scale_y_continuous(sec.axis = sec_axis(~.*800, name = "Tot.Genes Missed")) +
  ylab("Sd disease gene relative Rank") +
  xlab("Max. distance to primary reaction") +
  ggtitle("Disease gene St.dev. per patient & parameter combination") +
  labs(colour = "", fill = "Patients") +
  guides(colour = guide_legend(order = 1),
         fill = guide_legend(order = 2)) +
  facet_grid(Max_rxn ~ Z_threshold, labeller = labeller(Max_rxn = Rxn_labs, Z_threshold = Thresh_labs))
p <- pretty_plot(p)
ggsave(paste0(outdir_name,"/Sd_Rel_Rank_Per_Patient_",var_name,"_",sub_name,".png"), plot = p,
       width = 300, height = 200, dpi=resolution, units = "mm")


##### Average relative rank per patient, reversed (1 = rank 1) ------------
p <- ggplot(data = DT_per_patient, aes(x = Step, y = Av_Rev_Rel_rank)) +
  geom_bar(position = "dodge", stat="identity", aes(fill = PatientID)) +
  # geom_line(data = DT_per_parameter, aes(x = Step, y = Missed/50, group = 1, colour = "Genes Missed")) +
  # scale_y_continuous(sec.axis = sec_axis(~.*50, name = "Tot.Genes Missed")) +
  ylab("Average reverse Relative Rank (1-((rank-1)/(Tot. genes in set-1))") +
  xlab("Max. distance to primary reaction") +
  ggtitle("Disease gene rel.rank per patient & parameter combination") +
  labs(colour = "", fill = "Patients") +
  guides(colour = guide_legend(order = 1),
         fill = guide_legend(order = 2)) +
  facet_grid(Max_rxn ~ Z_threshold, labeller = labeller(Max_rxn = Rxn_labs, Z_threshold = Thresh_labs))
p <- pretty_plot(p)
ggsave(paste0(outdir_name,"/Av_Rev_Rel_Rank_Per_Patient_",var_name,"_",sub_name,".png"), plot = p,
       width = 300, height = 200, dpi=resolution, units = "mm")


##### St. dev. of average relative (reversed) rank per patient ------------
p <- ggplot(data = DT_per_patient, aes(x = Step, y = Sd_Rev_Rel_rank)) +
  geom_bar(position = "dodge", stat="identity", aes(fill = PatientID)) +
  # geom_line(data = DT_per_parameter, aes(x = Step, y = Missed/800, group = 1, colour = "Genes Missed")) +
  # scale_y_continuous(sec.axis = sec_axis(~.*800, name = "Tot.Genes Missed")) +
  ylab("Sd disease gene reverse relative Rank") +
  xlab("Max. distance to primary reaction") +
  ggtitle("Disease gene reverse St.dev. per patient & parameter combination") +
  labs(colour = "", fill = "Patients") +
  guides(colour = guide_legend(order = 1),
         fill = guide_legend(order = 2)) +
  facet_grid(Max_rxn ~ Z_threshold, labeller = labeller(Max_rxn = Rxn_labs, Z_threshold = Thresh_labs))
p <- pretty_plot(p)
ggsave(paste0(outdir_name,"/Sd_Rev_Rel_Rank_Per_Patient_",var_name,"_",sub_name,".png"), plot = p,
       width = 300, height = 200, dpi=resolution, units = "mm")


##### Average standard deviation (normal rank) ----------------------------
p <- ggplot(data = DT_per_parameter, aes(x = Step)) +
  # geom_line(aes(y = Av.Sd, group = 1, colour = "St.dev.")) +
  geom_errorbar(aes(ymax = Av.Rank + Av.Sd, ymin = Av.Rank-Av.Sd)) +
  geom_line(aes(y = Av.Rank, group = 1, colour = "Rank")) +
  geom_line(aes(y = Missed/4, group = 1, colour = "Missed")) +
  # geom_ribbon(aes(x= as.numeric(Step), ymin = 1/2 * Min_tot.gen, ymax=1/2 * Max_tot.gen, fill = "Random_genes")) +
  scale_y_continuous(sec.axis = sec_axis(~.*4, name = "Tot.Genes Missed")) +
  # scale_y_continuous(limits = c(0,8)) +
  ylab("Av. rank disease gene (excl. missed genes)") +
  xlab("Max. distance to primary reaction") +
  ggtitle("Average St.dev. per patient & parameter combination") +
  # labs(colour = "", fill = "Patients") +
  # guides(colour = guide_legend(order = 1),
         # fill = guide_legend(order = 2)) +
  facet_grid(Max_rxn ~ Z_threshold, labeller = labeller(Max_rxn = Rxn_labs, Z_threshold = Thresh_labs))
p <- pretty_plot(p, secondary_y_axis = TRUE)
ggsave(paste0(outdir_name,"/Av_Rank_Sd_Per_Patient_",var_name,"_",sub_name,".png"), plot = p,
       width = 300, height = 200, dpi=resolution, units = "mm")

