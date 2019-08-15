###########################################################################
# Libraries ---------------------------------------------------------------
###########################################################################

library("rstudioapi")
library("Matrix.utils")
library("data.table")

# library("XLConnect") # Only for older results. Switched to .RData files
# library("xlsx")


###########################################################################
# Set variables -----------------------------------------------------------
###########################################################################

code_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
Z_thresholds <- c("-0.5, 1", "-1, 1.5", "-1.5, 2", "-3, 3")
max_rxns <- c(8, 10, 12, 15, 17, 19)
steps <- c(0,1,2,3,4,5)
date <- "2019-08-02"
seeds <- c(2341, 6734892, 83, 698, 991)
patients_not_done <- NULL # in format c("P38.1", "P39.1","P40.1")


###########################################################################
# Prepare data ------------------------------------------------------------
###########################################################################

load(paste0(code_dir,"/../Data/Crossomics_DBS_Marten_Training_inclProt_function.RData"))
# xls_data <- xlsx::read.xlsx(paste0(code_dir,"/../Data/Crossomics_DBS_Marten_Training_inclProt_function.xlsx"), sheetIndex = 1, colIndex = c(1:7), rowIndex = c(1:107), stringsAsFactors = FALSE)
# save(xls_data, file = paste0(code_dir,"/../Data/Crossomics_DBS_Marten_Training_inclProt_function.RData"))



# columns to paste together
xls_data$Patient <- unlist(lapply(xls_data$Patient.number, function(x) unlist(strsplit(x, split = "\\."))[1]))
cols <- c( 'Dataset' , 'Patient' )

# create a new column `x` with the three columns collapsed together
xls_unique <- xls_data[!duplicated(apply( xls_data[ , cols ] , 1 , paste , collapse = "_" )),]

# Fix numbering of some patients
fixed_patients <- unlist(lapply(strsplit(xls_unique$Patient[nchar(xls_unique$Patient) == 2], split = ""), function(x) paste0(x[1],"0",x[2])))
xls_unique$Patient[nchar(xls_unique$Patient) == 2] <- fixed_patients

DT <- NULL



###########################################################################
# Collate results to 1 data table -----------------------------------------
###########################################################################

for (i in 1:nrow(xls_unique)){
  patient <- xls_unique$Patient[i]
  prot_func <- xls_unique$Gene.product.function[i]
  if (length(patients_not_done) > 0){
    if(xls_unique$Patient.number[i] %in% patients_not_done) next
  }
  cat("Patient:",patient,"\n")
  dis_gene <- xls_unique$Gene[i]
  dis_gene <- unlist(strsplit(dis_gene, split = "; "))
  is_trans <- xls_unique$Gene.product.function[i] == "transporter"
  patient_folder <- paste0(date,"/", patient, "_",xls_unique$Dataset[i])
  path <- paste0(code_dir,"/../Results/",patient_folder)
  
  for(seed in seeds){

    # Check if patient was included in the seed, continue to next one otherwise
    if(!dir.exists(paste0(path,"/seed",seed))) next
    
    for(maxrxn in max_rxns){
      
      for(threshs in Z_thresholds){
        
        thresh <- unlist(strsplit(threshs, ", "))
        
        for(step in steps){
          
          # path2 <- paste0(path,"/maxrxn",maxrxn,"_thresh_n",thresh[1],"_p",thresh[2],"_step_",step,"/MSEA_results.xls")
          # # if(!file.exists(path2)) path2 <- paste0(path,"/maxrxn",maxrxn,"thresh_n",thresh[1],"_p",thresh[2],"_step_",step,"/MSEA_results.xls")
          # # MSEA_results <- read.xlsx(path2, sheetIndex = 1)
          # wb <- loadWorkbook(path2)
          # MSEA_results <- readWorksheet(wb, sheet = 1, startRow = 0, endRow = 0, startCol = 0, endCol = 0)
          
          path2 <- paste0(path,"/seed",seed,"/maxrxn",maxrxn,"_thresh_n",thresh[1],"_p",thresh[2],"_step_",step,"/MSEA_results.RData")
          load(path2)
          MSEA_results <- metSetResult
          
          # Take the best scoring (possible) disease gene (when more than 1 is present) or a 'last place' when it is absent
          if(any(dis_gene %in% MSEA_results$metabolite.set)){
            dis_pos <- min(unlist(lapply(dis_gene, function(x) grep(x, MSEA_results$metabolite.set))))
            p.value <- MSEA_results$p.value[dis_pos]
            if(p.value==1){
              rank <- nrow(MSEA_results)
            } else {
              rank <- frank(MSEA_results$p.value)[dis_pos]
            }
          } else {
            # dis_pos <- nrow(MSEA_results)
            rank <- nrow(MSEA_results)
            p.value <- 1
          }
          # dis_pos <- grep(dis_gene, MSEA_results$metabolite.set)
          
          gene_set <- nrow(MSEA_results)
          
          patient_DT <- data.table(Patient = patient,
                                   Gene = paste(dis_gene,collapse = ";"),
                                   Protein_function = prot_func,
                                   Transporter = is_trans,
                                   Position = rank,
                                   Total_genes = gene_set,
                                   # Rev_Pos_frac = 1-((dis_pos-1)/(gene_set-1)),
                                   Rev_Pos_frac = 1-((rank-1)/(gene_set-1)),
                                   # Pos_frac = dis_pos/gene_set,
                                   Pos_frac = rank/gene_set,
                                   Z_threshold = threshs,
                                   Step = step,
                                   Max_rxn = maxrxn,
                                   P.value = p.value,
                                   Seed = seed)
          
          ###########################################################################
          # Make datatable of all patients ------------------------------------------
          ###########################################################################
          
          DT <- rbind(DT, patient_DT)
          
        }
      }
    }
  }
}





# To be in separate script
###########################################################################
# Make plot of effectiveness of variables ---------------------------------
###########################################################################

# DT$Z_threshold_order <- factor(DT$Z_threshold, levels = Z_thresholds)
DT[, Z_threshold:=as.factor(Z_threshold)]
DT[, Protein_function:=as.factor(Protein_function)]

library(ggplot2)
# p <- ggplot(DT, aes(x=Z_threshold, y=Pos_frac, color=as.factor(Max_rxn))) +
p <- ggplot(DT, aes(x=Z_threshold, y=Pos_frac, color=Protein_function)) +
  geom_boxplot() + 
  geom_point(position = "jitter") +
  # geom_point(position = "jitter", aes(shape=as.factor(Step))) +
  scale_x_discrete(name ="Z-value threshold")
# p + facet_grid(. ~ Step)
p <- ggplot(DT, aes(x = factor(Step), y = Rev_Pos_frac, fill = Transporter)) +
  geom_boxplot(width=0.4, fatten = 3, outlier.size = 0.1) +
  geom_violin(fill = "white", alpha = 0.1) +
  geom_dotplot(binaxis = "y", stackdir = "centerwhole", binwidth = 0.015) +
  scale_x_discrete(name ="Z-value threshold") +
  theme_light()
# p <- ggplot(DT, aes(x = factor(Step), y = Position, fill = Transporter)) +
#   geom_boxplot(width=0.1) +
#   geom_violin(fill = "white", alpha = 0.1) +
#   geom_dotplot(binaxis = "y", stackdir = "centerwhole", binwidth = 1.5) +
#   scale_x_discrete(name ="Z-value threshold")
  
  
png(file = "/Users/mkerkho7/DIMS2_repo/Crossomics/TestResults/AllPars.png", width = 1200, height = 800, res=1200)
p + facet_grid(Max_rxn ~ Z_threshold, labeller = label_both)
dev.off()
ggsave("/Users/mkerkho7/DIMS2_repo/Crossomics/TestResults/AllPars.png", 
       width = 300, height = 200, dpi=600, units = "mm")




# Subset of rows ----------------------------------------------------------

# DT_small <- DT[Z_threshold == "-1, 1.5" & (Step == 3) & (Max_rxn == 10 | Max_rxn == 12)]
# No transporter disease genes
DT_small <- DT[Transporter != TRUE]
p <- ggplot(DT_small, aes(x = factor(Step), y = Rev_Pos_frac, fill = Gene)) +
  stat_summary(fun.data="mean_sdl",  fun.args = list(mult=1), 
               geom="pointrange", color = "red") +
  geom_boxplot(fill = "red", alpha = 0.1, width=0.2, fatten = 5) +
  geom_violin(fill = "white", alpha = 0.1) +
  geom_dotplot(binaxis = "y", stackdir = "centerwhole", binwidth = 0.008) +
  scale_x_discrete(name ="Z-value threshold") +
  theme_linedraw()
png(file = "/Users/mkerkho7/DIMS2_repo/Crossomics/TestResults/AllParsNoTransporters.png", width = 2000, height = 2000)
p + facet_grid(Max_rxn ~ Z_threshold)
dev.off()

png(file = "/Users/mkerkho7/DIMS2_repo/Crossomics/TestResults/ByProtFunc.png", width = 800, height = 400)
ggplot(DT, aes(x = factor(Step), y = Rev_Pos_frac, fill = Protein_function)) +
  # geom_dotplot(binaxis = "y", stackdir = "center", binwidth = 0.008) +
  geom_boxplot() +
  theme_dark()
dev.off()

DT_small$Z_threshold_order <- factor(DT_small$Z_threshold, levels = Z_thresholds)

q <- ggplot(DT_small, aes(x = factor(Step), y = Pos_frac)) +
  geom_dotplot(binaxis = "y", stackdir = "center", binwidth = 0.005) +
  scale_x_discrete(name ="Z-value threshold") +
  geom_violin(alpha = 0.1)


# library(ggrepel)

bw <- 0.005

built <- ggplot_build(q)
point.pos <- built$data[[1]]

# Order rows 
idx <- order(DT_small[,Pos_frac])
DT_small <- DT_small[idx]

# Get the dimensions of the target device in pixels
size <- dev.size(units = 'px')
# Get the range of x and y domain values
extent <- with(built$layout$panel_params[[1]], abs(c(diff(x.range), diff(y.range))))

DT_small$ytext <- point.pos$y
DT_small$xtext <- point.pos$x

dev.off()
png(file = "/Users/mkerkho7/DIMS2_repo/Crossomics/TestResults/testRplotsmallstep3.png", width = 500, height = 1200)
q <- ggplot(DT_small, aes(x = factor(Step), y = Pos_frac)) +
  geom_dotplot(binaxis = "y", binwidth = bw, stackdir = 'center') +
  geom_text_repel(
    aes(xtext, ytext, label = paste(DT_small$Patient,DT_small$Gene)),
    box.padding = unit(2 * size[1] * bw / extent[1], 'points'),
    color = 'red'
  ) +
  scale_x_discrete(name ="Step") +
  geom_violin(alpha = 0.1)
q + facet_grid(Max_rxn ~ .)
dev.off()


###########################################################################
# Correct / incorrect prioritisation datatable ----------------------------
###########################################################################
# gene correctly prioritised?
DT$Prioritised40 <- DT$Position <= 40
DT$Prioritised20 <- DT$Position <= 20
DT$Prioritised10 <- DT$Position <= 10
DT$Prioritised5 <- DT$Position <= 5

DT_noTrans <- DT[!DT$Transporter]

tmpDT <- DT

DT_prioritised <- NULL
for(maxrxn in max_rxns){
  for(threshs in Z_thresholds){
    thresh <- unlist(strsplit(threshs, ", "))
    for(step in steps){
      prior_frac40 <- mean(tmpDT[Step==step & Z_threshold==threshs & Max_rxn==maxrxn, Prioritised40])
      prior_frac20 <- mean(tmpDT[Step==step & Z_threshold==threshs & Max_rxn==maxrxn, Prioritised20])
      prior_frac10 <- mean(tmpDT[Step==step & Z_threshold==threshs & Max_rxn==maxrxn, Prioritised10])
      prior_frac5 <- mean(tmpDT[Step==step & Z_threshold==threshs & Max_rxn==maxrxn, Prioritised5])
      prior_av <- mean(tmpDT[Step==step & Z_threshold==threshs & Max_rxn==maxrxn, Rev_Pos_frac])
      missed_dis <- mean(tmpDT[Step==step & Z_threshold==threshs & Max_rxn==maxrxn, Position==Total_genes])
      sd_prior40 <- sd(tmpDT[Step==step & Z_threshold==threshs & Max_rxn==maxrxn, Prioritised40])
      sd_prior20 <- sd(tmpDT[Step==step & Z_threshold==threshs & Max_rxn==maxrxn, Prioritised20])
      sd_prior10 <- sd(tmpDT[Step==step & Z_threshold==threshs & Max_rxn==maxrxn, Prioritised10])
      sd_prior5 <- sd(tmpDT[Step==step & Z_threshold==threshs & Max_rxn==maxrxn, Prioritised5])
      sd_av <- sd(tmpDT[Step==step & Z_threshold==threshs & Max_rxn==maxrxn, Rev_Pos_frac])
      sd_missed <- sd(tmpDT[Step==step & Z_threshold==threshs & Max_rxn==maxrxn, Position==Total_genes])
      varDT <- data.table(Z_threshold = threshs,
                          Step = step,
                          Max_rxn = maxrxn,
                          Prior_frac40 = prior_frac40,
                          Prior_frac20 = prior_frac20,
                          Prior_frac10 = prior_frac10,
                          Prior_frac5 = prior_frac5,
                          Sd_prior40 = sd_prior40,
                          Sd_prior20 = sd_prior20,
                          Sd_prior10 = sd_prior10,
                          Sd_prior5 = sd_prior5,
                          Prior_av = prior_av,
                          
                          Missed_frac = missed_dis,
                          Sd_missed = sd_missed
                          )
      DT_prioritised <- rbind(DT_prioritised, varDT)
    }
  }
}


###########################################################################
# Combination plot of correctly prioritised and missed genes --------------
###########################################################################
library(RColorBrewer)
my_palette = rev(brewer.pal(5, "Greens"))[c(2:5)]

DT_prioritised[,"signif40"] <- ifelse(DT_prioritised[,"Prior_frac40"] == max(DT_prioritised[,"Prior_frac40"]),1,0)
DT_prioritised[,"signif20"] <- ifelse(DT_prioritised[,"Prior_frac20"] == max(DT_prioritised[,"Prior_frac20"]),1,0)
DT_prioritised[,"signif10"] <- ifelse(DT_prioritised[,"Prior_frac10"] == max(DT_prioritised[,"Prior_frac10"]),1,0)
DT_prioritised[,"signif5"] <- ifelse(DT_prioritised[,"Prior_frac5"] == max(DT_prioritised[,"Prior_frac5"]),1,0)


png(file = "/Users/mkerkho7/DIMS2_repo/Crossomics/TestResults/PriorMissed.png", width = 1600, height = 800)
p <- ggplot(DT_prioritised) +
  geom_line(aes(x = Step, y = Prior_frac40, colour = my_palette[4])) +
  geom_point(aes(x = Step, y = Prior_frac40, colour = my_palette[4]), size=0.5) +
  # geom_errorbar(aes(x = Step, ymax = Prior_frac40 + Sd_prior40, ymin = Prior_frac40 - Sd_prior40), position = "dodge") +
  geom_line(aes(x = Step, y = Prior_frac20, colour = my_palette[3])) +
  geom_point(aes(x = Step, y = Prior_frac20, colour = my_palette[3]), size=0.5) +
  geom_line(aes(x = Step, y = Prior_frac10, colour = my_palette[2])) +
  geom_point(aes(x = Step, y = Prior_frac10, colour = my_palette[2]), size=0.5) +
  geom_line(aes(x = Step, y = Prior_frac5, colour = my_palette[1])) +
  geom_point(aes(x = Step, y = Prior_frac5, colour = my_palette[1]), size=0.5) +
  geom_line(aes(x = Step, y = Missed_frac, color = "Red")) +
  geom_point(aes(x = Step, y = Missed_frac, colour = "Red"), size=0.5) +
  theme_dark() + 
  # scale_color_discrete(name = "Series", labels = c("Prior 40", "Prior 20", "Prior 10", "Prior 5", "Missed"))+ 
  scale_color_manual(name = "Frac. Series", 
                     labels = c("Prior 5", "Prior 10", "Prior 20", "Prior 40", "Missed/p=1"),
                     values=c(my_palette[1], my_palette[2], my_palette[3], my_palette[4],'RED'))
p + facet_grid(Max_rxn ~ Z_threshold, labeller = label_both) +
  geom_point(data = DT_prioritised[DT_prioritised$signif40 ==1, ],aes(x=Step, y=Prior_frac40),shape = "*", size=8, show.legend = FALSE, colour = "black") +
  geom_point(data = DT_prioritised[DT_prioritised$signif20 ==1, ],aes(x=Step, y=Prior_frac20),shape = "*", size=8, show.legend = FALSE, colour = "black") +
  geom_point(data = DT_prioritised[DT_prioritised$signif10 ==1, ],aes(x=Step, y=Prior_frac10),shape = "*", size=8, show.legend = FALSE, colour = "black") +
  geom_point(data = DT_prioritised[DT_prioritised$signif5 ==1, ],aes(x=Step, y=Prior_frac5),shape = "*", size=8, show.legend = FALSE, colour = "black")
dev.off()


###########################################################################
# Make plot of correctly prioritised genes --------------------------------
###########################################################################

png(file = "/Users/mkerkho7/DIMS2_repo/Crossomics/TestResults/CorrectPrior.png", width = 800, height = 400)
p <- ggplot(DT_prioritised, aes(x = factor(Step), y = Prior_frac)) +
  # geom_dotplot(binaxis = "y", stackdir = "center", binwidth = 0.008) +
  geom_boxplot() +
  # geom_errorbar(aes(ymax = Prior_frac + Sd, ymin = Prior_frac - Sd),
  #               position = "dodge")
  theme_dark() 
p + facet_grid(Max_rxn ~ Z_threshold)
dev.off()


###########################################################################
# Make plot of Disease genes not picked up --------------------------------
###########################################################################

png(file = "/Users/mkerkho7/DIMS2_repo/Crossomics/TestResults/MissedGene.png", width = 800, height = 400)
p <- ggplot(DT_missedgene, aes(x = factor(Step), y = Frac_missed)) +
  # geom_dotplot(binaxis = "y", stackdir = "center", binwidth = 0.008) +
  geom_boxplot() +
  # geom_errorbar(aes(ymax = Frac_missed + Sd, ymin = Frac_missed - Sd),
  #               position = "dodge")
  theme_dark() 
p + facet_grid(Max_rxn ~ Z_threshold)
dev.off()
