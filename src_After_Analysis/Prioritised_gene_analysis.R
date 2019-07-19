library("rstudioapi")
library("Matrix.utils")
# library("xlsx")
library("XLConnect")
code_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
# Z_neg_thresholds <- c(-0.5, -1, -1.5, -3, -10)
# Z_pos_thresholds <- c(1, 1.5, 2, 3, 10)
Z_thresholds <- c("-0.5, 1", "-1, 1.5", "-1.5, 2", "-3, 3")
max_rxns <- c(6, 8, 10, 12,15)
steps <- c(0,1,2,3,4)

load(paste0(code_dir,"/../Data/Crossomics_DBS_Marten_Training.RData"))


# columns to paste together
xls_data$Patient <- unlist(lapply(xls_data$Patient.number, function(x) unlist(strsplit(x, split = "\\."))[1]))
cols <- c( 'Dataset' , 'Patient' )

# create a new column `x` with the three columns collapsed together
xls_unique <- xls_data[!duplicated(apply( xls_data[ , cols ] , 1 , paste , collapse = "_" )),]

fixed_patients <- unlist(lapply(strsplit(xls_unique$Patient[nchar(xls_unique$Patient) == 2], split = ""), function (x) paste0(x[1],"0",x[2])))
xls_unique$Patient[nchar(xls_unique$Patient) == 2] <- fixed_patients

DT <- NULL

# patients_not_done <- c("P9.1")

for (i in 1:nrow(xls_unique)){
  patient <- xls_unique$Patient[i]
  # if (xls_unique$Patient.number[i] %in% patients_not_done) next
  cat("Patient:",patient,"\n")
  dis_gene <- xls_unique$Gene[i]
  dis_gene <- unlist(strsplit(dis_gene, split = "; "))
  patient_folder <- paste0("2019-07-18/", patient, "_",xls_unique$Dataset[i])
  path <- paste0(code_dir,"/../Results/",patient_folder)
  
  for(maxrxn in max_rxns){
    
    for(threshs in Z_thresholds){
      
      thresh <- unlist(strsplit(threshs, ", "))
      
      for(step in steps){
        
        path2 <- paste0(path,"/maxrxn",maxrxn,"_thresh_n",thresh[1],"_p",thresh[2],"_step_",step,"/MSEA_results.xls")
        # if(!file.exists(path2)) path2 <- paste0(path,"/maxrxn",maxrxn,"thresh_n",thresh[1],"_p",thresh[2],"_step_",step,"/MSEA_results.xls")
        # MSEA_results <- read.xlsx(path2, sheetIndex = 1)
        wb <- loadWorkbook(path2)
        MSEA_results <- readWorksheet(wb, sheet = 1, startRow = 0, endRow = 0, startCol = 0, endCol = 0)
        
        # Take the best scoring (possible) disease gene (when more than 1 is present)
        if(any(dis_gene %in% MSEA_results$metabolite.set)){
          dis_pos <- min(unlist(lapply(dis_gene, function(x) grep(x, MSEA_results$metabolite.set))))
        } else {
          dis_pos <- nrow(MSEA_results)
        }
        # dis_pos <- grep(dis_gene, MSEA_results$metabolite.set)

        gene_set <- nrow(MSEA_results)
        
        patient_DT <- data.table(Patient = patient,
                                 Gene = paste(dis_gene,collapse = ";"),
                                 Position = dis_pos,
                                 Total_genes = gene_set,
                                 # Pos_frac = 1-((dis_pos-1)/(gene_set-1)),
                                 Pos_frac = dis_pos/gene_set,
                                 Z_threshold = threshs,
                                 Step = step,
                                 Max_rxn = maxrxn)
        
        ###########################################################################
        # Make datatable of all patients ------------------------------------------
        ###########################################################################
        
        DT <- rbind(DT, patient_DT)

      }
    }
  }
}





# To be in separate script
###########################################################################
# Make plot of effectiveness of variables ---------------------------------
###########################################################################

DT$Z_threshold_order <- factor(DT$Z_threshold, levels = Z_thresholds)

library(ggplot2)
p <- ggplot(DT, aes(x=factor(Z_threshold_order), y=Pos_frac, color=as.factor(Max_rxn))) +
  geom_boxplot() + 
  geom_point(position = "jitter") +
  # geom_point(position = "jitter", aes(shape=as.factor(Step))) +
  scale_x_discrete(name ="Z-value threshold")
# p + facet_grid(. ~ Step)
p <- ggplot(DT, aes(x = factor(Step), y = Pos_frac)) +
  geom_dotplot(binaxis = "y", stackdir = "centerwhole", binwidth = 0.005) +
  scale_x_discrete(name ="Z-value threshold") +
  geom_violin(alpha = 0.1)

png(file = "/Users/mkerkho7/DIMS2_repo/Crossomics/TestResults/testRplots.png", width = 2000, height = 2000)
p + facet_grid(Max_rxn ~ Z_threshold_order)
dev.off()



# Subset of rows ----------------------------------------------------------

DT_small <- DT[Z_threshold == "-1, 1.5" & (Step == 3) & (Max_rxn == 10 | Max_rxn == 12)]
DT_small$Z_threshold_order <- factor(DT_small$Z_threshold, levels = Z_thresholds)

q <- ggplot(DT_small, aes(x = factor(Step), y = Pos_frac)) +
  geom_dotplot(binaxis = "y", stackdir = "center", binwidth = 0.005) +
  scale_x_discrete(name ="Z-value threshold") +
  geom_violin(alpha = 0.1)


library(ggrepel)

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


