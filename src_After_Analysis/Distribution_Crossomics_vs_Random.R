# Random chance to get a gene in the top 10



library(ggplot2)
library(data.table)
library(stringr)

# Date of data
date <- "2019-12-10"

best_pars <- c("4;-3,3;15","4;-3,3;19","5;-3,3;15","4;-3,3;17","5;-3,3;19")

code_dir <- paste0(dirname(rstudioapi::getActiveDocumentContext()$path),"/../Results/")
outdir_name <- paste0(dirname(rstudioapi::getActiveDocumentContext()$path),"/../Plots/",date)


#####
# Distributions of best parameter combinations
# DT_per_parameter_tra_val <- readRDS(paste0(code_dir,date,"/MSEA_DT_per_parameter_tra_val.RDS"))
DT <- readRDS(paste0(code_dir,date,"/MSEA_DT_compiled_all.RDS"))
DT[,c("Prioritised15","Prioritised10","Prioritised05","Prioritised02") := 
     list(Position <= 15, Position <= 10, Position <= 5, Position <= 2)]

DT_distribution <- data.table()
DT_distribution[, c("Step", "Z_threshold", "Max_rxn", "Seed", "Prior.frac10") := 
                  DT[Include == TRUE, list(
                    mean(Prioritised10)
                  ), by = .(Step, Z_threshold, Max_rxn, Seed)]]

DT_distribution[, Par.ID := do.call(paste, c(list(Step, Z_threshold, Max_rxn), sep = ";"))]
DT_distribution[, Par.ID := factor(str_replace(Par.ID, " ", ""))]
DT_distribution[, Include := DT_distribution$Par.ID %in% best_pars]


#####
# Distribution by random chance
# In general, there are 201 genes in the list and 1 of them is considered the disease gene
# distribution <- NULL
gene <- NULL
# success <- NULL
top10 <- list()

for(j in c(1:97)){
  gene[j] <- sample(201, 1)
}

for(j in c(1:1000)){
  top10[[j]] <- sample(201, 10, replace = FALSE)
}

sampleseed <- c(1:1000)

DT_random <- data.table(Random_gene = rep(gene, each = 1000), Random_seed =rep(c(1:1000), 97))
for(i in c(1:1000)){
  DT_random[Random_seed == i, Prioritised10 := Random_gene %in% top10[[i]]]
}
DT_distribution_random <- data.table()
DT_distribution_random[,c("Random_seed", "Prior.frac10", "Par.ID") := DT_random[, list(
  mean(Prioritised10),
  "Random"), by = .(Random_seed)]]


#####
# Combined real and simulated prioritization
DT_combined <- rbind(DT_distribution_random, DT_distribution[Include == TRUE,], fill = TRUE)
DT_combined$Par.ID <- factor(DT_combined$Par.ID, 
                             levels = c("Random",best_pars), ordered = TRUE)

#####
# Visualization
p <- ggplot(data = DT_combined, aes(x = Par.ID, y = Prior.frac10)) +
  geom_boxplot() +
  labs(y = "Average correct prioritization", x = "Parameter combination")

ggsave(paste0(outdir_name,"/Cross-omics_vs_random_performance.svg"), plot = p,
       width = 200, height = 50, units = "mm")
  
                    