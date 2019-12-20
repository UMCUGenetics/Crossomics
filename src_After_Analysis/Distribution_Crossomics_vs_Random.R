# Random chance to get a gene in the top 10

# In general, there are 201 genes in the list and 1 of them is considered the disease gene

library(ggplot2)

# Date of data
date <- "2019-12-10"
outdir_name <- paste0(dirname(rstudioapi::getActiveDocumentContext()$path),"/../Plots/",date)


distribution <- NULL
for(k in c(1:10)){
  gene <- NULL
  success <- NULL
  top10 <- list()
  
  for(j in c(1:97)){
    gene[j] <- sample(201, 1)
  }
  
  for(j in c(1:1000)){
    top10[[j]] <- sample(201, 10, replace = FALSE)
  }
  
  success <- unlist(lapply(top10, function(x) gene %in% x))
  
  distribution <- c(distribution, mean(success))
}

distribution_df <- data.frame("nr_success" = distribution)

steps <- unique(DT_per_parameter$Step)
threshs <- unique(DT_per_parameter$Z_threshold)
maxrxns <- unique(DT_per_parameter$Max_rxn)

p <- ggplot(distribution_df) + 
  geom_density(aes(x = nr_success, y=..scaled.., colour = "Random"))
for(st in steps){
  for(thr in threshs){
    for(mx in maxrxns){
      unique_values <- unique(DT_per_parameter[Validation == FALSE & Step == 0 & Z_threshold == thr & Max_rxn == mx, Prior.frac10])
      if(length(unique_values) == 1){
        next()
      }
      p <- p + geom_density(data = DT_per_parameter[Validation == FALSE & Step == st & Z_threshold == thr & Max_rxn == mx, ], 
                            aes(x = Prior.frac10, y=..scaled.., colour = "Cross-omics"))
    }
  }
}


p <- p + theme_light() +
  scale_x_continuous("Ratio in top 10", breaks = seq(0,0.7,0.1), limits = c(0, 0.7)) +
  ylab("Distribution") +
  labs(color = "Method") 
  # ggtitle("Comparison performance random vs Cross-omics")

ggsave(paste0(outdir_name,"/Cross-omics_vs_random_performance.png"), plot = p,
       width = 200, height = 50, dpi=300, units = "mm")
  
  # geom_vline(xintercept = DT_per_parameter[Validation == FALSE, Av.prior.frac10], colour = 'red') +

