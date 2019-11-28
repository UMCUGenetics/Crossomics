# Little script to check the status of all gene-metabolite sets

setwd("/Users/mkerkho7/DIMS2_repo/Crossomics/Data/2019-08-12/")

maxrxns <- c(8, 10, 12, 15, 17, 19)
steps <- c(0, 1, 2, 3, 4, 5)

outdir_name <- paste0(dirname(rstudioapi::getActiveDocumentContext()$path),"/../Plots/",Sys.Date())
if (!file.exists(outdir_name))dir.create(outdir_name, recursive = TRUE)

gn_fi_list <- list.files("../mss/")
gn_fi_list <- sub(".RData",".RDS", gn_fi_list)


library(data.table)
library(ggplot2)

nr_dt_rows <- length(maxrxns)*length(steps)*length(gn_fi_list)

char_vec <- rep("NA", nr_dt_rows)
int_vec <- rep(0L, nr_dt_rows)

DT <- data.table(Gene = char_vec,
                 Nr_mets = int_vec,
                 Step = int_vec,
                 Maxrxn = int_vec)

i <- 0

for(mx in maxrxns){
  for(st in steps){
    subdir <- paste0("maxrxn",mx,"/mss_",st,"_HMDBtranslated/")
    for(gn in gn_fi_list){
      i <- i + 1
      met_set <- tryCatch(readRDS(paste0(subdir,gn)), error = function(e) paste("file",gn,"not found for Step",st,"Maxrxn",mx,"\n"))
      hmdb_codes <- met_set[grep("HMDB", met_set[,"hmdb"]),"hmdb"]
      nr_unq_mets <- uniqueN(hmdb_codes)
      DT[i, c("Gene","Nr_mets","Step","Maxrxn") := list(gn, nr_unq_mets, st, mx)]
    }
  }
}
DT[, Step := as.factor(Step)]
DT[, Maxrxn := as.factor(Maxrxn)]
DT[, Nr_mets_log := log10(Nr_mets + 1)]


# DT_par <- data.table(Median = rep(0L, length(maxrxns)*length(steps)),
#                      Min_ran = rep(0L, length(maxrxns)*length(steps)),
#                      Max_ran = rep(0L, length(maxrxns)*length(steps)))
# DT[row_nr, c("Median","Min_ran","Max_ran", "Step", "Max_rxn") := list(
#   DT[, median(Nr_mets), by = .(Step, Maxrxn)], patient, paste(dis_gene,collapse = ";"), dataset, train_val_set, dbs,  as.integer(set_size), as.character(prot_func), as.character(is_trans), 
#                rank, gene_set, 1-((rank-1)/(gene_set-1)), rank/gene_set, threshs, step, maxrxn,
#                as.double(p.value), seed)]

# DT[, median(Nr_mets), by = Step ]
# DT[, range(Nr_mets), by = .(Step, Maxrxn)]


DT[, Empty := sum(Nr_mets) == 0, by = Gene]
empty_genes <- unique(DT[Empty == TRUE, Gene])

DT <- DT[ Empty == FALSE, ]



q <- position_dodge(0.8)
p <- ggplot(DT, aes(y = Nr_mets_log, x = Step, fill = Maxrxn)) +
  geom_point(aes(fill = Maxrxn), size = 0.5, shape = 21, position = position_jitterdodge(), alpha = 0.3) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA) +
  # stat_summary(geom = "crossbar", width = 0.6, fatten=0, color="red", position = q,
  #              fun.data = function(x){c(y=median(x), ymin=median(x), ymax=median(x))}) +
  scale_y_continuous(breaks = round(seq(min(DT$Nr_mets_log), ceiling(max(DT$Nr_mets_log)), by = 0.5),1)) +
  labs(y = "Log10(Nr_mets + 1)", x = "Distance to primary reaction", fill = "Extension stringency")
  # scale_fill_gradient(name = "n occurrences", trans="pseudo_log")

ggsave(paste0(outdir_name,"/Met_set_overview.png"), plot = p,
       width = 300, height = 100, dpi=100, units = "mm")


