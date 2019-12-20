# Little script to check the status of all gene-metabolite sets

setwd("/Users/mkerkho7/DIMS2_repo/Crossomics/Data/2019-08-12/")

maxrxns <- c(8, 10, 12, 15, 17, 19)
steps <- c(0, 1, 2, 3, 4, 5)

date <- "2019-12-10"

outdir_name <- paste0(dirname(rstudioapi::getActiveDocumentContext()$path),"/../Plots/",date)
if (!file.exists(outdir_name))dir.create(outdir_name, recursive = TRUE)

gn_fi_list <- list.files("../mss/")
gn_fi_list <- sub(".RData",".RDS", gn_fi_list)


library(data.table)
library(ggplot2)

nr_dt_rows <- length(maxrxns)*length(steps)*length(gn_fi_list)

char_vec <- rep("NA", nr_dt_rows)
int_vec <- rep(0L, nr_dt_rows)

DT_met_sets <- data.table(Gene = char_vec,
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
      DT_met_sets[i, c("Gene","Nr_mets","Step","Maxrxn") := list(gn, nr_unq_mets, st, mx)]
    }
  }
}
DT_met_sets[, Step := as.factor(Step)]
DT_met_sets[, Maxrxn := as.factor(Maxrxn)]
DT_met_sets[, Nr_mets_log := log2(Nr_mets + 1)]

DT_met_sets[, Empty := sum(Nr_mets) == 0, by = Gene]
empty_genes <- unique(DT_met_sets[Empty == TRUE, Gene])

DT_met_sets <- DT_met_sets[ Empty == FALSE, ]
DT_met_sets[, Gene := str_replace(Gene, pattern = ".RDS", replacement = "")]

saveRDS(DT_met_sets, file = paste0(code_dir,"/../Results/",date,"/Metabolite_set_sizes.RDS"))

q <- position_dodge(0.8)
p1 <- ggplot(DT_met_sets, aes(y = Nr_mets_log, x = Step, fill = Maxrxn)) +
  geom_point(aes(fill = Maxrxn), size = 0.5, shape = 21, position = position_jitterdodge(), alpha = 0.3) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA) +
  # stat_summary(geom = "crossbar", width = 0.6, fatten=0, color="red", position = q,
  #              fun.data = function(x){c(y=median(x), ymin=median(x), ymax=median(x))}) +
  scale_y_continuous(breaks = round(seq(min(DT_met_sets$Nr_mets_log), ceiling(max(DT_met_sets$Nr_mets_log)), by = 1),1)) +
  # labs(y = "Log2(Nr_mets + 1)", x = "Distance to primary reaction", fill = "Extension stringency")
  labs(y = bquote(~Log[2]~"([mets in set] + 1)"), x = "Distance to primary reaction", fill = "Extension stringency")
  # scale_fill_gradient(name = "n occurrences", trans="pseudo_log")

ggsave(paste0(outdir_name,"/Met_set_overview.png"), plot = p,
       width = 300, height = 100, dpi=100, units = "mm")

# And for the disease genes (need to get a vector of the genes themselves somewhere else)
p2 <- ggplot(DT_met_sets[Gene %in% disease_genes], aes(y = Nr_mets_log, x = Step, fill = Maxrxn)) +
  geom_point(aes(fill = Maxrxn), size = 0.5, shape = 21, position = position_jitterdodge(), alpha = 0.3) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA) +
  # stat_summary(geom = "crossbar", width = 0.6, fatten=0, color="red", position = q,
  #              fun.data = function(x){c(y=median(x), ymin=median(x), ymax=median(x))}) +
  scale_y_continuous(breaks = round(seq(min(DT_met_sets[Gene %in% disease_genes]$Nr_mets_log), ceiling(max(DT_met_sets[Gene %in% disease_genes]$Nr_mets_log)), by = 1),1)) +
  # labs(y = "Log2(Nr_mets + 1)", x = "Distance to primary reaction", fill = "Extension stringency")
  labs(y = bquote(~Log[2]~"([mets in set] + 1)"), x = "Distance to primary reaction", fill = "Extension stringency")
# scale_fill_gradient(name = "n occurrences", trans="pseudo_log")

facet_wrap(c(p1, p2))

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}
mylegend<-g_legend(p1)

gridExtra::grid.arrange(p1, p2, nrow = 2)

p3 <- gridExtra::grid.arrange(arrangeGrob(p1 + theme(legend.position="none"),
                                          p2 + theme(legend.position="none"),
                                          nrow=2),
                              mylegend, nrow=2,heights=c(10, 1))
ggsave(paste0(outdir_name,"/Met_set_overview_disease_genes.png"), plot = p2,
       width = 300, height = 100, dpi=100, units = "mm")


