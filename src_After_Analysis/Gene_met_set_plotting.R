# Plotting script for the output of Gene_met_set_stats.R

library(ggplot2)
library(data.table)
library(rstudioapi)



date <- "2019-12-10"

code_dir <- paste0(dirname(rstudioapi::getActiveDocumentContext()$path),"/../Results/")
outdir_name <- paste0(dirname(rstudioapi::getActiveDocumentContext()$path),"/../Plots/",date)
if (!file.exists(outdir_name))dir.create(outdir_name, recursive = TRUE)

DT_met_sets <- readRDS(file = paste0(code_dir, date,"/Metabolite_set_sizes.RDS"))
disease_genes <- scan(paste0(code_dir,"/../Data/disease_genes.txt"), character(), quote = "", quiet = TRUE)


#### Special legend function
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}


#### Plotting

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

ggsave(paste0(outdir_name,"/Met_set_overview.svg"), plot = p1,
       width = 300, height = 100, units = "mm")

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

ggsave(paste0(outdir_name,"/Met_set_overview_disease_genes.svg"), plot = p2,
       width = 300, height = 100, units = "mm")

facet_wrap(c(p1, p2))


mylegend<-g_legend(p1)

gridExtra::grid.arrange(p1, p2, nrow = 2)

p3 <- gridExtra::grid.arrange(arrangeGrob(p1 + theme(legend.position="none"),
                                          p2 + theme(legend.position="none"),
                                          nrow=2),
                              mylegend, nrow=2,heights=c(10, 1))



