library("BridgeDbR")
library("plotrix")

path = "./results/metabolite_sets_step_0,1,2,3_1.0_filter_1.1/mss_0"
fileName = "./results/metabolite_sets_step_0,1,2,3_1.0_filter_1.1/dist_WG_step_0.RData"

files_done = list.files(path)
# files_todo = list.files("./results/metabolite_sets_step_0,1,2,3_1.0_filter_1.1/mss_WG_step_0")

# no_set = length(files_todo) - length(files_done)

dist = NULL

for (i in 1:length(files_done)){
  load(paste(path,files_done[i],sep="/"))
  
  metaboliteSet=result_mets_0 # <=================================!!!!!!!!!!!!!!!!!!!!!!!!

  # message( dim(as.data.frame(retVal$mets))[1] )
  
  gene_in=strsplit(files_done[i], split = "." , fixed=TRUE)[[1]][1]
  
  # A lot of missing HMDB identifiers!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (!is.null(metaboliteSet)){
    if (is.null(dim(metaboliteSet))){
      metaboliteSet=data.frame(t(metaboliteSet))
    }
  }

  ##################################################
  # temporarily work around to be fixed in findMetabolicEnvironment
  metaboliteSet = as.matrix(metaboliteSet)
  index = which(is.na(metaboliteSet[,"hmdb"]))
  if (length(index)>0) metaboliteSet[index,"hmdb"] = "character(0)"
  ##################################################

  index = which(metaboliteSet[,"hmdb"] == "character(0)")

  # pressent in BridgeDB?
  index.sub = which(metaboliteSet[index,"kegg"] != "character(0)")
  kegg_id = metaboliteSet[index[index.sub],"kegg"]

  mapper = loadDatabase("./BridgeDB/metabolites_20150717.bridge")
  hmdb = getSystemCode("HMDB")
  kegg = getSystemCode("KEGG Compound")

  if (length(kegg_id)>0){
    for (k in 1:length(kegg_id)){
      if (!is.null(unlist(map(mapper, kegg, kegg_id[k], hmdb)[1]))) metaboliteSet[index[index.sub[k]],"hmdb"] = unlist(map(mapper, kegg, kegg_id[k], hmdb)[1])
    }
  }

  index = which(metaboliteSet[,"hmdb"] == "character(0)")
  index.sub = which(metaboliteSet[index,"chebi"] != "character(0)")
  chebi_id = metaboliteSet[index[index.sub],"chebi"]

  chebi = getSystemCode("ChEBI")

  if (length(chebi_id)>0){
    for (k in 1:length(chebi_id)){
      if (!is.null(unlist(map(mapper, chebi, chebi_id[k], hmdb)[1]))) metaboliteSet[index[index.sub[k]],"hmdb"] = unlist(map(mapper, chebi, chebi_id[k], hmdb)[1])
    }
  }

  index = which(metaboliteSet[,"hmdb"] == "character(0)")
  if (length(index)>0) metaboliteSet = metaboliteSet[-index,,drop=FALSE]
   
  message(paste("gene: ", gene_in))
  message(paste("metaboliteSet: ", dim(metaboliteSet)[1]))
  
  # dist = c(dist,dim(metaboliteSet)[1])
  dist=c(dist, length(unique(metaboliteSet[,"hmdb"])))
}

save(dist, file=fileName)

############################################################################
width=1024
height=768
path="./results/metabolite_sets_step_0,1,2,3_1.0_filter_1.1/mss_1"
load("./results/metabolite_sets_step_0,1,2,3_1.0_filter_1.1/dist_WG_step_1.RData")

color = "red"

files_done = list.files(path)
# files_todo = list.files("./results/mss_WG_step_0")
# no_set = length(files_todo) - length(files_done)

# dist = c(rep(0, no_set),dist)

f = factor(dist)
t = table(factor(f))

max(as.numeric(names(t)))
# 334
head(as.vector(t))

x=as.numeric(names(t))
y=as.vector(t)

max = max((y[-1])) + 50
gap = y[1] - (max + 600)

y = ifelse(y > max, y - gap, y)
plot(x, y, col=color, pch=16, xlab= '# metabolites in set', ylab='# genes', yaxt="n", main= paste("# genes =",
    length(dist), "\n", "with set =", sum(y[-1])), type="l", lwd=2)

yat <- pretty(y)
yat <- yat[yat!=max]
ylab <- ifelse(yat>max, yat+gap, yat)
axis(2,at=yat, labels=ylab)
axis.break(2,max,style="slash") 
legend("topright", legend = c("step = 1"), pch=c(16), col=c(color))

############################################################################
############################################################################
############################################################################
############################################################################
dist=NULL
load("./results/metabolite_sets_step_0,1,2,3_1.0_filter_1.1/dist_WG_step_0.RData")
index=which(dist==0)
dist=dist[-index]
# files_done = list.files("./results/metabolite_sets_step_0,1,2,3_1.0_filter_1.1/mss_0")
# no_set = length(files_todo) - length(files_done)
# dist = c(rep(0, no_set),dist)
f = factor(dist)
t = table(factor(f))
x0=as.numeric(names(t))
y0=as.vector(t)

dist=NULL
load("./results/metabolite_sets_step_0,1,2,3_1.0_filter_1.1/dist_WG_step_1.RData")
index=which(dist==0)
dist=dist[-index]
# files_done = list.files("./results/metabolite_sets_step_0,1,2,3_1.0_filter_1.1/mss_1")
# no_set = length(files_todo) - length(files_done)
# dist = c(rep(0, no_set),dist)
f = factor(dist)
t = table(factor(f))
x1=as.numeric(names(t))
y1=as.vector(t)

dist=NULL
load("./results/metabolite_sets_step_0,1,2,3_1.0_filter_1.1/dist_WG_step_2.RData")
index=which(dist==0)
dist=dist[-index]
# files_done = list.files("./results/metabolite_sets_step_0,1,2,3_1.0_filter_1.1/mss_2")
# files_todo = list.files("./results/metabolite_sets_step_0,1,2,3_1.0_filter_1.1/mss_WG_step_0")
# no_set = length(files_todo) - length(files_done)
# dist = c(rep(0, no_set),dist)
f = factor(dist)
t = table(factor(f))
x2=as.numeric(names(t))
y2=as.vector(t)

dist=NULL
load("./results/metabolite_sets_step_0,1,2,3_1.0_filter_1.1/dist_WG_step_3.RData")
index=which(dist==0)
dist=dist[-index]
# files_done = list.files("./results/metabolite_sets_step_0,1,2,3_1.0_filter_1.1/mss_3")
# no_set = length(files_todo) - length(files_done)
# dist = c(rep(0, no_set),dist)
f = factor(dist)
t = table(factor(f))
x3=as.numeric(names(t))
y3=as.vector(t)


# max = max(c(y0[-1],y1[-1],y2[-1])) + 20
# gap = min(c(y0[1],y1[1],y2[1])) - (max + 300)
# 
# y0 = ifelse(y0 > max, y0 - gap, y0)
# y1 = ifelse(y1 > max, y1 - gap, y1)
# y2 = ifelse(y2 > max, y2 - gap, y2)

# plot(x2, y2, col='blue', pch=16, xlab= '# metabolites in set', ylab='# genes', yaxt="n",
#      main=paste("Total number of genes =",length(files_todo), "\n", "of which", sum(y0[-1]), "form metabolite sets"))
# 
# xlim=75
# ylim=400
# lines(x=c(0:xlim), y=rep(ylim, length(c(0:xlim))), col='black', lty="dashed")
# lines(x=rep(xlim, length(c(0:ylim))), y=c(0:ylim), col='black', lty="dashed")
# lines(x=c(0:xlim), y=rep(0, length(c(0:xlim))), col='black', lty="dashed")
# lines(x=rep(0, length(c(0:ylim))), y=c(0:ylim), col='black', lty="dashed")
# 
# points(x1, y1, col='red', pch=17)
# points(x0, y0, col='green', pch=18)
# 
# plot(x2, y2, col='blue', pch=16, xlab= '', ylab='', xlim=c(0,xlim), ylim=c(0,ylim), cex=1.2)
# points(x1, y1, col='red', pch=17, xlim=c(0,xlim), ylim=c(0,ylim), cex=1.2)
# points(x0, y0, col='green', pch=18, xlim=c(0,xlim), ylim=c(0,ylim), cex=1.2)

############################################################################################
library("Cairo")

# cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5

CairoPNG(filename = "mets_sets_distribution.png")
  plot(x0[1:25], y0[1:25], col='green', xlab= 'Set size ', ylab='Number of sets', type="l", lwd=2, cex.main=1.4, cex.lab=1.5,
     main="Distribution of metabolites over sets, set size range: 1-25", axes = FALSE, xlim = c(1,25))
  axis(side = 1, at = c(1,5,10,15,20,25), cex.axis=1.2)
  axis(side = 2, at = c(1,200,400,600,800,983), cex.axis=1.2)
  
  lines(x1[1:25], y1[1:25], col='red', lwd=2)
  lines(x2[1:25], y2[1:25], col='blue', lwd=2)
  lines(x3[1:25], y3[1:25], col='dark grey', lwd=2)
  
  legend("topright", legend = c("primary", "step = 1", "step = 2", "step = 3"), pch=16, col=c("green","red","blue","dark gray"), cex = 1.2)
dev.off()

CairoPNG(filename = "mets_sets_distribution_2.png")
  plot(x0, y0, col='green', xlab= 'Set size', ylab='Number of sets', type="l", lwd=2, cex.main=1.4, cex.lab=1.5,
       main="All sets (3737), maximum set size: 275", axes = FALSE,xlim = c(1,275))
  axis(side = 1, at = c(1,50,100,150,200,250,275), cex.axis=1.2)
  axis(side = 2, at = c(1,200,400,600,800,983), cex.axis=1.2)
  
  lines(x1, y1, col='red', lwd=2)
  lines(x2, y2, col='blue', lwd=2)
  lines(x3, y3, col='dark grey', lwd=2)
  
  legend("topright", legend = c("primary", "step = 1", "step = 2", "step = 3"), pch=16, col=c("green","red","blue","dark gray"), cex = 1.2)
dev.off()
 

# , pch=16, yaxt="n",xaxt="n"

# main=paste("Total number of genes =",length(dist), "\n", "of which", sum(y0[-1]), "form metabolite sets")


# # , type="l", lwd=2
# plot(x0, y0, col='green', pch=16, xlab= '# of metabolites in set', ylab='# of sets', sex=1, type="b", lwd=2, xlim=c(1,150),
#      main="Distribution of metabolites over sets")

# paste("Total number of coding genes = 19029", "\n", "of which", sum(y0), "form metabolite sets")


# lines(x1, y1, col='red', type="b", sex=1, pch=16,lwd=2)
# lines(x2, y2, col='blue', type="b",  sex=1, pch=16,lwd=2)
# lines(x3, y3, col='dark grey', type="b",  sex=1, pch=16,lwd=2)

# points(x1, y1, col='red', pch=16)
# points(x2, y2, col='blue', pch=16)
# points(x3, y3, col='dark grey', pch=16)

# lines(x1, y1, col='red', lwd=2)
# lines(x2, y2, col='blue', lwd=2)
# lines(x3, y3, col='dark grey', lwd=2)


# xlim=50
# ylim=860
# lines(x=c(0:xlim), y=rep(ylim, length(c(0:xlim))), col='black', lty="dashed")
# lines(x=rep(xlim, length( c(0:ylim))), y=c(0:ylim), col='black', lty="dashed")
# lines(x=c(0:xlim), y=rep(0, length(c(0:xlim))), col='black', lty="dashed")
# lines(x=rep(0, length(c(0:ylim))), y=c(0:ylim), col='black', lty="dashed")
# 
# 
# plot(x2[-1], y2[-1], col='blue', pch=16, xlab= '', ylab='', xlim=c(0,xlim), ylim=c(0,ylim), cex=1.2, type="l",lwd=2)
# lines(x1[-1], y1[-1], col='red', pch=17, xlim=c(0,xlim), ylim=c(0,ylim),lwd=2)
# lines(x0[-1], y0[-1], col='green', pch=18, xlim=c(0,xlim), ylim=c(0,ylim),lwd=2)
#######################################################################################

# yat <- pretty(y2)
# yat <- yat[yat!=max]
# ylab <- ifelse(yat>max, yat+gap, yat)
# axis(2,at=yat, labels=ylab)
# axis.break(2,max,style="slash") 
# legend("topright", legend = c("primary", "step = 1", "step = 2", "step = 3"), pch=c(16:18), col=c("green","red","blue","dark gray"))
