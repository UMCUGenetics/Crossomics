library("epade")

base_col="light yellow"
red="red"

library("XLConnect")
# ls("package:XLConnect")

# PAH, CBS
wb = loadWorkbook("./results/crossomics_SP_PAH_CBS_step1_not_normalized_1/MSEA_overview.xls")
data = readWorksheet(wb, sheet = 1)

colscheme = rbind(c(red,rep(base_col,9)),
                  c(rep(base_col,1),red,rep(base_col,8)),
                  c(rep(base_col,2),red,rep(base_col,7)),
                  c(rep(base_col,1),red,rep(base_col,8)),
                  c(red,rep(base_col,9)),
                  c(red,rep(base_col,9)))

index=c(1:(dim(data)[1]/3))*3-1
tmp=matrix((as.numeric(as.matrix(data[index,c(1:10)]))), nrow = length(index), ncol = 10)
genes=c(rep("PAH",5),"CBS")

wb = loadWorkbook("./results/crossomics_SP_MAT1A_MTHFR_step1_not_normalized_2/MSEA_overview.xls")
data = readWorksheet(wb, sheet = 1)

colscheme = rbind(colscheme,
                  c(rep(base_col,4),red,rep(base_col,5)),
                  c(red,rep(base_col,9)),
                  c(rep(base_col,2),red,rep(base_col,7)),
                  c(rep(base_col,2),red,rep(base_col,7)))

index=c(1:(dim(data)[1]/3))*3-1
tmp=rbind(tmp, matrix((as.numeric(as.matrix(data[index,c(1:10)]))), nrow = length(index), ncol = 10))
genes=c(genes, rep("MAT1A",2),rep("MTHFR",2))

wb = loadWorkbook("./results/crossomics_SP_GCDH_IVD_step1_not_normalized_2/MSEA_overview.xls")
data = readWorksheet(wb, sheet = 1)

colscheme = rbind(colscheme,
                  c(red,rep(base_col,9)),
                  c(rep(base_col,1),red,rep(base_col,8)),
                  c(rep(base_col,1),red,rep(base_col,8)),
                  c(red,rep(base_col,9)))

index=c(1:(dim(data)[1]/3))*3-1
tmp=rbind(tmp, matrix((as.numeric(as.matrix(data[index,c(1:10)]))), nrow = length(index), ncol = 10))
genes=c(genes, rep("GCDH",2),rep("IVD",2))

wb = loadWorkbook("./results/crossomics_SP_MCT1_step1_normalized_1/MSEA_overview.xls")
data = readWorksheet(wb, sheet = 1)

colscheme = rbind(colscheme,
                  c(red,rep(base_col,9)))

index=c(1:(dim(data)[1]/3))*3-1
tmp=rbind(tmp, matrix((as.numeric(as.matrix(data[index,c(1:10)]))), nrow = length(index), ncol = 10))
genes=c(genes, "SLC16A1")

wb = loadWorkbook("./results/crossomics_SPIV_more_noise/MSEA_overview.xls")
data = readWorksheet(wb, sheet = 1)

colscheme = rbind(colscheme,
                  c(red,red,rep(base_col,8)),
                  c(red,red,red,rep(base_col,7)),
                  c(red,red,red,rep(base_col,7)),
                  c(rep(base_col,3),red,rep(base_col,6)),
                  c(rep(base_col,3),red,rep(base_col,6)),
                  c(rep(base_col,4),red,rep(base_col,5)),
                  c(rep(base_col,1),red,rep(base_col,8)),
                  c(rep(base_col,1),red,rep(base_col,8)),
                  c(rep(base_col,2),red,rep(base_col,7)),
                  c(rep(base_col,1),red,rep(base_col,8)),
                  c(rep(base_col,1),red,rep(base_col,8)),
                  c(rep(base_col,3),red,rep(base_col,6)),
                  c(rep(base_col,4),red,rep(base_col,5)),
                  c(rep(base_col,2),red,rep(base_col,7)),
                  c(red,rep(base_col,9)),
                  c(rep(base_col,2),red,red,rep(base_col,6)),
                  c(red,rep(base_col,9)),
                  c(rep(base_col,1),red,rep(base_col,8)))

index=c(1:(dim(data)[1]/3))*3-1
tmp=rbind(tmp, matrix((as.numeric(as.matrix(data[index,c(1:10)]))), nrow = length(index), ncol = 10))
genes=c(genes, "MCCC1","AMT","AMT","CYP27A1","CYP27A1","DHCR7","DHCR7","FAH","ACADVL","ACADM","ACADM","MUT","MUT","OTC","OTC","ABCG5","SLC7A7","CPT1A")


####################
tmp = genes = colscheme = NULL
wb = loadWorkbook("./results/crossomics/MSEA_overview.xls")
data = readWorksheet(wb, sheet = 1)

colscheme = rbind(colscheme,
                  c(red,rep(base_col,3),red,rep(base_col,5)),
                  c(base_col,red,red,red,rep(base_col,6)),
                  c(red,red,red,rep(base_col,7)),
                  c(red,rep(base_col,9)),
                  c(red,rep(base_col,9)),
                  c(rep(base_col,2),red,rep(base_col,7)),
                  c(rep(base_col,1),red,rep(base_col,8)),
                  c(rep(base_col,2),red,rep(base_col,7)),
                  c(rep(base_col,1),red,rep(base_col,8)),
                  c(red,rep(base_col,9)),
                  c(red,rep(base_col,9)),
                  c(base_col,red,rep(base_col,8)),
                  c(rep(base_col,4),red,rep(base_col,5)),
                  c(rep(base_col,1),red,rep(base_col,8)),
                  c(rep(base_col,3),red,rep(base_col,6)),
                  c(rep(base_col,2),red,red,rep(base_col,6)),
                  c(red,rep(base_col,9)),
                  c(red,rep(base_col,9)))

index=c(1:(dim(data)[1]/3))*3-1
tmp=rbind(tmp, matrix((as.numeric(as.matrix(data[index,c(1:10)]))), nrow = length(index), ncol = 10))
genes=c(genes, "MCCC1","AMT","AMT","CYP27A1","CYP27A1","DHCR7","DHCR7","FAH","ACADVL","ACADM","ACADM","MUT","MUT","OTC","OTC","ABCG5","SLC7A7","CPT1A")
###################

tmp = genes = colscheme = NULL
wb = loadWorkbook("./results/crossomics_Fisher_weighted_step=1_-1=Z=1.5_a=1.2_whole_set/MSEA_overview_all_ranked.xls")
# wb = loadWorkbook("./results/crossomics_new_sets_step_1_thresh_-1.5,2_metabolite_sets_1.0_filter_1.1_manuscript/MSEA_overview_all.xls")
# wb = loadWorkbook("./results/crossomics/MSEA_overview_ALL.xls")
data = readWorksheet(wb, sheet = 1)

dim(data)
data=data[,1:28]
# data=data[,1:30]

# 14 X 1
# 9 x 2
# 4 x 3
# 4 x 4
# 1 x 7
# 1 X 8

# colscheme = rbind(c(red, rep(base_col,27)),#PAH 1
#                   c(base_col, red, rep(base_col,26)),#PAH 2
#                   c(red, rep(base_col,27)),#PAH 1
#                   c(red, rep(base_col,27)),#PAH 1
#                   c(base_col, red, rep(base_col,26)),#PAH 2
#                   c(red, rep(base_col,27)),#CBS 1
#                   c(base_col, red, rep(base_col,26)),#MAT1A 2
#                   c(base_col, red, rep(base_col,26)),#MAT1A 2
#                   c(rep(base_col,2), red, rep(base_col,25)),#MTHFR 3
#                   c(base_col, red, rep(base_col,26)),#MTHFR 2
#                   c(base_col, red, rep(base_col,26)),#GCDH 2
#                   c(rep(base_col,15), red, rep(base_col,12)),#GCDH 16
#                   c(rep(base_col,10), red, rep(base_col,17)),#IVD 11
#                   c(rep(base_col,20), red, rep(base_col,7)),#IVD 21
#                   c(rep(base_col,18), red, rep(base_col,9)),#SLC16A1 19
#                   c(rep(base_col,7), red, rep(base_col,20)), # MCCC1 8
#                   c(rep(base_col,2), red, rep(base_col,25)), # "AMT" 3
#                   c(base_col, red, rep(base_col,26)), # "AMT" 2
#                   c(base_col, red, rep(base_col,26)), # "CYP27A1" 2
#                   c(red, rep(base_col,27)), # ,"CYP27A1" 1
#                   c(rep(base_col,3), red, rep(base_col,24)), # ,"DHCR7" 4
#                   c(rep(base_col,2), red, rep(base_col,25)), # ,"DHCR7" 3
#                   c(rep(base_col,9), red, rep(base_col,18)), # ,"FAH" 10
#                   c(rep(base_col,2), red, rep(base_col,25)), # ,"ACADVL" 3
#                   c(red, rep(base_col,27)), # ,"ACADM" 1
#                   c(red, rep(base_col,27)), # ,"ACADM" 1
#                   c(rep(base_col,27), red), # ,"MUT" 100
#                   c(rep(base_col,27), red), # "MUT" 100
#                   c(base_col, red, rep(base_col,26)), # "OTC" 2
#                   c(red, rep(base_col,27)), # ,"OTC" 1
#                   c(rep(base_col,27), red), # ,"ABCG5" 100
#                   c(rep(base_col,27), red), # ,"AGA" 100
#                   c(rep(base_col,27), red), # ,"SLC7A7" 100
#                   c(base_col, red, rep(base_col,26))) # ,"CPT1A" 2


# colscheme = rbind(c(base_col, red, rep(base_col,28)),#PAH 2
#                   c(base_col, red, rep(base_col,28)),#PAH 2
#                   c(red, rep(base_col,29)),#PAH 1
#                   c(red, rep(base_col,29)),#PAH 1
#                   c(red, rep(base_col,29)),#PAH 1
#                   c(red, rep(base_col,29)),#CBS 1
#                   c(red, rep(base_col,29)),#MAT1A 1
#                   c(red, rep(base_col,29)),#MAT1A 1
#                   c(rep(base_col,3), red,rep(base_col,26)),#MTHFR 4
#                   c(rep(base_col,2), red, rep(base_col,27)),#MTHFR 3
#                   c(base_col, red, rep(base_col,28)),#GCDH 2
#                   c(base_col, red, rep(base_col,28)),#GCDH 2
#                   c(red, rep(base_col,29)),#IVD 1
#                   c(red, rep(base_col,29)),#IVD 1
#                   c(base_col, red, rep(base_col,28)),#SLC16A1 2
#                   c(rep(base_col,6), red, rep(base_col,23)), # MCCC1 7
#                   c(red, rep(base_col,29)), # "GCSH" 1
#                   c(base_col, red, rep(base_col,28)), # "GCSH" 2
#                   c(red, rep(base_col,29)), # "CYP27A1" 1
#                   c(rep(base_col,2), red, rep(base_col,27)), # ,"CYP27A1" 3
#                   c(base_col, red, rep(base_col,28)), # ,"DHCR7" 2
#                   c(red, rep(base_col,29)), # ,"DHCR7" 1
#                   c(base_col, red, rep(base_col,28)), # ,"FAH" 2
#                   c(base_col, red, rep(base_col,28)), # ,"ACADVL" 2
#                   c(red, rep(base_col,29)), # ,"ACADM" 1
#                   c(red, rep(base_col,29)), # ,"ACADM" 1
#                   c(rep(base_col,3), red, rep(base_col,26)), # ,"MUT" 4
#                   c(rep(base_col,4), red, rep(base_col,25)), # "MUT" 5
#                   c(rep(base_col,4), red, rep(base_col,25)), # "OTC" 5
#                   c(red, rep(base_col,29)), # ,"OTC" 1
#                   c(rep(base_col,7), red, rep(base_col,22)), # ,"ABCG5" 8
#                   c(base_col, red, rep(base_col,28)), # ,"SLC7A7" 2
#                   c(red, rep(base_col,29))) # ,"CPT1A" 1

colscheme = rbind(c(red, rep(base_col,29)),
                  c(red, rep(base_col,29)),
                  c(red, rep(base_col,29)),
                  c(red, rep(base_col,29)),
                  c(red, rep(base_col,29)),
                  c(red, rep(base_col,29)),
                  c(red, rep(base_col,29)),
                  c(red, rep(base_col,29)),
                  c(red, rep(base_col,29)),
                  c(red, rep(base_col,29)),
                  c(red, rep(base_col,29)),
                  c(red, rep(base_col,29)),
                  c(red, rep(base_col,29)),
                  c(red, rep(base_col,29)),
                  c(red, rep(base_col,29)),
                  c(base_col, red, rep(base_col,28)),
                  c(base_col, red, rep(base_col,28)),
                  c(base_col, red, rep(base_col,28)),
                  c(base_col, red, rep(base_col,28)),
                  c(base_col, red, rep(base_col,28)),
                  c(base_col, red, rep(base_col,28)),
                  c(base_col, red, rep(base_col,28)),
                  c(base_col, red, rep(base_col,28)),
                  c(base_col, red, rep(base_col,28)),
                  c(base_col, red, rep(base_col,28)),
                  c(rep(base_col,2), red, rep(base_col,27)),
                  c(rep(base_col,2), red, rep(base_col,27)),
                  c(rep(base_col,3), red,rep(base_col,26)),
                  c(rep(base_col,3), red,rep(base_col,26)),
                  c(rep(base_col,4), red, rep(base_col,25)),
                  c(rep(base_col,4), red, rep(base_col,25)),
                  c(rep(base_col,6), red, rep(base_col,23)),
                  c(rep(base_col,7), red, rep(base_col,22)))


index=c(1:((dim(data)[1]+1)/3))*3-2
tmp=matrix((as.numeric(as.matrix(data[index,c(1:28)]))), nrow = length(index), ncol = 28)
# genes1=c(rep("PAH",5),"CBS", rep("MAT1A",2),rep("MTHFR",2),rep("GCDH",2),rep("IVD",2),"SLC16A1","MCCC1","GCSH","GCSH","CYP27A1","CYP27A1","DHCR7","DHCR7","FAH","ACADVL","ACADM","ACADM","MUT","MUT","OTC","OTC","ABCG5","SLC7A7","CPT1A")

# tmp=matrix((as.numeric(as.matrix(data[index,c(1:28)]))), nrow = length(index), ncol = 28)
# genes1=c(rep("PAH",5),"CBS", rep("MAT1A",2),rep("MTHFR",2),rep("GCDH",2),rep("IVD",2),"SLC16A1","MCCC1","GCSH","GCSH","CYP27A1","CYP27A1","DHCR7","DHCR7","FAH","ACADVL","ACADM","ACADM","MUT","MUT","OTC","OTC","ABCG5","AGA","SLC7A7","CPT1A")

tmp[which(tmp==0)]=NA

genes=paste0("P",c(3,4,5,6,7,8,13,14,17,19,22,25,26,30,33,1,2,11,12,15,18,21,23,24,32,10,20,9,27,28,29,16,31))


# genes=paste0("P",c(1:34))

bar3d.ade(t(tmp), wall=1, xlab="patients", ylab="p-value", zlab="rank", col=t(colscheme), main = "Enrichment results for training set",
          xticks=genes, xw=0.8, zw=0.5,
          zticks = c(1:28))
# bar3d.ade(t(tmp), wall=1, xlab="patients", ylab="p-value", zlab="rank", col=t(colscheme), main = "Cross-omics results for test data set",
#           xticks=genes, xw=0.8, zw=0.5,
#           zticks = c(1:28))


#, alpha=1
# , xw=1, zw=0.2, axes=FALSE
cbind(genes,genes1)

#############################
