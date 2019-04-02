###################################################
### chunk number 1: options
###################################################
options(width=50, digits=2)


###################################################
### chunk number 2: ALL
###################################################
library("ALL")
data("ALL")


###################################################
### chunk number 3: bcrabl
###################################################
bcell = grep("^B", as.character(ALL$BT)) # B-cell ALL
moltyp = which(as.character(ALL$mol.biol) 
               %in% c("NEG", "BCR/ABL"))
bcrneg = ALL[, intersect(bcell, moltyp)]
bcrneg$mol.biol = factor(bcrneg$mol.biol) # drop unused levels


###################################################
### chunk number 4: bcrabl-result
###################################################
dim(bcrneg)


###################################################
### chunk number 5: nsFilter
###################################################
library("genefilter")
bcrneg_filt1 = nsFilter(bcrneg, var.cutoff=0.5)$eset
dim(bcrneg_filt1)


###################################################
### chunk number 6: gsc-filter
###################################################
library(GSEABase)
gsc <- GeneSetCollection(bcrneg_filt1,
                         setType=KEGGCollection())


###################################################
### chunk number 7: gsc-egs
###################################################
gsc
gsc[[2]]


###################################################
### chunk number 8: gsc-filter-2
###################################################
ok <- sapply(geneIds(gsc), length)>10
gsc <- gsc[ok]
length(gsc)
uids <- unique(unlist(geneIds(gsc)))
bcrneg_filt2 <- bcrneg_filt1[uids,]
dim(bcrneg_filt2)


###################################################
### chunk number 9: rowtttest
###################################################
rtt <- rowttests(bcrneg_filt2, "mol.biol")
rttStat <- rtt$statistic
names(rttStat) <- featureNames(bcrneg_filt2)
head(rttStat)


###################################################
### chunk number 10: zstats
###################################################
zCalc <- function(ids, tStat) {
  sum(tStat[ids]) / sqrt(length(ids))
}
z <- sapply(geneIds(gsc), zCalc, tStat=rttStat)
names(z) <- names(gsc)
head(z)


###################################################
### chunk number 11: qqnorm
###################################################
qqnorm(z)
qqline(z)


###################################################
### chunk number 12: outlier
###################################################
z[which.min(z)]


###################################################
### chunk number 13: kegg-id
###################################################
keggId <- names(z[which.min(z)])
keggGS <- gsc[[keggId]]
keggES <- bcrneg_filt2[keggGS,]
library(Category)
getPathNames(keggId)


###################################################
### chunk number 14: keggInfo
###################################################
KEGGmnplot(keggId,
           bcrneg_filt2, 
           annotation(bcrneg_filt2),
           bcrneg_filt2$mol.biol, 
           pch=16, col="darkblue")


###################################################
### chunk number 15: olap
###################################################
overlap <- gsc[["04512"]] & gsc[["04510"]]
length(GSEABase::geneIds(overlap))


