setwd("Z:/Users/Daniel/Support/scripts/R_code/SCDE2018")
# linux path: setwd("/workspace/ws03/IData/Users/Daniel/Support/scripts/R_code/SCDE2018/")
library("scde")
##### Get parameters #####
Sys.setlocale("LC_COLLATE", "C")
parameters<-read.table("input/parameters.txt", sep="\t", header=T, row.names=1, as.is=T);
parameterNames<-row.names(parameters)
combineAllGroupsPrior<-as.logical(parameters["CombineAllGroupsPriorBool", 1])
cleanCounts<-as.logical(parameters["CleanCountsBool", 1])
groupstring<-as.character(parameters["Group", 1])
batchstring<-as.character(parameters["Batch", 1])
cleanCounts<-as.logical(parameters["CleanCountsBool", 1])
cpuNumber<-strtoi(parameters["CpuNumber", 1])
minLibSize<-strtoi(parameters["MinLibSize", 1])
minGeneReadCount<-strtoi(parameters["MinGeneReadCount", 1])
minGeneCellCount<-strtoi(parameters["MinGeneCellCount", 1])
linearFit<-as.logical(parameters["LinearFitBool", 1])
threshouldSegmentation<-as.logical(parameters["ThreshouldSegmentationBool", 1])
minNoFail<-strtoi(parameters["MinNoFail", 1])
minCountThreshold<-strtoi(parameters["MinCountThreshold", 1])
zeroCountThreshold<-strtoi(parameters["ZeroCountThreshold", 1])
zeroLambda<-as.numeric(parameters["ZeroLambda", 1])
minSizeEntries<-strtoi(parameters["MinSizeEntries", 1])
maxPairs<-strtoi(parameters["MaxPairs", 1])
minPairsPerCell<-strtoi(parameters["MinPairsPerCell", 1])
localThetaFit<-as.logical(parameters["LocalThetaFitBool", 1])
filterUnCorrelatedCells<-as.logical(parameters["FilterUnCorrelatedCellsBool", 1])
gePriorGridResolution<-strtoi(parameters["GEPriorGridResolution", 1])
diffNumberRandomizations<-strtoi(parameters["DiffNumberRandomizations", 1])
caseGroup <- as.character(parameters["Compare", 1])
controlGroup <- as.character(parameters["CompareTo", 1])
### load design table ###
designTable<-read.table("input/design.txt", sep="\t", header=T, row.names=1, na.strings="", quote="", comment.char=""); # this potiential will change column names
groupinput<- NULL
if(groupstring!=""){
  groupinput <- factor(designTable[, groupstring])
}
batchinput<-NULL
if(batchstring!=""){
  batchinput <- factor(designTable[, batchstring])
}
### Load data matrix ###
datamatrix<-as.matrix(read.table("input/datamatrix.txt", sep="\t", header=T, row.names=1, na.strings = "."));
### clean up the dataset ###
if(cleanCounts){
  cd <- clean.counts(datamatrix, min.lib.size=minLibSize, min.reads = minGeneReadCount, min.detected = minGeneCellCount)
  allcells<-colnames(datamatrix)
  cleancells<-colnames(cd)
  cellindices<-match(cleancells, allcells)
  groupinput<-groupinput[cellindices]
} else {
  cd<-datamatrix
}
# calculate models
priorgroup<-NULL
if(!combineAllGroupsPrior){
  priorgroup <- groupinput
}
o.ifm <- scde.error.models(counts = cd, groups = priorgroup, n.cores = cpuNumber, threshold.segmentation = threshouldSegmentation, save.crossfit.plots = FALSE, save.model.plots = FALSE, verbose = 0,
                           linear.fit=linearFit, min.nonfailed=minNoFail, min.count.threshold=minCountThreshold, zero.count.threshold=zeroCountThreshold, zero.lambda=zeroLambda,min.size.entries= minSizeEntries,
                           max.pairs=maxPairs, min.pairs.per.cell=minPairsPerCell,local.theta.fit=localThetaFit)
# filter out cells that don't show positive correlation with
# the expected expression magnitudes (very poor fits)
if(filterUnCorrelatedCells){
  valid.cells <- o.ifm$corr.a > 0
  o.ifm <- o.ifm[valid.cells, ]
}
# estimate gene expression prior
o.prior <- scde.expression.prior(models = o.ifm, counts = cd, length.out = gePriorGridResolution, show.plot = FALSE)
# run differential expression tests on all genes.
ediff <- scde.expression.difference(o.ifm, cd, o.prior, groups  =  groupinput, batch = batchinput, n.randomizations  =  diffNumberRandomizations, n.cores  =  cpuNumber, verbose  =  0)
p.values <- 2*pnorm(abs(ediff$Z),lower.tail=F) # 2-tailed p-values
p.values.adj <- 2*pnorm(abs(ediff$cZ),lower.tail=F) # Adjusted to control for FDR
lblfc <- ediff$lb
lfc <- ediff$mle
ublfc <- ediff$ub
if(levels(groupinput)[1]==controlGroup){
  ### reverse the comparison direction
  lblfc<- -ediff$ub
  lfc<- -ediff$mle
  ublfc<- -ediff$lb
}
de <- cbind(lblfc,lfc,ublfc, p.values, p.values.adj)
colnames(de) <- c("LB","Estimate[log2FC]","UB","RawPvalue", "AdjPvalue")
rownames(de) <- rownames(ediff)
# get expression magntiude estimates
o.fpm <- scde.expression.magnitude(o.ifm, counts = cd)
### export results ###
write.table(de, file = "output/SCDE.txt", row.names = T, col.names = NA, sep = "\t", quote = F)
write.table(o.fpm, file = "output/fitExpression.txt", sep="\t", quote=F, col.names=NA, row.names=T)
