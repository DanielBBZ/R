#install DESeq2 package if it's not installed yet
#source("https://bioconductor.org/biocLite.R")
#biocLite("DESeq2")

library("DESeq2")
 
setwd("//192.168.1.174/ws03/IData/Users/Daniel/Support/TestCase/R/Benchmark")

### Load data matrix ###
datamatrix <- read.table("MicroArrayData.txt", sep="\t", header=T, row.names=1)

#---------------- test1. DESeq2 Analysis for modle with A + B (balanced design)------------------------#
### load design tables ###
designT1 <- read.csv("design1.txt", sep="\t", header=T, row.names=1)
designT1$A = factor(designT1$A,levels=c("a", "b"))
designT1$B = factor(designT1$B,levels=c("a", "b"))

dds1<-DESeqDataSetFromMatrix(countData=round(datamatrix), colData=designT1, design= ~ A + B)

result1 <- DESeq(dds1, minReplicatesForReplace=7, modelMatrixType="standard")

resultsNames(result1)

test1 <- results(result1, contrast=c("A", "b", "a"))

write.table(test1, "test1.txt", sep="\t", row.names=TRUE, col.names = TRUE, quote=FALSE)

#---------------- test2. DESeq2 Analysis for modle with A + B + A:B (balanced design) ------------------------#
### use the same data matrix with test1 ###
### use the same design table with test2 (balanced) ###
designT2 <- read.csv("design2.txt", sep="\t", header=T, row.names=1)
designT2$A = factor(designT2$A,levels=c("a", "b"))
designT2$B = factor(designT2$B,levels=c("a", "b"))

dds2<-DESeqDataSetFromMatrix(countData=round(datamatrix), colData=designT2, design= ~ A + B + A:B)

result2 <- DESeq(dds2, minReplicatesForReplace=7, modelMatrixType="standard")

resultsNames(result2)

test2A <- results(result2, contrast=c(0,0,1,0))
test2B <- results(result2, contrast=c(0,0,1,1))

write.table(test2A, "test2A.txt", sep="\t", row.names=TRUE, col.names = TRUE, quote=FALSE)
write.table(test2B, "test2B.txt", sep="\t", row.names=TRUE, col.names = TRUE, quote=FALSE)


#----------------test3: DESeq2 Analysis for modle with A + B + A:B (balanced design) ------------------------#
### use the same data matrix with test1 ###
### load design file ###
designT3 <- read.csv("design3.txt", sep="\t", header=T, row.names=1)
designT3$A = factor(designT3$A,levels=c("a", "b"))
designT3$B = factor(designT3$B,levels=c("a", "b", "c"))

dds3<-DESeqDataSetFromMatrix(countData=round(datamatrix), colData=designT3, design= ~ A + B + A:B)

result3 <- DESeq(dds3, minReplicatesForReplace=7, modelMatrixType="standard")

resultsNames(result3)

test3_A <- results(result3, contrast=c(0,0,1,0,0,0))
test3_B <- results(result3, contrast=c(0,0,0,1,0,1))

write.table(test3_A, "test3_A.txt", sep="\t", row.names=TRUE, col.names = TRUE, quote=FALSE)
write.table(test3_B, "test3_B.txt", sep="\t", row.names=TRUE, col.names = TRUE, quote=FALSE)

#----------------test4: DESeq2 Analysis for modle with A + B + A:B (unbalanced design) ------------------------#
### use the same data matrix with test1 ###
### load design file ###
designT4 <- read.csv("design4.txt", sep="\t", header=T, row.names=1)
designT4$A = factor(designT4$A,levels=c("a", "b"))
designT4$B = factor(designT4$B,levels=c("a", "b", "c"))

# make the names identical
#colnames(datamatrix) <- rownames(designT4)

dds4<-DESeqDataSetFromMatrix(countData=round(datamatrix), colData=designT4, design= ~ A + B + A:B)

result4 <- DESeq(dds4, minReplicatesForReplace=7, modelMatrixType="standard")

resultsNames(result4)

test4_A <- results(result4, contrast=c(0,0,1,-1,0,0))
test4_B <- results(result4, contrast=c(0,0,1,-1,1,-1))

write.table(test4_A, "test4_A_new.txt", sep="\t", row.names=TRUE, col.names = TRUE, quote=FALSE)
write.table(test4_B, "test4_B_new.txt", sep="\t", row.names=TRUE, col.names = TRUE, quote=FALSE)
