#install DESeq2 package if it's not installed yet
#source("https://bioconductor.org/biocLite.R")
#biocLite("DESeq2")

#Purpose; run DESeq in R and AS, compare the results, create scatter plot and show R2 and P value:

library("DESeq2")
library("ggplot2")

setwd("//192.168.1.174/ws03/IData/Users/Daniel/Test/TestProject/DESeq2")

### Load data matrix ###
datamatrix <- read.table("AbsoluteReads_final.txt", sep="\t", header=T, row.names=1)

#---------------- test1. DESeq2 Analysis for modle with Day only------------------------#
### load design tables ###
designT1 <- read.csv("AbsoluteReads_final_Design.txt", sep="\t", header=T, row.names=1)
designT1$Day = factor(designT1$Day,levels=c("0", "3"))

dds1<-DESeqDataSetFromMatrix(countData=round(datamatrix), colData=designT1, design= ~ Day)

result1 <- DESeq(dds1, minReplicatesForReplace=7, modelMatrixType="standard")

resultsNames(result1)

test1 <- results(result1, contrast=c("Day", "3", "0"))

write.table(test1, "test1.txt", sep="\t", row.names=TRUE, col.names = TRUE, quote=FALSE)

#------------------compare with AS result-------------------------------------#
### Load data matrix ###
ASresult <- read.table("ASresult.txt", sep="\t", header=T, row.names=1)

#generate fold change for R output
test1['FoldChange'] <- ifelse(test1$log2FoldChange >0, 2^test1$log2FoldChange, 0-2^(-test1$log2FoldChange))


#-------------------scatter plot--------------------------------------------#
test1result <- as.data.frame(test1@listData)
row.names(test1result) <- test1@rownames

colnames(test1result)
#colnames(test1result) <- c("R18Mean","R18log2FoldChange","lfcSE","stat","R18RawPValue", "R18adjPvalue", "R18FoldChange")
colnames(test1result) <- c("R10Mean","R10log2FoldChange","lfcSE","stat","R10RawPValue", "R10adjPvalue", "R10FoldChange")


colnames(ASresult)
colnames(ASresult) <- c("ASlog2FoldChange", "ASFoldChange","ASRawPvalue", "ASadjPvalue","ASGroupMax")


#combine two tables
newdata <- merge(test1result,ASresult, by="row.names")



colnames(newdata)

#another way to convert to numeric
cols = c(9, 10, 11, 12)
newdata[,cols] = apply(newdata[,cols], 2, function(x) as.numeric(as.character(x)))


#remove rows with NA value
row.has.na <- apply(newdata, 1, function(x){any(is.na(x))})
sum(row.has.na)
newdata <- newdata[!row.has.na,]

str(newdata)



myscatter <- function(x, y, data, title){
  plot(x, y, main=title,
       xla="DESeq10", ylab="AS", pch=19)
  abline(fit <- lm(y ~ x, data=data), col='red')
  my.p = summary(fit)$coefficients[2,4]
  my.R2 = summary(fit)$adj.r.squared
  
  rp = vector('expression',2)
  rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                     list(MYVALUE = format(my.R2,dig=3)))[2]
  rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                     list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
  legend('topleft', legend = rp, bty = 'n')
}

attach(newdata)
par(mfrow=c(2,2))
myscatter(R10FoldChange,ASFoldChange, newdata, "FoldChange" )
myscatter(R10log2FoldChange,ASlog2FoldChange, newdata, "log2FoldChang")
myscatter(R10RawPValue,ASRawPvalue, newdata, "RawPvalue" )
myscatter(R10adjPvalue,ASadjPvalue, newdata,"adjPvalue" )

write.table(newdata, "R10_AS.txt", sep="\t", col.names = TRUE, quote=FALSE)

#one way to convert factor to numeric
for(i in c(1,ncol(newdata))) {
  newdata[,i] <- as.numeric(as.character(newdata[,i]))
}



#--------------load R18_AS.txt and compare with each other--------------------#

R10 <- test1result[,c(1,2,5,6,7)]
R18_AS <- read.table("R18_AS.txt", sep="\t", header=T)
row.names(R18_AS) <- R18_AS[,1]

#combine two tables
alldata <- merge(R10,R18_AS, by="row.names")

attach(alldata)
par(mfrow=c(2,2))
myscatter2(R18FoldChange,R10FoldChange, alldata, "FoldChange" )
myscatter2(R18log2FoldChange,R10log2FoldChange, alldata, "log2FoldChang")
myscatter2(R18RawPValue,R10RawPValue, alldata, "RawPvalue" )
myscatter2(R18adjPvalue,R10adjPvalue, alldata,"adjPvalue" )


myscatter2 <- function(x, y, data, title){
  
  plot(x, y, main=title,
       xla="DESeq18", ylab="DESeq10", pch=19)
  abline(fit <- lm(y ~ x, data=data), col='red')
  my.p = summary(fit)$coefficients[2,4]
  my.R2 = summary(fit)$adj.r.squared
  
  rp = vector('expression',2)
  rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                     list(MYVALUE = format(my.R2,dig=3)))[2]
  rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                     list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
  legend('topleft', legend = rp, bty = 'n')
}

















df <- data.frame(x = c(1:100))
df$y <- 2 + 3 * df$x + rnorm(100, sd = 40)

df <- newdata[,c("R18FoldChange", "ASFoldChange")]
df <- newdata[,c("R18RawPValue", "ASRawPvalue")]
colnames(df) <-c("x", "y")


lm_eqn = function(m) {
  
  l <- list(a = format(coef(m)[1], digits = 2),
            b = format(abs(coef(m)[2]), digits = 2),
            r2 = format(summary(m)$r.squared, digits = 3));
  
  if (coef(m)[2] >= 0)  {
    eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2,l)
  } else {
    eq <- substitute(italic(y) == a - b %.% italic(x)*","~~italic(r)^2~"="~r2,l)    
  }
  
  as.character(as.expression(eq));                 
}

p <- ggplot(data = df, aes(x = x, y = y)) +
  geom_smooth(method = "lm", se=FALSE, color="black", formula = y ~ x) +
  geom_point()

p1 = p + geom_text(aes(x = 25, y = 300, label = lm_eqn(lm(y ~ x, df))), parse = TRUE)
p1


df.labs <- data.frame(x = 25, y = 300, label = lm_eqn(newdata))
p <- p + geom_text(data = df.labs, aes(x = x, y = y, label = label), parse = TRUE)
