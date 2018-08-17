
#Purpose: compare two tables, create scatter plot and show R2 and P value:


library("ggplot2")

setwd("E:/R/comparison/")

### Load data matrix ###
AS1ng <- read.table("AS1ng.txt", sep="\t", header=T, row.names=1)
GG1ng <- read.table("GG1ng.txt", sep="\t", header=T, row.names=1)

colnames(AS1ng)
colnames(GG1ng)

#combine two tables
newdata <- merge(AS1ng,GG1ng, by="row.names")
colnames(newdata)

#-----------------------------optional part-----------------------------
#one way to convert factor to numeric
for(i in c(1,ncol(newdata))) {
  newdata[,i] <- as.numeric(as.character(newdata[,i]))
}

#another way to convert to numeric
cols = c(9, 10, 11, 12)
newdata[,cols] = apply(newdata[,cols], 2, function(x) as.numeric(as.character(x)))


#remove rows with NA value
row.has.na <- apply(newdata, 1, function(x){any(is.na(x))})
sum(row.has.na)
newdata <- newdata[!row.has.na,]
#-----------------------------optional part-----------------------------


str(newdata)


myscatter <- function(x, y, data, title){
  plot(x, y, main=title,
       xla="AS", ylab="GG", pch=19)
  abline(fit <- lm(y ~ x, data=data), col='red')
  my.p = summary(fit)$coefficients[2,4]
  my.R2 = summary(fit)$adj.r.squared
  my.Coefficients = summary(fit)$coefficients[2,1]
  
  rp = vector('expression',3)
  rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                     list(MYVALUE = format(my.R2,dig=3)))[2]
  rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                     list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
  rp[3] = substitute(expression(italic(Cf) == MYTHIRDVALUE), 
                     list(MYTHIRDVALUE = format(my.Coefficients, digits = 3)))[2]
  legend('topleft', legend = rp, bty = 'n')
}

attach(newdata)
par(mfrow=c(2,2))
myscatter(AS_1ngS1_1,GG_1ngS1_1, newdata, "1ngS1_1")
myscatter(AS_1ngS1_2,GG_1ngS1_2, newdata, "1ngS1_2")
myscatter(AS_1ngS1_3,GG_1ngS1_3, newdata, "1ngS1_3")
myscatter(AS_1ngS1_4,GG_1ngS1_4, newdata, "1ngS1_4")



#write.table(newdata, "combined.txt", sep="\t", col.names = TRUE, quote=FALSE)





#--------------load another table and compare--------------------#
### Load data matrix ###
AlignReport <- read.table("AlignReport_GGAS.txt", sep="\t", header=T, row.names=1)

str(AlignReport)

colnames(AlignReport)

attach(AlignReport)
par(mfrow=c(1,2))
myscatter(AS.Total_read,GG.Total_read, AlignReport, "Total_read No.")
myscatter(AS.Uniquely_mapped_read,GG.Uniquely_mapped_read, AlignReport, "Uniquely_mapped_read No.")
myscatter(AS.Non.uniquely_mapped_read,GG.Non.uniquely_mapped_read, AlignReport, "Non.uniquely_mapped_read No.")
myscatter(AS.Unmapped_read,GG.Unmapped_read, AlignReport, "Unmapped_read No.")




#--------------log2 transform the data--------------------#

AS1nglog2<-log2(AS1ng+0.01)
GG1nglog2<-log2(GG1ng+0.01)
#combine two tables
newdatalog2 <- merge(AS1nglog2,GG1nglog2, by="row.names")
colnames(newdatalog2)


attach(newdatalog2)
par(mfrow=c(2,2))
myscatter(AS_1ngS1_1,GG_1ngS1_1, newdatalog2, "1ngS1_1")
myscatter(AS_1ngS1_2,GG_1ngS1_2, newdatalog2, "1ngS1_2")
myscatter(AS_1ngS1_3,GG_1ngS1_3, newdatalog2, "1ngS1_3")
myscatter(AS_1ngS1_4,GG_1ngS1_4, newdatalog2, "1ngS1_4")

