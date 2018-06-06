# Take a table matrix as input, output the 95% confidence of interval for the median, 
#applicable both to row or column

# ci.median function according to the R package asbio
ci.median<-function(x,conf=.95){
  n<-nrow(as.matrix(x))
  if(qbinom((1-conf)/2,n,0.5)==0)stop("CI not calculable")
  L<- qbinom((1-conf)/2,n,0.5)
  U<-n-L+1
  if(L>=U)stop("CI not calculable")
  order.x<-sort(x)
  #res<-list()
  #res$head<-paste(paste(as.character(conf*100),"%",sep=""),c("Confidence interval for population median"))
  #res$ci<-c(median=median(x),lower=order.x[L],upper=order.x[n-L+1])
  #res$ends<-c("Estimate",paste(as.character(c((1-conf)/2,1-((1-conf)/2))*100),"%",sep=""))
  #res$coverage<-1-(2*pbinom(q=L-1,n,0.5))
  #class(res)<-"ci"
  #res
  res <-c()
  res <-c(median=median(x),lower=order.x[L],upper=order.x[n-L+1])
  res
}


#setwd("//192.168.1.174/ws03/IData/Users/Daniel/Support/TestCase/R/Benchmark")

### Load data matrix ###
#datamatrix <- read.table("MicroArrayData.txt", sep="\t", header=T, row.names=1)

#x<-rnorm(20)

byrow <- "FALSE"


if (as.logical(byrow)){
  myci <- t(apply(datamatrix, 1, ci.median))
} else {
  myci <- t(apply(datamatrix, 2, ci.median))
}

colnames(myci) <- c("median", "lower", "higher")

head(myci)
