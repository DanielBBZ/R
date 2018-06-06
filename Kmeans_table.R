### Load data matrix ###
dat <- read.table("TsneScores.txt", sep="\t", header=T, row.names=1)


Kmeanslow <- 3
Kmeanshigh <- 3

set.seed(10)

i <- Kmeanslow

while (i <= Kmeanshigh) {
  clusters <- kmeans(dat[,1:2], i)
  
  newclass <- paste("KMeans",i,sep = "_")
  # Save the cluster number in the dataset as column 'class'
  
  dat$class <- as.factor(clusters$cluster)
  
  colnames(dat)[which(names(dat) == "class")] <- newclass
  i <- i + 1}

TSNE_Kmeans <- dat





attach(mydata)
plot(V1, V2, main="T-sNE color plot",
     xla="R", ylab="AS", pch=19)



library(ggplot2)

qplot(V1, V2, colour = class, data = mydata)


colnames(mydata)[which(names(mydata) == "class")] <- "newColumnName"
