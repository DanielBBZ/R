setwd("//192.168.1.174/ws03/IData/Users/Daniel/Support/TestCase/HierarchicalCluster_R")

### Load data matrix ###
datamatrix <- read.table("Table.txt", sep="\t", header=T, row.names=1)
d <- dist(datamatrix)
hc <- hclust(d)

plot(hc)


datamatrixT <- read.table("TransposeTable.txt", sep="\t", header=T, row.names=1)
dT <- dist(datamatrixT, method = "manhattan")
hcT <- hclust(dT)
plot(hcT)
