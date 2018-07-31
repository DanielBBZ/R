

setwd("//input/PCA/")

dat <- read.table("data.txt", sep="\t", header=T, row.names=1, check.names = FALSE)
des <- read.table("design.txt", sep="\t", header=T, row.names=1, check.names = FALSE)


#principal component analysis
prin_comp <- prcomp(t(dat), scale = TRUE)


## make a scree plot
#pca.var <- prin_comp$sdev^2
#pca.var.per <- round(pca.var/sum(pca.var)*100, 1)
#barplot(pca.var.per, main="Scree Plot", xlab="Principal Component", ylab="Percent Variation")




#output PC loading table:
LoadingTable <- cbind(prin_comp$sdev, (prin_comp$sdev)^2, (prin_comp$sdev)^2/sum((prin_comp$sdev)^2))
colnames(LoadingTable) <- c("stdev", "variance", "proportion")
rownames(LoadingTable) <- paste("PC", 1:nrow(LoadingTable), sep="")

#output PCtable with first three principle components
PCtable <- prin_comp$x[,1:3]
PCtable <- cbind(PCtable, des)
colnames(PCtable)[1] <- paste("PC1 ", round(LoadingTable[1,3]*100, digits = 2), "%", sep="")
colnames(PCtable)[2] <- paste("PC2 ", round(LoadingTable[2,3]*100, digits = 2), "%", sep="")
colnames(PCtable)[3] <- paste("PC3 ", round(LoadingTable[3,3]*100, digits = 2), "%", sep="")

  
## get the loading scores for each gene to the first three pc.
VariableLoading <- prin_comp$rotation[,1:3]
colnames(VariableLoading)[1:3] <- colnames(PCtable)[1:3]


gr <- factor(des[,"Group1"])
gr2 <- factor(des[,"Group2"])
pca3d(prin_comp,  components = 1:3, group=gr2, show.group.labels = TRUE)



#2D plot
## now make a fancy looking plot that shows the PCs and the variation:
library(ggplot2)
pca.data <- data.frame(Sample=rownames(pca$x),
                       X=pca$x[,1],
                       Y=pca$x[,2])
pca.data
ggplot(data=pca.data, aes(x=X, y=Y, label=Sample)) +
  geom_text() +
  xlab(paste("PC1 - ", pca.var.per[1], "%", sep="")) +
  ylab(paste("PC2 - ", pca.var.per[2], "%", sep="")) +
  theme_bw() +
  ggtitle("My PCA Graph")

## get the name of the top 10 measurements (genes) that contribute
## most to pc1.
loading_scores <- pca$rotation[,1]
gene_scores <- abs(loading_scores) ## get the magnitudes
gene_score_ranked <- sort(gene_scores, decreasing=TRUE)
top_10_genes <- names(gene_score_ranked[1:10])

top_10_genes ## show the names of the top 10 genes

pca$rotation[top_10_genes,1] ## show the scores (and +/- sign)



#3D plot
library(pca3d)
data(metabo)

pca <- prcomp(metabo[,-1], scale. = TRUE)
gr <- factor(metabo[,1])
summary(gr)

pca3d(pca, group=gr, show.centroids=TRUE, show.group.labels = TRUE)
pca3d(pca, group=gr, fancy=TRUE)

snapshotPCA3d(file=first_plot.png)


pca3d( pca, group= metabo[,1],
       show.centroids=TRUE, bg= "white",
       axes.color= "black", new= TRUE )
