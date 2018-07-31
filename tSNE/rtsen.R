# load the Rtsne package

library(Rtsne)
library(ggplot2)

# loading data into R
mydata <- read.table("Z:\\Users\\Daniel\\Test\\SCData\\10X_GSE99457\\ProxGFP.txt", header=TRUE, sep="\t", row.names=1)

mydata <- read.table("Z:\\Users\\Daniel\\Test\\SCData\\GSE81904\\CountData\\AScount_SubSub.txt", header=TRUE, sep="\t", row.names=1)
mydata <- read.table("Z:\\Users\\Daniel\\Test\\SCData\\GSE84498\\Landcount\\LandCountGene.txt", header=TRUE, sep="\t", row.names=1)
mydata <- read.table("Z:\\Users\\Daniel\\Test\\SCData\\GSE84498\\Landcount\\GEOCountnew.txt", header=TRUE, sep="\t", row.names=1)

genename <- row.names(mydata)
samplename <- colnames(mydata)

#Preprocess data:
newcount <- as.matrix(mydata)
row.names(newcount) <- NULL #remove row name
trans_newcount <- t(newcount) #transpose matrix, so each row will be a sample
row.names(trans_newcount) <- NULL  #remove row name

#trans_newcount <- trans_newcount[-1,] #remove first row of the dataframe

#trans_newcount <- as.numeric(trans_newcount)
set.seed(2)
# run Rtsne with default parameters
rtsne_out <- Rtsne(trans_newcount)
rtsne_out2 <- Rtsne(trans_newcount)
rtsne_out3 <- Rtsne(trans_newcount)


## getting the two dimension matrix
d_tsne_1 = as.data.frame(rtsne_out$Y)
d_tsne_1 = as.data.frame(rtsne_out2$Y) 
d_tsne_1 = as.data.frame(rtsne_out3$Y) 

## plotting the results without clustering
ggplot(d_tsne_1, aes(x=V1, y=V2)) +  
  geom_point(size=0.25) +
  guides(colour=guide_legend(override.aes=list(size=6))) +
  xlab("") + ylab("") +
  ggtitle("t-SNE") +
  theme_light(base_size=20) +
  theme(axis.text.x=element_blank(),
        axis.text.y=element_blank()) +
  scale_colour_brewer(palette = "Set2")




## keeping original data
d_tsne_1_original=d_tsne_1

## Creating k-means clustering model, and assigning the result to the data used to create the tsne
fit_cluster_kmeans=kmeans(scale(d_tsne_1), 6)  
d_tsne_1_original$cl_kmeans = factor(fit_cluster_kmeans$cluster)

## Creating hierarchical cluster model, and assigning the result to the data used to create the tsne
fit_cluster_hierarchical=hclust(dist(scale(d_tsne_1)))

## setting 3 clusters as output
d_tsne_1_original$cl_hierarchical = factor(cutree(fit_cluster_hierarchical, k=6)) 



#Plotting the cluster models onto t-SNE output
#Now time to plot the result of each cluster model, based on the t-SNE map.

plot_cluster=function(data, var_cluster, palette)  
{
  ggplot(data, aes_string(x="V1", y="V2", color=var_cluster)) +
    geom_point(size=0.25) +
    guides(colour=guide_legend(override.aes=list(size=2))) +
    xlab("") + ylab("") +
    ggtitle("") +
    theme_light(base_size=20) +
    theme(axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          legend.direction = "horizontal", 
          legend.position = "bottom",
          legend.box = "horizontal") + 
    scale_colour_brewer(palette = palette) 
}


plot_k=plot_cluster(d_tsne_1_original, "cl_kmeans", "Accent")  
plot_h=plot_cluster(d_tsne_1_original, "cl_hierarchical", "Set1")

## and finally: putting the plots side by side with gridExtra lib...
library(gridExtra)  
grid.arrange(plot_k, plot_h,  ncol=2)  



row.names(d_tsne_1) <- samplename
write.table(d_tsne_1, "Z:\\Users\\Daniel\\Test\\SCData\\GSE84498\\Landcount\\tSNE_GEOCount_sub2.txt", sep="\t", row.names=TRUE, col.names = TRUE, quote=FALSE)

write.table(d_tsne_1, "Z:\\Users\\Daniel\\Test\\SCData\\GSE81904\\CountData\\tSNE_SubSub.txt", sep="\t", row.names=TRUE, col.names = TRUE, quote=FALSE)


# plot the output of Rtsne into d:\\barneshutplot.jpg file of 2400x1800 dimension
jpeg("Z:\\Users\\Daniel\\Support\\scripts\\R_code\\RtSNE\\barneshutplot2.jpg", width=2400, height=1800)
plot(rtsne_out2$Y, t='n', main="BarnesHutSNE")
text(rtsne_out2$Y, labels=colnames(mydata))