setwd("//192.168.1.174/ws03/IData/Users/Daniel/Test/TestProject/PCA/")

dat <- read.table("data.txt", sep="\t", header=T, row.names=1, check.names = FALSE)
#des <- read.table("design.txt", sep="\t", header=T, row.names=1, check.names = FALSE)

tsne <- read.table("NewTsneScores.txt", sep="\t", header=T, row.names=1, check.names = FALSE)


library(ggplot2)

plot1 = ggplot(tsne, aes(V1, V2)) +
  geom_point(aes(colour = t(dat)[,"Probe1"])) +
  ggtitle("Probe1")

plot2 = ggplot(tsne, aes(V1, V2)) +
  geom_point(aes(colour = t(dat)[,"Probe2"])) +
  ggtitle("Probe2")


plot3 = ggplot(tsne, aes(V1, V2)) +
  geom_point(aes(colour = t(dat)[,"Probe3"])) +
  ggtitle("Probe3")

plot4 = ggplot(tsne, aes(V1, V2)) +
  geom_point(aes(colour = t(dat)[,"Probe4"])) +
  ggtitle("Probe4")

multiplot(plot1, plot2, plot3, plot4, cols = 2)


install.packages("rlist")
library(rlist)

genelist <- c("Probe1", "Probe2", "Probe3")
plot_list <- list()
colorscale <- list()
for (value in genelist){
  p = ggplot(PCtable, aes(PC1, PC2)) +
               geom_point(aes(colour = t(dat)[,value])) +
               ggtitle(value)
  list.append(plot_list, p)
}
  
#grid.arrange(plot_list[[1]], plot_list[[2]], plot_list[[3]], ncol=2, nrow = 2)

multiplot(plot_list, cols=2)

# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
