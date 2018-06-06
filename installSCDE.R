
require(devtools)
install_version("flexmix", version = "2.3-13", repos = "http://cran.us.r-project.org")


# install SCDE package from Github

install.packages("devtools")
library("devtools")
install_github('hms-dbmi/scde')

library("scde")

#data(es.mef.small)
#cd2 <- clean.counts(es.mef.small, min.lib.size=1000, min.reads = 1, min.detected = 1)

setwd("//192.168.1.174/ws03/IData/Users/Daniel/Support/scripts/R_code/RSCDE/")

Count <- read.table("Z:\\Users\\Daniel\\Support\\scripts\\R_code\\RSCDE\\Count.txt", header=TRUE, sep="\t", row.names=1)
Design <- read.table("Z:\\Users\\Daniel\\Support\\scripts\\R_code\\RSCDE\\Design.txt", header=TRUE, sep="\t", row.names=1)


setwd("//192.168.1.174/ws03/IData/Users/Daniel/Test/SCData/10X_GSE99457/FullProcess/SCDE/")
Count <- read.table("Prox1-GFP.SingleCellCounts.txt", header=TRUE, sep="\t", row.names=1)
Design <- read.table("Design.txt", header=TRUE, sep="\t", row.names=1)


Design$DiseaseState = factor(Design$DiseaseState, levels=c("Control", "Type1D_09", "Type1D_22", 
                                                           "Type1D_40", "Type1D_50", "Type1D_77", "Type1D_79"))
Design$Disease = factor(Design$Disease, level = c('Control', 'Type1D'))


#table(Design$Disease)
#table(Design$DiseaseState)

# clean up the dataset
Count<-apply(Count,2,function(x) {storage.mode(x) <- 'integer'; x})
Count.df <- as.data.frame(Count)
cd <- clean.counts(Count.df, min.lib.size=1000, min.reads = 1, min.detected = 1)


# calculate models
o.ifm <- scde.error.models(counts = cd, n.cores = 1, 
                           threshold.segmentation = TRUE, save.crossfit.plots = FALSE, 
                           save.model.plots = FALSE, verbose = 1)





library(scde)
data(pollen)
cd <- clean.counts(pollen)
## change verbose to 1 to get immediate feedback
knn <- knn.error.models(cd, k = ncol(cd)/4, n.cores = 1, min.count.threshold = 2, min.nonfailed = 5, max.model.plots = 10, verbose=1)
