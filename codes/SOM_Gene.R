## Load the packages
install.packages("caret")
install.packages("som") 
install.packages("devtools")
install_github("ggbiplot2", "vqv")

library(som) 
require(caret)
library(devtools)
library(ggbiplot)


## Set the path
setwd("../")
getwd()
rm(list=ls())


## Load the data
jalt <- read.csv("./data/DEX_genes.csv", head=TRUE, sep=",")
head(jalt)
#jalt <- jalt[,6:12]
jalt <- jalt[,6:12]
nrow(jalt)


## Cluster genes using SOM algorithm
jalt.n <- normalize(jalt)
head(jalt.n)
nrow(jalt.n)
#foo <- som(jalt.n, xdim=1, ydim=3, topol="hexa", rlen=100, alpha=c(0.008,0.007))
foo <- som(jalt.n, xdim=1, ydim=3, topol="hexa", rlen=100)
plot(foo)
jalt.cluster <- data.frame(foo$visual$y)


## Change the cluster names
colnames(jalt.cluster) <- c("cluster")
jalt.cluster$cluster[jalt.cluster$cluster == "0"] <- "Cluster 3"
jalt.cluster$cluster[jalt.cluster$cluster == "1"] <- "Cluster 1"
jalt.cluster$cluster[jalt.cluster$cluster == "2"] <- "Cluster 2"
head(jalt.cluster,50)


## Count the number of genes in each cluster
length(which(jalt.cluster == "Cluster 1"))
length(which(jalt.cluster == "Cluster 2"))
length(which(jalt.cluster == "Cluster 3"))


## Perform PCA
## Load the data
jalt.n <- data.frame(jalt.n)
colnames(jalt.n) <- c("P3_m", "P4_p", "P5_p", "P6_p", "P4_d", "P5_d", "P6_d")
jalt.cluster <- jalt.cluster[, 1]
head(jalt.n)


# log transform 
#log.jalt <- log10(jalt+1)


# apply PCA - scale. = TRUE is highly 
# advisable, but default is FALSE. 
jalt.pca <- prcomp(jalt.n, scale = TRUE)


# print method
print(jalt.pca)
plot(jalt.pca, type = "l")
summary(jalt.pca)

  
# Plot PCA with clusters markers from SOM
g <- ggbiplot(jalt.pca,obs.scale = 1, 
         var.scale=1,groups=jalt.cluster)
g <- g + theme(panel.background = element_blank())
g <- g + theme(axis.line.x = element_line(color="black", size = 0.5),
            axis.line.y = element_line(color="black", size = 0.5))
plot(g)         


## Save the result image
ggsave("./results/myPCA.pdf")