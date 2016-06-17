## Load the packages
source("https://bioconductor.org/biocLite.R")
biocLite("edgeR")
install.packages("statmod")
install.packages("locfit")
install.packages("xlsx")
install.packages("rJava")

library(rJava)
library(edgeR)
library(statmod)
library(locfit)


## Set the path
setwd("./")
getwd()
rm(list=ls())


## Load the data
rawdata <- read.csv("./data/pnas.1402835111.sd01.csv", head=TRUE, sep=",")
head(rawdata)
cds <- DGEList(counts=rawdata[,2:64], genes=rawdata[,1])


## Estimate library size and do normilization
cds$samples$lib.size <- colSums(cds$counts)
cds <- calcNormFactors(cds, method="TMM")
cds$samples


## Extract information for each species
lyc <- cds[,1:21]
pen <- cds[,22:42]
hab <- cds[,43:63]


## Filtering out lowly expressed genes
(nrow(lyc)+nrow(pen)+nrow(hab))/20
keep <- rowSums(lyc$counts)/sum(lyc$samples$lib.size)*1e6 >= 5
lyc <- lyc[keep,]

keep <- rowSums(pen$counts)/sum(pen$samples$lib.size)*1e6 >= 5
pen <- pen[keep,]

keep <- rowSums(hab$counts)/sum(hab$samples$lib.size)*1e6 >= 5
hab <- hab[keep,]
(nrow(lyc)+nrow(pen)+nrow(hab))/20


## Extract normilized counts for following tests
lyc_genes <- lyc$genes
lyc_nd <- cpm(lyc, normalized.lib.sizes=TRUE)

pen_genes <- pen$genes
pen_nd <- cpm(pen, normalized.lib.sizes=TRUE)

hab_genes <- hab$genes
hab_nd <- cpm(hab, normalized.lib.sizes=TRUE)


## Calculate the coefficient of variation for each species
mean <- apply(lyc_nd, 1, mean)
sd <- apply(lyc_nd, 1, sd)
cof <- sd/mean
merged = rbind(sd, mean, cof)
merged <- t(merged)
merged <- cbind(lyc_genes, merged, lyc_nd)
head(merged)
write.csv(merged, "./results/lyc_dev.csv") 


mean <- apply(pen_nd, 1, mean)
sd <- apply(pen_nd, 1, sd)
cof <- sd/mean
merged = rbind(sd, mean, cof)
merged <- t(merged)
merged <- cbind(pen_genes, merged, pen_nd)
head(merged)
write.csv(merged, "./results/pen_dev.csv") 


mean <- apply(hab_nd, 1, mean)
sd <- apply(hab_nd, 1, sd)
cof <- sd/mean
merged = rbind(sd, mean, cof)
merged <- t(merged)
merged <- cbind(hab_genes, merged, hab_nd)
head(merged)
write.csv(merged, "./results/hab_dev.csv") 


## read the newly formed three datasets
lycF <- read.csv("./results/lyc_dev.csv", head=TRUE, sep=",", colClasses=c("NULL",rep(NA,25)))
penF <- read.csv("./results/pen_dev.csv", head=TRUE, sep=",", colClasses=c("NULL",rep(NA,25)))
habF <- read.csv("./results/hab_dev.csv", head=TRUE, sep=",", colClasses=c("NULL",rep(NA,25)))


## change to specific gene names
lycF$genes <- gsub('Solyc', 'lyc_Solyc', lycF$genes)
penF$genes <- gsub('Solyc', 'pen_Solyc', penF$genes)
habF$genes <- gsub('Solyc', 'hab_Solyc', habF$genes)


## change the column names
colnames(lycF) <- c("genes", "sd", "mean", "cof", "P3_m_1", "P3_m_2", "P3_m_3", "P4_d_1", "P4_d_2,", "P4_d_3", "P4_p_1", 
                    "P4_p_2", "P4_p_3", "P5_d_1", "P5_d_2", "P5_d_3", "P5_p_1", "P5_p_2", "P5_p_3", 
                    "P6_d_1", "P6_d_2", "P6_d_3", "P6_p_1", "P6_p_2", "P6_p_3")
colnames(penF) <- c("genes", "sd", "mean", "cof", "P3_m_1", "P3_m_2", "P3_m_3", "P4_d_1", "P4_d_2,", "P4_d_3", "P4_p_1", 
                    "P4_p_2", "P4_p_3", "P5_d_1", "P5_d_2", "P5_d_3", "P5_p_1", "P5_p_2", "P5_p_3", 
                    "P6_d_1", "P6_d_2", "P6_d_3", "P6_p_1", "P6_p_2", "P6_p_3")
colnames(habF) <- c("genes", "sd", "mean", "cof", "P3_m_1", "P3_m_2", "P3_m_3", "P4_d_1", "P4_d_2,", "P4_d_3", "P4_p_1", 
                    "P4_p_2", "P4_p_3", "P5_d_1", "P5_d_2", "P5_d_3", "P5_p_1", "P5_p_2", "P5_p_3", 
                    "P6_d_1", "P6_d_2", "P6_d_3", "P6_p_1", "P6_p_2", "P6_p_3")
head(lycF)


## merge three datasets together and sorted by BCV
merged = rbind(lycF, penF, habF)
head(merged)
tail(merged)
sorted <- merged[order(merged$cof, decreasing=TRUE), ]
write.csv(sorted, "./results/all_genes.csv") 


## print out the target most variable genes
cutoff <- nrow(sorted)/20
quantile = head(sorted, cutoff)
head(quantile)


## average the replicate values
quantile2 <- within(quantile2, P3_m <- apply(quantile[5:7], 1, mean))
quantile2 <- within(quantile2, P4_p <- apply(quantile[11:13], 1, mean))
quantile2 <- within(quantile2, P5_p <- apply(quantile[17:19], 1, mean))
quantile2 <- within(quantile2, P6_p <- apply(quantile[23:25], 1, mean))
quantile2 <- within(quantile2, P4_d <- apply(quantile[8:10], 1, mean))
quantile2 <- within(quantile2, P5_d <- apply(quantile[14:16], 1, mean))
quantile2 <- within(quantile2, P6_d <- apply(quantile[20:22], 1, mean))
head(quantile2)
quantile2 <- quantile2[,-c(5:25)] 


## write the final results
write.csv(quantile, "./results/DEX_genes.csv") 
write.csv(quantile2, "./data/DEX_genes.csv") 

