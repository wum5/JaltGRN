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
group <- factor(c(rep("lyc", 21), rep("pen",21), rep("hab",21)))
cds <- DGEList(counts=rawdata[,2:64], genes=rawdata[,1], group = group)


## Order transcripts by counts number
o <- order(rowSums(cds$counts), decreasing=TRUE)
cds <- cds[o,]
nrow(cds)


## Extract information for each species
lyc <- cds[,1:21]
pen <- cds[,22:42]
hab <- cds[,43:63]


## Filtering out lowly expressed genes
keep <- rowSums(cpm(lyc)>1) >= ncol(lyc)/2
lyc <- lyc[keep,]
nrow(lyc)

keep <- rowSums(cpm(pen)>1) >= ncol(pen)/2
pen <- pen[keep,]
nrow(pen)

keep <- rowSums(cpm(hab)>1) >= ncol(hab)/2
hab <- hab[keep,]
nrow(hab)


## Estimate library size and do normilization
lyc$samples$lib.size <- colSums(lyc$counts)
lyc <- calcNormFactors(lyc, method="TMM")
lyc$samples

pen$samples$lib.size <- colSums(pen$counts)
pen <- calcNormFactors(pen, method="TMM")
pen$samples

hab$samples$lib.size <- colSums(hab$counts)
hab <- calcNormFactors(hab, method="TMM")
hab$samples


## Calculate dispersion for lyco
lyc <- estimateDisp(lyc)
lyc$common.dispersion
lyc_bcv <- cbind(lyc$genes,lyc$tagwise.dispersion,lyc$counts)


## Calculate dispersion for penn
pen <- estimateDisp(pen)
pen$common.dispersion
pen_bcv <- cbind(pen$genes,pen$tagwise.dispersion,pen$counts)


## Calculate dispersion for harb
hab <- estimateDisp(hab)
hab$common.dispersion
hab_bcv <- cbind(hab$genes,hab$tagwise.dispersion,hab$counts)


## Write down the output files
write.csv(lyc_bcv, "./data/lyc_bcv.csv") 
write.csv(pen_bcv, "./data/pen_bcv.csv")
write.csv(hab_bcv, "./data/hab_bcv.csv")


## read the newly formed three datasets
lycF <- read.csv("./data/lyc_bcv.csv", head=TRUE, sep=",", colClasses=c("NULL",rep(NA,23)))
penF <- read.csv("./data/pen_bcv.csv", head=TRUE, sep=",", colClasses=c("NULL",rep(NA,23)))
habF <- read.csv("./data/hab_bcv.csv", head=TRUE, sep=",", colClasses=c("NULL",rep(NA,23)))


## change to specific gene names
lycF$genes <- gsub('Solyc', 'lyc_Solyc', lycF$genes)
penF$genes <- gsub('Solyc', 'pen_Solyc', penF$genes)
habF$genes <- gsub('Solyc', 'hab_Solyc', habF$genes)


## change the column names
colnames(lycF) <- c("genes", "bcv", "P3_m_1", "P3_m_2", "P3_m_3", "P4_d_1", "P4_d_2,", "P4_d_3", "P4_p_1", 
                    "P4_p_2", "P4_p_3", "P5_d_1", "P5_d_2", "P5_d_3", "P5_p_1", "P5_p_2", "P5_p_3", 
                    "P6_d_1", "P6_d_2", "P6_d_3", "P6_p_1", "P6_p_2", "P6_p_3")
colnames(penF) <- c("genes", "bcv", "P3_m_1", "P3_m_2", "P3_m_3", "P4_d_1", "P4_d_2,", "P4_d_3", "P4_p_1", 
                    "P4_p_2", "P4_p_3", "P5_d_1", "P5_d_2", "P5_d_3", "P5_p_1", "P5_p_2", "P5_p_3", 
                    "P6_d_1", "P6_d_2", "P6_d_3", "P6_p_1", "P6_p_2", "P6_p_3")
colnames(habF) <- c("genes", "bcv", "P3_m_1", "P3_m_2", "P3_m_3", "P4_d_1", "P4_d_2,", "P4_d_3", "P4_p_1", 
                    "P4_p_2", "P4_p_3", "P5_d_1", "P5_d_2", "P5_d_3", "P5_p_1", "P5_p_2", "P5_p_3", 
                    "P6_d_1", "P6_d_2", "P6_d_3", "P6_p_1", "P6_p_2", "P6_p_3")
colnames(lycF)


## merge three datasets together and sorted by BCV
merged = rbind(lycF, penF, habF)
head(merged)
tail(merged)
sorted <- merged[order(merged$bcv, decreasing=TRUE), ]
write.csv(sorted, "./data/all_genes.csv") 


## print out the target most variable genes
cutoff <- nrow(sorted)/20
quantile = head(sorted, cutoff)
write.csv(quantile, "./data/DEX_genes.csv") 

