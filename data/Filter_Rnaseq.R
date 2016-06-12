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


## Load the data
rawdata <- read.csv("pnas.1402835111.sd01.csv", head=TRUE, sep=",")
head(rawdata)
group <- factor(c(rep("lyc", 21), rep("pen",21), rep("hab",21)))
cds <- DGEList(counts=rawdata[,2:64], genes=rawdata[,1], group = group)


## Order transcripts by counts number
o <- order(rowSums(cds$counts), decreasing=TRUE)
cds <- cds[o,]
nrow(cds)


## Filtering out lowly expressed genes
keep <- rowSums(cpm(cds)>1) > ncol(cds)/2
cds <- cds[keep,]
nrow(cds)


## Estimate library size and do normilization
cds$samples$lib.size <- colSums(cds$counts)
cds <- calcNormFactors(cds, method="TMM")
cds$samples


## Extract information for each species
lyc <- cds[,1:21]
pen <- cds[,22:42]
hab <- cds[,43:63]


## Calculate dispersion for lyco
lyc <- estimateDisp(lyc)
lyc$common.dispersion
lyc_bcv <- cbind(lyc$genes,lyc$tagwise.dispersion)


## Calculate dispersion for penn
pen <- estimateDisp(pen)
pen$common.dispersion
pen_bcv <- cbind(pen$genes,pen$tagwise.dispersion)


## Calculate dispersion for harb
hab <- estimateDisp(hab)
hab$common.dispersion
hab_bcv <- cbind(hab$genes,hab$tagwise.dispersion)


## Write down the output files
write.csv(lyc_bcv, "./lyc_bcv.csv") 
write.csv(pen_bcv, "./pen_bcv.csv")
write.csv(hab_bcv, "./hab_bcv.csv")


## read the newly formed three datasets
lycF <- read.csv("lyc_bcv.csv", head=TRUE, sep=",", colClasses=c("NULL",NA,NA))
penF <- read.csv("pen_bcv.csv", head=TRUE, sep=",", colClasses=c("NULL",NA,NA))
habF <- read.csv("hab_bcv.csv", head=TRUE, sep=",", colClasses=c("NULL",NA,NA))


## change to specific gene names
lycF$genes <- gsub('Solyc', 'lyc_Solyc', lycF$genes)
penF$genes <- gsub('Solyc', 'pen_Solyc', penF$genes)
habF$genes <- gsub('Solyc', 'hab_Solyc', habF$genes)


## change the column names
colnames(lycF) <- c("genes", "bcv")
colnames(penF) <- c("genes", "bcv")
colnames(habF) <- c("genes", "bcv")
colnames(lycF)


## merge three datasets together and sorted by BCV
merged = rbind(lycF, penF, habF)
head(merged)
tail(merged)
sorted <- merged[order(merged$bcv, decreasing=TRUE), ]
write.csv(quantile, "./all_genes.csv") 


## print out the target most variable genes
cutoff <- nrow(sorted)/20
quantile = head(sorted, cutoff)
write.csv(quantile, "./DEX_genes.csv") 

