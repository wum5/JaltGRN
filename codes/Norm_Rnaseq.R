setwd("./")
rawdata <- read.csv("pnas.1402835111.sd01.csv", head=TRUE, sep=",")
# head(rawdata)


## Load the packages
source("https://bioconductor.org/biocLite.R")
biocLite("edgeR")
install.packages("statmod")
install.packages("locfit")
install.packages("xlsx")

library(edgeR)
library(statmod)
library(locfit)
library(xlsx)

## Load the data
ly <- DGEList(counts=rawdata[,2:22], genes=rawdata[,1])
pn <- DGEList(counts=rawdata[,23:43], genes=rawdata[,1])
hb <- DGEList(counts=rawdata[,44:64], genes=rawdata[,1])



## Lyco's Data
## Order transcripts by counts number
nrow(ly)
o <- order(rowSums(ly$counts), decreasing=TRUE)
ly <- ly[o,]
nrow(ly)

## Filtering out lowly expressed genes
keep <- rowSums(cpm(ly)>1) > 10
ly <- ly[keep,]
#summary(cpm(ly))
nrow(ly)

## Estimate library size
ly$samples$lib.size <- colSums(ly$counts)
ly <- calcNormFactors(ly, method="TMM")
ly$samples

group <- factor(c(rep("P3m",3), rep("P4d", 3), rep("P4p", 3), 
		rep("P5d", 3), rep("P5p", 3), rep("P6d", 3), rep("P6p", 3)))
data.frame(Sample=colnames(ly), group)
design <- model.matrix(~group)
rownames(design) <- colnames(ly)

ly <- estimateDisp(ly, design, robust=TRUE)
ly$common.dispersion
#plotBCV(ly)

fit1 <- glmQLFit(ly, design)
ltr1 <- glmQLFTest(fit1, coef=2:7)
top1 <- topTags(ltr1, n=30000)
write.csv(top1, "./Lyco_DEX.csv") 



## Penn's Data
## Order transcripts by counts number
nrow(pn)
o <- order(rowSums(pn$counts), decreasing=TRUE)
pn <- pn[o,]

## Filtering out lowpn expressed genes
keep <- rowSums(cpm(pn)>1) > 10
pn <- pn[keep,]
#summary(cpm(pn))
nrow(pn)

## Estimate library size
pn$samples$lib.size <- colSums(pn$counts)
pn <- calcNormFactors(pn, method="TMM")
pn$samples

group <- factor(c(rep("P3m",3), rep("P4d", 3), rep("P4p", 3), 
                  rep("P5d", 3), rep("P5p", 3), rep("P6d", 3), rep("P6p", 3)))
data.frame(Sample=colnames(pn), group)
design <- model.matrix(~group)
rownames(design) <- colnames(pn)

pn <- estimateDisp(pn, design, robust=TRUE)
pn$common.dispersion
#plotBCV(pn)

fit2 <- glmQLFit(pn, design)
ltr2 <- glmQLFTest(fit2, coef=2:7)
top2 <- topTags(ltr2, n=30000)
write.csv(top2, "./Penn_DEX.csv") 




## Harb's Data
## Order transcripts by counts number
nrow(hb)
o <- order(rowSums(hb$counts), decreasing=TRUE)
hb <- hb[o,]

## Filtering out lowhb expressed genes
keep <- rowSums(cpm(hb)>1) > 10
hb <- hb[keep,]
#summary(cpm(hb))
nrow(hb)

## Estimate library size
hb$samples$lib.size <- colSums(hb$counts)
hb <- calcNormFactors(hb, method="TMM")
hb$samples

group <- factor(c(rep("P3m",3), rep("P4d", 3), rep("P4p", 3), 
                  rep("P5d", 3), rep("P5p", 3), rep("P6d", 3), rep("P6p", 3)))
data.frame(Sample=colnames(hb), group)
design <- model.matrix(~group)
rownames(design) <- colnames(hb)

hb <- estimateDisp(hb, design, robust=TRUE)
hb$common.dispersion
#plotBCV(hb)

fit3 <- glmQLFit(hb, design)
ltr3 <- glmQLFTest(fit3, coef=2:7)
top3 <- topTags(ltr3, n=30000)
write.csv(top3, "./Harb_DEX.csv") 



### Merge the three datasets
lyco <- read.csv("Lyco_DEX.csv", head=TRUE, sep=",", colClasses=c("NULL",NA,NA,NA))
penn <- read.csv("Penn_DEX.csv", head=TRUE, sep=",", colClasses=c("NULL",NA,NA,NA))
harb <- read.csv("Harb_DEX.csv", head=TRUE, sep=",", colClasses=c("NULL",NA,NA,NA))
lyco$genes <- gsub('Solyc', 'lyc_Solyc', lyco$genes)
penn$genes <- gsub('Solyc', 'pen_Solyc', penn$genes)
harb$genes <- gsub('Solyc', 'har_Solyc', harb$genes)
merged = rbind(lyco, penn, harb)
head(merged)
tail(merged)
sorted <- merged[order(merged$FDR), ]
cutoff <- nrow(sorted)/20
quantile = head(sorted, cutoff)
write.csv(quantile, "./DEX_genes.csv") 


