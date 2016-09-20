## Gene co-expression network (Fig.2A in Ichihashi et al., PNAS 2014)

## Meng Wu's notes:
## The authors selected the cut-off values based on 99.99999999% quantiles of normal distribution of transformed z scores (un-adjusted)
## from the Spearman Pearson Correlation values. The the correlations were calculated on the data from all the three species.
## Also, the input data are the normalized read counts of genes in the Cluster 2 from their SOM clustering results.

setwd("/Users/mengwu/Documents/Leaf_GRN/Fig2A")
library(reshape)
library(gplots)
library(RColorBrewer)

## read in data
counts=read.csv ("dev.genes.mm.csv", row.names=1) # this data is normalized read counts using EdgeR 2.11
dim(counts) # [1] 539   63
samples=colnames(counts)
genes=rownames(counts)

### check histogram
Y.u=t(log(counts+(1e-6)))
par(mfrow=c(2,1))
hist(as.matrix(counts), breaks=50, main="edge R normalized counts")
hist(as.vector(Y.u), breaks=50, main="logarithm of the normalized counts")
par(mfrow=c(1,1))

### make correlation table
corr.u=cor(Y.u)
dim(corr.u)  # [1] 539  539
corr.s=cor(Y.u, method="spearman") # Spearman's rank correlation coefficient is more robust to outlier and nonlinear relationships.
hist(corr.s[upper.tri(corr.s)]) # masked values above diagonal

### fishers-z transformation (to make the sample correlation more comparable), 
n=63 ### this is the sample size 
z.s= sqrt(n-3)*0.5*log((1+corr.u)/(1-corr.u))
summary(z.s[upper.tri(z.s)])
hist(z.s[upper.tri(z.s)]) # looks normal

### cut off 
thre.z=qnorm(0.9999999999)  ## normal quanitle 
adjcent.z=abs(z.s)>thre.z  ## symmetric ajacency matrix: 1: there is an edge; 0 : there is no edge 
diag(adjcent.z)=0  ## genes do not connect themselves in the network
rownames(adjcent.z)=rownames(corr.u)
colnames(adjcent.z)=colnames(corr.u)
sum(adjcent.z)/2 ## 18982 edges 

#save adjacency matrix
write.csv(adjcent.z, "pseudo_adjacency_matrix.csv")

### community detection methods on the subgraph
## set up graph
index=rowSums(adjcent.z)>0
weight.adjcent.z=adjcent.z[index,index]
library(igraph)
g.temp=graph.adjacency(weight.adjcent.z, mode="undirected", diag=FALSE)

## fastgreedy.community
community.fastgreedy=fastgreedy.community(g.temp)
community.fastgreedy
# Graph community structure calculated with the fast greedy algorithm
# Number of communities (best split): 15 
# Modularity (best split): 0.116111 (above 0.3 is a good indicator of significant community structure in a network)

table(community.fastgreedy$membership) #size of each cluster
#    1   2   3   4   5   6   7   8   9  10  11  12  13  14  15 
#  114 149  58  32  46   6   4   4   3   3   2   2   2   2   2 

## betweenness and hub
## betweenness: The vertex and edge betweenness are (roughly) defined by the number of geodesics (shortest paths) going through a vertex or an edge.
hist(betweenness(g.temp))
b <- betweenness(g.temp, normalized=TRUE)

## extract No. of edges
df.z.g=rowSums(weight.adjcent.z)
hub <- df.z.g

## community
c <- community.fastgreedy$membership

key <- cbind(b, hub, c)
write.csv(data.frame(key),"key.genes.csv")

### visualization 
v.size=rep(2,length(V(g.temp)))
plot(g.temp, vertex.label=NA, vertex.size=v.size)
V(g.temp)$color <- "gray57"
V(g.temp)[community.fastgreedy$membership==1]$color <- "mediumturquoise"
V(g.temp)[community.fastgreedy$membership==2]$color <- "yellow3"
V(g.temp)[community.fastgreedy$membership==3]$color <- "lightpink2"
V(g.temp)[community.fastgreedy$membership==4]$color <- "deepskyblue4"
V(g.temp)[community.fastgreedy$membership==5]$color <- "black"
V(g.temp)[df.z.g>120]$color <- "darkred" # hub genes
v.label=rep("",length(V(g.temp)))
#v.label=V(g.temp)$name  # if you want to put gene name
#v.size[V(g.temp)$name %in% "BOPa"]=4 # if you want to change size of specific nodes
V(g.temp)$shape <- "circle"
pdf("Leaf_gene_coexpression_network_ver4name.pdf", useDingbats=FALSE) 
plot(g.temp, vertex.size=v.size, vertex.frame.color=NA, vertex.label=v.label, 
     vertex.label.cex=0.05,edge.color="gray57", edge.width =0.4, layout=layout.fruchterman.reingold(g.temp))
dev.off()

### plot Hierarchial Clustering (by Meng)
counts <- read.csv("dev.genes.mm.csv", row.names=1)
colnames(counts)<-c(rep("H3",3), rep("H4A",3), rep("H4B",3), rep("H5A",3), rep("H5B",3), rep("H6A",3), rep("H6B",3),
                    rep("L3",3), rep("L4A",3), rep("L4B",3), rep("L5A",3), rep("L5B",3), rep("L6A",3), rep("L6B",3),
                    rep("P3",3), rep("P4A",3), rep("P4B",3), rep("P5A",3), rep("P5B",3), rep("P6A",3), rep("P6B",3))
dim(counts) # [1] 539   63
samples=colnames(counts)
genes=rownames(counts)
Y.u=log(counts+1,2)
dd=as.dist(1-cor(t(Y.u)))

pdf("Hierarchial Clustering.pdf", useDingbats=FALSE) 
hc.complete=hclust(dd, method="complete")
plot(hc.complete, main="Complete Linkage", xlab="", ylab="", sub="", cex=.1)
groups <- cutree(hc.complete, 10)
rect.hclust(hc.complete, k=10, border="red")
dev.off()

pdf("Heat Map.pdf", useDingbats=FALSE) 
x <- as.matrix(Y.u)
dim(x) # [1] 539   63
## Row clustering (adjust here distance/linkage methods to what you need!)
hr <- hclust(as.dist(1-cor(t(Y.u), method="pearson")), method="complete")
## Column clustering (adjust here distance/linkage methods to what you need!)
hc <- hclust(as.dist(1-cor(Y.u, method="spearman")), method="complete")
heatmap.2(x, keysize=1.2, col=brewer.pal(11,"RdBu"), labRow = FALSE, Rowv=as.dendrogram(hr), margins = c(3.5, 1.5),
          Colv=as.dendrogram(hc), scale="row", density.info="none", trace="none", RowSideColors=as.character(groups))
dev.off()

