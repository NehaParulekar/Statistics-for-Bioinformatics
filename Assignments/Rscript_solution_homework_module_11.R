# Solutions to Module 11

#problem 1 

# Clustering analysis on the "CCND3 Cyclin D3" gene expression values of the Golub et al. (1999) data

#loading the data

data(golub, package = "multtest")
gol.fac <- factor(golub.cl,levels=0:1, labels= c("ALL","AML"))
grep("CCND3 Cyclin D3",golub.gnames[,2])

# 1(a) Conduct hierarchical clustering using single linkage and Ward linkage. Plot the cluster dendrogram for both fit. Get two clusters from each of the methods. Use function table() to compare the clusters with the two patient groups ALL/AML. Which linkage function seems to work better here?
clus.data <- golub[1042,]

# Conduct hierarchical clustering using single linkage and Ward linkage
# Plot cluster dendrogram
hclust.single<- hclust(dist(clus.data, method="euclidian"), method="single") 
plot(hclust.single,labels=gol.fac) 
hclust.ward <- hclust(dist(clus.data, method="euclidian"), method="ward.D2")
plot(hclust.ward,labels=gol.fac)

# dislpay graphs
par(mfrow=c(1,2))
plot(hclust.single,labels=gol.fac,cex=0.4) 
plot(hclust.ward,labels=gol.fac,cex=0.4)

# Get two clusters from each of the methods
clust.2single <- cutree(hclust.single,2)
clust.2ward <- cutree(hclust.ward,2)

# Use function table() to compare the clusters with the two patient groups ALL/AML
table(gol.fac, clust.2single)
table(gol.fac, clust.2ward)


# 1(b) Use k-means cluster analysis to get two clusters. Use table() to compare the two clusters with the two patient groups ALL/AML.
# Use k-means cluster analysis to get two clusters
clusters.kmean <- kmeans(clus.data, centers=2)
table(gol.fac, clusters.kmean$cluster)


#1(d) Find the two cluster means from the k-means cluster analysis. Perform a bootstrap on the cluster means. Do the confidence intervals for the cluster means overlap? Which of these two cluster means is estimated more accurately?
clusters.kmean$centers
initial <- clusters.kmean$centers
n <- length(clusters.kmean$cluster)
nboot<-1000
boot.cl <- matrix(NA,nrow=nboot, ncol=2) 
for (i in 1:nboot){
  dat.star <- clus.data[sample(1:n,replace=TRUE)]
  cl <- kmeans(dat.star, centers=initial)
  boot.cl[i,] <- cl$centers
}
apply(boot.cl,2,mean)
quantile(boot.cl[,1],c(0.025,0.975))
quantile(boot.cl[,2],c(0.025,0.975))


#1(e) Produce a plot of K versus SSE, for K=1, ., 30. How many clusters does this plot suggest?
K <- (1:30)
SSE<- rep(NA,length(K)) 
for (k in K) {SSE[k] <- kmeans(clus.data, centers=k, nstart = 10)$tot.withinss}
plot(K, SSE, type='o', xaxt='n'); axis(1,at = K, las=2)

# Problem 2 

#Cluster analysis on part of Golub data.
#2(a) Select the oncogenes and antigens from the Golub data. (Hint: Use grep() ).
golub.oncogene <- agrep("^oncogene",golub.gnames[,2])
golub.antigen <- agrep("^antigen",golub.gnames[,2])
golub.oncoangn <- unique(c(golub.onco, golub.angn))
gene.factor <- factor(c(rep("oncogene",length(golub.onco)),rep("antigen",length(golub.angn))))


#2(b) On the selected data, do clustering analysis for the genes (not for the patients). Using K-means and K-medoids with K=2 to cluster the genes. Use table() to compare the resulting two clusters with the two gene groups oncogenes and antigens for each of the two clustering analysis.
library(cluster)
data <- data.frame(golub[golub.oncoangn,])
# k-means
clustermea <- kmeans(data, centers=2)
#k-medeoids
clustermed <- pam(data, k=2)
#oncogenes comparision
oncocomp <- table(gene.factor, clustermea$cluster)
#antigens comparision
antigcomp <- table(gene.factor, clustermed$cluster)
oncocomp 
antigcomp


#2(c) Use appropriate tests (from previous modules) to test the marginal independence in the two by two tables in (b). Which clustering method provides clusters related to the two gene groups?
# using fisher test
fisher.test(oncocomp)
fisher.test(antigcomp)

#(d) Plot the cluster dendrograms for this part of golub data with single linkage and complete linkage, using Euclidean distance.
hclust.single <- hclust(dist(data, method="euclidian"), method="single") 
# Eucledian, Single
plot(hclust.single,labels=gene.factor, cex=0.62) 
hclust.complete <- hclust(dist(data, method="euclidian"), method="complete")
# Eucledian, Complete 
plot(hclust.complete,labels=gene.factor, cex=0.62)


# Problem 3 

# 3(a) Using k-means clustering, produce a plot of K versus SSE, for K=1,., 30. How many clusters appears to be there?

#loading the data
install.packages('ISLR')
library(ISLR)
ncidata<-NCI60$data
ncilabs<-NCI60$labs

#  Produce a plot of K versus SSE, for K=1,..., 30
K <- (1:30) 
SSE.1 <- rep(NA,length(K)) 
for (k in K) {SSE.1[k]<-kmeans(ncidata, centers=k, nstart = 10)$tot.withinss}
plot(K, SSE.1, type='o', xaxt='n'); axis(1,at = K, las=2)


# (b) Do K-medoids clustering (K=7) with 1-correlation as the dissimilarity measure on the data. Compare the clusters with the cell lines. Which type of cancer is well identified in a cluster? Which type of cancer is not grouped into a cluster? According to the clustering results, which types of cancer are most similar to ovarian cancer?

library(cluster)
clusters.7km <- pam(as.dist(1-cor(t(ncidata))),k=7)
table(factor(ncilabs), clusters.7km$cluster)
