#homework solution module 12


# Problem 1

# load the data 
library("ALL"); 
data(ALL);data <- exprs(ALL)

# 1(a)Define an indicator variable ALL.fac such that ALL.fac=1 for T-cell patients and ALL.fac=2 for B-cell patients

ALL.fac <- factor(ALL$BT %in% c("B","B1","B2","B3","B4"), labels=c("1","2"))

# 1(b)Plot the histograms for the first three genes' expression values in one row.

gene.1 <- data[1,] 
gene.2 <- data[2,] 
gene.3 <- data[3,]
par(mfrow=c(1,3))
hist(gene.1, main = "Gene 1")
hist(gene.2, main = "Gene 2")
hist(gene.3, main = "Gene 3")

# 1(c) Plot the pairwise scatterplots for the first five genes.

gene.4 <- data[4,]
gene.5 <- data[5,]
pairs(cbind(gene.1, gene.2, gene.3, gene.4, gene.5))

# 1(d)Do a 3D scatterplot for the genes "39317_at", "32649_at" and "481_at", and color according to ALL.fac (give different colors for B-cell versus T-cell
# patients). Can the two patient groups be distinguished using these three genes?

install.packages("scatterplot3d")
require(scatterplot3d)
par(mfrow=c(1,1))
X <- rbind(exprs(ALL[c("39317_at","32649_at","481_at"),]))
scatterplot3d(t(X),color=ALL.fac)

# 1(e) Do K-means clustering for K=2 and K=3 using the three genes in (d). Compare the resulting clusters with the two patient groups. Are the two groups discovered by the clustering analysis?

cluster1 <- kmeans(t(X),centers=2,nstart=10)
table(ALL.fac,cluster1$cluster)
cluster2 <- kmeans(t(X),centers=3,nstart=10)
table(ALL.fac,cluster2$cluster)

# 1(f) Carry out the PCA on the ALL data set with scaled variables. What proportion of variance is explained by the first principal component? By the second principal component?

P.ALL <- prcomp(data, scale=TRUE)
summary(P.ALL)

# 1(g)Do a biplot of the first two principal components. Observe the pattern for the loadings. What info is the first principal component summarizing?

par(mfrow=c(1,1))
biplot(P.ALL, cex=0.5)

# 1(h)For the second principal component PC2, print out the three genes with biggest loadings and the three genes with smallest loadings. 

gene.order<- order(P.ALL$x[,2], decreasing=T) 
dimnames(data)[[1]][[gene.order[1]]]
dimnames(data)[[1]][[gene.order[2]]]
dimnames(data)[[1]][[gene.order[3]]]
dimnames(data)[[1]][[gene.order[12623]]]
dimnames(data)[[1]][[gene.order[12624]]]
dimnames(data)[[1]][[gene.order[12625]]]

# 1(i) Find the gene names and chromosomes for the gene with biggest PC2 value and the gene with smallest PC2 value.

library("annotate")
library("hgu95av2.db")
genename<-get("481_at", env=hgu95av2GENENAME)
chromosomes<-get("481_at", env=hgu95av2CHRLOC)
genename
chromosomes
genenamelow<-get("39317_at", env=hgu95av2GENENAME)
chromosomeslow<-get("39317_at", env=hgu95av2CHRLOC)
genenamelow
chromosomeslow



# problem 2 

# 2(a) Create a data set consisting of the first four numerical variables in the iris data set (That is, to drop the last variable Species which is categorical). Then make a scaled data set that centers each of the four variables (columns) to have mean zero and variance one.

iris.data <- iris[1:4]
mean <- mean(iris.data[,1])
sd <- sd(iris.data[,1])
Sepal.Length <- NULL
for (i in 1:150){Sepal.Length[i] <- (iris.data[i,1]-mean)/sd}

mean <- mean(iris.data[,2])
sd <- sd(iris.data[,2])
Sepal.Width <- NULL
for (i in 1:150){Sepal.Width[i] <- (iris.data[i,2]-mean)/sd}

mean <- mean(iris.data[,3])
sd <- sd(iris.data[,3])
Petal.Length <- NULL
for (i in 1:150){Petal.Length[i] <- (iris.data[i,3]-mean)/sd}

mean <- mean(iris.data[,4])
sd <- sd(iris.data[,4])
Petal.Width <- NULL
for (i in 1:150){Petal.Width[i] <- (iris.data[i,4]-mean)/sd}

scaled.data <- cbind(Sepal.Length, Sepal.Width, Petal.Length, Petal.Width)
xxx=data.frame(Sepal.Length, Sepal.Width, Petal.Length, Petal.Width)

# 2(b) Calculate the correlations between the columns of the data sets using the cor() function. Show that these correlations are the same for scaled and the unscaled data sets.

cor.scaled <- cor(scaled.data)
cor.unscaled <- cor(iris.data)
cor.scaled
cor.unscaled
all.equal(cor.scaled, cor.unscaled)

# 2(c) Calculate the Euclidean distances between the columns of the scaled data set using dist() function. Show that the squares of these Euclidean distances are proportional to the (1-correlation)s. What is the value of the proportional factor here?

euclidian.dist <- dist(t(scaled.data), method="euclidean")
euclidian.sqdist<-euclidian.dist^2
euclidian.sqdist
cor.data<-as.dist(1-cor(scaled.data))
cor.data
prop.factor<-euclidian.sqdist/cor.data
prop.factor

# 2(d) Show the outputs for doing PCA on the scaled data set and on the unscaled data set. (Apply PCA on the two data sets with option "scale=FALSE". Do NOT use option "scale=TRUE", which will scale data no matter which data set you are using.) Are they the same?

pca.unscaled <- prcomp(iris.data, scale=FALSE)
pca.unscaled
pca.scaled <- prcomp(scaled.data, scale=FALSE)
pca.scaled

# 2(e) What proportions of variance are explained by the first two principle components in the scaled PCA and in the unscaled PCA?

summary(pca.unscaled) 
summary(pca.scaled) 

# 2(f) Find a 90% confidence interval on the proportion of variance explained by the second principal component.

data <- scaled.data
p <- ncol(data)
n <- nrow(data)
nboot<-1000
sdevs <- array(dim=c(nboot,p))
for (i in 1:nboot) {
  dat.star <- data[sample(1:n,replace=TRUE),]
  sdevs[i,] <- prcomp(dat.star)$sdev
}
quantile(sdevs[,2], c(0.05,0.95))

