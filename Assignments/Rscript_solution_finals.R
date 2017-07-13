# solutions to Finals Module

# Problem 1 

# downloading the file DataPois in the working directory and loading it

getwd()
y <- as.numeric(t(read.table(file="DataPois.txt", header=TRUE)))

# 1(a) What is the sample size n? What is the sample mean Y ?

#sample size n is 
str(y)
#MEan of the sample 
print("mean of th esample is:")
mean(y)

# 1(b) Find the value of MLE ^??? on this data set using numerical method.

llh <- function(x) - sum(log(dpois(y, lambda=exp(x))))
print("the vaue of theta is:")
optim(0, llh)$par


# 1(c)Test the null hypothesis that ??? ???1 at level 0.05, using a bootstrap confidence interval.

n<-length(y)
nboot<-1000
boot.xbar <- rep(NA, nboot)
for (i in 1:nboot) {
  data.star <- y[sample(1:n,replace=TRUE)]
  llh<-function(x) sum(log(dpois(data.star, lambda=exp(x))))
  nllh<- function(x)-llh(x)
  boot.xbar[i]<-optim(1, nllh)$par
}
quantile(boot.xbar,c(0.025,0.975))


# Problem 2

# 2(a) Delete the cancer types with only one or two cases ("K562A-repro", etc.). Keep only the cancer types with more than 3 cases.

library(ISLR)
nci.data<- NCI60$data
nci.labs<- NCI60$labs
data <- NULL
z=1
for (i in 1:64){
  if (sum(nci.labs==(nci.labs[i])) <= 3){
    data[z] <- i
    z=z+1   
  }
}
new.data <- nci.data[-data,]
new.labs <- nci.labs[-data]

# 2(b) Analyze the expression values of the first gene in the data (first column). Does the first gene express differently in different types of cancers? If so, in which pairs of cancer types does the first gene express differently?

gene<- new.data[,1]
# Perform the anova
anova(lm(gene ~ new.labs))
# Performing pairwise t test with fdr adjustment.
pairwise.t.test(gene, new.labs, p.adjust.method = 'fdr')

# 2(c) Check the model assumptions for analysis in part (b). Is ANOVA analysis appropriate here

#loading the library
library(lmtest)
# using shapiro to test for normal
shapiro.test(residuals(lm(gene ~ new.labs)))
# bptest
bptest(lm(gene ~ new.labs), studentize = FALSE)

# 2(d) Apply ANOVA analysis to each of the 6830 genes. At FDR level of 0.05, how many genes express differently among different types of cancer patients?

anova <- apply(new.data, 2, function(x) anova(lm(x ~ new.labs))$Pr[1])
p.fdr <- p.adjust(p=anova, method="fdr")
sum(p.fdr<0.05)


# Problem 3 

# 3(a) Make pairwise scatterplots for all variables in the data set. Which variables appears to be linearly correlated with the life expectancy based on the
# scatterplots?

# loading the data 
data(state)
head(state.x77)
# collecting the information and plotting the scatter plot
data <- as.data.frame(state.x77[,1:8])
pairs(data) 

# 3(b) Conduct a regression analysis different from the example analysis in module 9. We regress the life expectancy on three variables: the per capita income
# (Income), the illiteracy rate (Illiteracy) and mean number of frost days (Frost).

req.data <- as.data.frame(state.x77[,c('Life Exp', 'Income', 'Illiteracy', 'Frost')]) 
names(req.data)<- c('Life.Exp', 'Income', 'Illiteracy', 'Frost')
reg.analysis <- lm(Life.Exp~Income+Illiteracy+Frost, data=req.data)
summary(reg.analysis)

# 3 (c)Find delete-one-cross-validated mean square errors.

n <- dim(req.data)[1]
error <- rep(NA, n)
for (i in 1:n) {
  m <- req.data[-i,]
  reg.ana <- lm(Life.Exp~Income+Illiteracy+Frost, data=m)
  coeffs <- reg.ana$coefficients
  pred <- coeffs[1] + sum(coeffs[-1] * req.data[i,-1])
  error[i] <- (req.data[i,1] - pred)^2
}
error <- mean(error)
print(error)


# Problem 4

# 4(a) Select gene expression data for only the B-cell patients. The analysis in following parts will only use these gene expression data on the B-cell patients.

#load the data 
library(ALL)
data(ALL)
ALL.B <- ALL[,which(ALL$BT %in% c("B","B1","B2","B3","B4"))]
B.data <- exprs(ALL.B)
str(B.data)

# 4(b) Select only those genes whose coefficient of variance (i.e., standard deviation divided by the mean) is greater than 0.2. How many genes are selected?

# take mean for the data
mean <- apply(data,1,mean)
# take standard deviation for the data
std.dev <- apply(data,1,sd)
#genes whose coefficient of variance is greater than 0.2
coeff.var<- subset(std.dev/mean,std.dev/mean>0.2)
new.data<- data[names(coeff.var),]
# no. of genes whose coefficient of variance
sum(std.dev/mean > 0.2)
row.names(new.data)

# 4(d) Conduct a hierarchical clustering analysis with filtered genes in (b). (For uniformity in grading, we ask everyone to use the filter in (b). It may not be your best filter in (c).) How do the clusters compare to the B-stages? How does do the clusters compare to the molecule biology types (in variable ALL$mol.biol)? Provide the confusion matrices of the comparisons, with 4 clusters.

hc <- hclust(dist(t(new.data)),method="complete")
cluster<-cutree(hc,k=4)
table(ALL$BT[1:95,drop=T],cluster)
table(ALL$mol.bio[1:95],cluster)

# 4(e) Draw two heatmaps for the expression data in (d), one for each comparison. Using colorbars to show the comparison types (B-stages or molecule biology types). The clusters reflect which types better: B-stages or molecule biology types?

library(gplots)
color.map <- function(T) {
  if (T=="B") "red"
  else if(T=="B1") "darkblue"
  else if(T=="B2") "yellow"
  else if(T=="B3") "green"
  else "purple"
}
clmy <- as.character(ALL$BT[1:95]) 
patientcolors <- unlist(lapply(clmy, color.map))
heatmap.2(data.new,
          col=greenred(75),
          scale="row",
          dendrogram="column",
          ColSideColors=patientcolors,
          margin=c(3, 12), cexRow=0.5,
          key=FALSE, trace="none", labRow=NA, labCol=NA)

color.map <- function(T) {
  if (T=="ALL1/AF4") "red"
  else if(T=="BCR/ABL") "darkblue"
  else if(T=="E2A/PBX1") "yellow"
  else if(T=="NEG") "green"
  else if(T=="NUP-98") "brown4"
  else "purple"
}
clmy <- as.character(ALL$mol.biol[1:95]) 
patientcolors <- unlist(lapply(clmy, color.map))
dev.off()
heatmap.2(data.new,
          col=greenred(75),
          scale="row",
          dendrogram="column",
          ColSideColors=patientcolors,
          margin=c(3, 12), cexRow=0.5,
          key=FALSE, trace="none", labRow=NA, labCol=NA)

#4(f)We focus on predicting the B-cell differentiation in the following analysis. We merge the last two categories "B3" and "B4", so that we are studying 3 classes: "B1", "B2" and "B34".

clases <- ALL[,which(ALL$BT %in% c("B1","B2","B3","B4"))]$BT
levels(clases)[levels(clases)=="B3"] <- "B34"
levels(clases)[levels(clases)=="B4"] <- "B34"
clases <- droplevels(clases)
library(limma)
linear.model <- model.matrix(~ 0 + factor(levels))
colnames(linear.model) <- c("B1","B2","B34") 
fit <- lmFit(ALL[,which(ALL$BT %in% c("B1","B2","B3","B4"))], linear.model)
fit <- eBayes(fit)
cont.ma <- makeContrasts(B1-B2,B2-B34,B1-B34, levels=clases)
fit1 <- contrasts.fit(fit, cont.ma)
fit1 <- eBayes(fit1)
dim(topTable(fit1, number=Inf, p.value=0.05, adjust.method="fdr"))[[1]]
a = topTable(fit1, number=Inf, p.value=0.05, adjust.method="fdr")
topTable(fit1, number=Inf, p.value=0.05, adjust.method="fdr")[1:10,]
names <- row.names(a)
names

#4(g) Fit SVM and the classification tree on these selected genes in part (f), evaluate their performance with delete-one-cross-validated misclassification rate.

library("hgu95av2.db");library(ALL);data(ALL)
ALLB <-ALL[,ALL$BT %in% c("B1","B2","B3", "B4")];
ALLB$BT <- levels;

pano<- apply(exprs(ALLB), 1, function(x) anova(lm(x ~ ALLB$BT))$Pr[1])
names<- featureNames(ALL)[pano<0.000001] 
symb <-mget(names, env = hgu95av2SYMBOL) 
ALLBTnames <-ALLB[names, ] 
probedat <-as.matrix(exprs(ALLBTnames))
row.names(probedat)<-unlist(symb)

# For SVM

library(e1071)
Bclasses <- factor(ALLBTnames$BT)
exprgene <- t(probedat)
svmest <- svm(Bclasses~exprgene, type = "C-classification", kernel = "linear") #Do svm for classification of Bclasses based on exprgene, use linear kernel.
svmpred <-predict(svmest, exprgene)
table(svmpred, Bclasses)

Bclasses <- factor(ALLBTnames$BT)
exprgene<- t(probedat)
n<-length(Bclasses)
mcr.raw<-rep(NA, n) 
for (i in 1:n) {
  svmest<- svm(exprgene[-i,], Bclasses[-i], type = "C-classification", kernel = "linear")
  svmpred <-predict(svmest, t(exprgene[i,]))
  mcr.raw[i] <-mean(svmpred!=Bclasses[i])
}

mcr.cv<-mean(mcr.raw)
mcr.cv


# For Classification tree

require(rpart)
class <- factor(ALLBTnames$BT)
ctr <- rpart(class~., data = data.frame(t(probedat)))
plot(ctr, branch = 0, margin = 0.1);
text(ctr, digits=3)
rpartpred <- predict(ctr, type="class")
table(rpartpred, class)

class <- factor(ALLBTnames$BT)
exprsgene <- data.frame(t(probedat))
n <- length(class)
p <- 1
nsim <- 1000
tree.raw <- rep(NA,nsim)
for(i in 1:nsim) {
  testID <- sample(n,p,replace=FALSE)
  tr.est <- rpart(class[-testID]~., data=exprsgene[-testID,])
  tr.pred <- predict(ctr, newdata = exprsgene[testID,], type = "class")
  tree.raw[i] <- mean(tr.pred!=class[testID])
}
tree.cv <- mean(tree.raw)
tree.cv


#4(h) We select the genes passing both filters in (b) and (f). How many genes are selected? Redo part (g) on these genes passing both filters.

sum(row.names(new.data) %in% names)
as.numeric(row.names(new.data) %in% names)
new.genes <- row.names(new.data)[row.names(new.data) %in% names]
new.genes

# Problem 5 

install.package('glm')
library(glm)

pois.data = read.table(file ='DataPoisreg.txt',header= TRUE)
attach(poisdata)
m1 = glm(y~x,family=poisson(), data=pois.data)
summary(m1)
b0hat = m1$coeff[[1]]
b1hat = m1$coeff[[2]]

