# solutions to homework 13

# Problem 1

# 1(a) Define an indicator variable IsB such that IsB=TRUE for B-cell patients and IsB=FALSE for T-cell patients.

#loading the data
library("ALL")
data(ALL)
install.packages(c('rpart','ROCR','VGAM'))
library("hgu95av2.db")
library(ROCR) 

#defining indicator variable
IsB <- factor(ALL$BT %in% c("B","B1","B2","B3","B4"))

# 1(b) Use two genes "39317_at" and "38018_g_at" to fit a classification tree for IsB. Print out the confusion matrix. Plot ROC curve for the tree.

genes <- c("39317_at", "38018_g_at")
#get the gene names
gene.names <- mget(genes, env = hgu95av2SYMBOL)
#keep only the gene passing filter
ALL.BTnames <- ALL[genes, ]
#save expression values in 'values'
values <- as.matrix(exprs(ALL.BTnames))
#change row names to gene names
row.names(values) <- unlist(gene.names)
#load library
require(rpart) 
#The true classes: B,B1,B2,B3,B4 stages
stages <- factor(IsB)
#Fit the tree on data set 'values'. The transpose is needed since we are classifying the patients (rows).
tree.fit <-rpart(stages~., data = data.frame(t(values)))
#Plot the tree with V-shaped branch (=0), leave margins for later text additions.
plot(tree.fit, branch=0, margin=0.1, cex=0.7)
#Add text for the decision rules on the tree. use 3 digits for numbers
text(tree.fit, digits=3)
#predicted classes from the tree
predicted.part<- predict(tree.fit, type="class")
#confusion matrix compare predicted versus true classes
table(predicted.part, stages)
predicted.prob <- predict(tree.fit,type="prob")[,2]
#using predicted.prob to predict the response (which need values of T/F, thus we apply ==???B1???)
pred <- prediction(predicted.prob, IsB=="TRUE")
#compute tpr and fpr for pred
perf <- performance(pred, "tpr", "fpr" )
#Plot tpr versus fpr, i.e., ROC curve
plot(perf)


# 1(c) Find its empirical misclassification rate (mcr), false negative rate (fnr) and specificity. Find the area under curve (AUC) for the ROC curve.
# misclassification = 
print("the empirical misclassification rate is")
(11+2)/128
# fpr =
print("the false positive rate is")
2/(31+2)
# fnr =
print("the false negetive rate is")
11/(84+11)
# Specificity = tnr =
print("the specificity is")
31/(2+31)
# AUC = 0.922807
print("the area under curve AUC is:")
performance(pred,"auc")


# 1(d) Use 10-fold cross-validation to estimate its real false negative rate (fnr). What is your estimated fnr?

# loading the data 
library("ALL")
data("ALL")
probe<-c("39317_at","38018_g_at")
exprs.data<-exprs(ALL)[probe,]
IsB<-c(rep("TRUE",sum(ALL$BT %in% c("B","B1","B2","B3","B4"))),rep("FALSE",sum(ALL$BT %in% c("T","T1","T2","T3","T4"))))
data.lgr<-data.frame(IsB,t(exprs.data))
require(caret)
data.lgr<-data.frame(IsB,t(exprs.data))
n<-dim(data.lgr)[1]
index<-1:n
K<-10
flds<-createFolds(index,k=K)
fnr.true.raw<-rep(NA,K)
for(i in 1:K){
  testID<-flds[[i]]
  data.tr<-data.lgr[-testID,]
  data.test<-data.lgr[testID,]
  fit.lgr <- glm(IsB~., family=binomial(link='logit'), data=data.tr)
  pred.prob<-predict(fit.lgr,newdata=data.test,type="response")
  pred.B <- (pred.prob>0.5)
  fnr.true.raw[i] <- (sum(pred.B != "TRUE" & IsB[testID] == "TRUE"))/(sum(IsB[testID] == "TRUE"))
}
fnr.true<-mean(fnr.true.raw)
fnr.true

#1(e) Do a logistic regression, using genes "39317_at" and "38018_g_at" to predict IsB. Find an 80% confidence interval for the coefficient of gene "39317_at".
group.names <- c("39317_at", "38018_g_at")
exprs.data <- exprs(ALL)[group.names,]
data.lgr <- data.frame(IsB, t(exprs.data))
#Fit the generalized linear model, for binomial (0/1 responses) with logit link.That means logistic regression.
fit.data <- glm(IsB~., family=binomial(link='logit'), data=data.lgr)
#get prediction probability for all cases in data set
probability <- predict(fit.data, data=data.lgr$exprs.data, type="response")
pred.B1 <- factor(probability> 0.5, levels=c(TRUE,FALSE), labels=c("B","not_B"))
#Substitute names in "pred.B1" variable.
IsB1 <- factor(IsB, levels=c(TRUE,FALSE), labels=c("B","not_B"))
#Compare classified  classes in a confusion matrix
table(pred.B1, IsB1)

confint(fit.data, level=0.8)

#1(f) Use n-fold cross-validation to estimate misclassification rate (mcr) of the logistic regression classifier. What is your estimated mcr?

#loading the data
install.packages('caret')
require(caret);
ALLB <- ALL[,which(ALL$BT %in% c("B","B1","B2","B3","B4"))]
data.lgr <- data.frame(IsB,t(expr.data))
n <- dim(data.lgr)[1]
index <- 1:n 
K <- n
folds <- createFolds(index, k=K)
mcr.cv.raw <- rep(NA, K)
for (i in 1:K) {
  test.ID <- folds[[i]] 
  data.tr <- data.lgr[-test.ID,]
  data.test <- data.lgr[test.ID,] 
  fit.log <- glm(IsB~., family=binomial(link='logit'), data=data.tr)
  predicted.prob <- predict(fit.log, newdata=data.test, type="response")
  pred <- (predicted.prob> 0.5)
  mcr.cv.raw[i] <- sum(pred!=data.test$IsB)/length(pred)
}
mcr.cv <- mean(mcr.cv.raw) 
mcr.cv 

#1(g) Conduct a PCA on the scaled variables of the whole ALL data set (NOT just the two genes used above). We do this to reduce the dimension in term of genes (so this PCA should be done on the transpose of the matrix of expression values). To simply our future analysis, we use only the first K principal components (PC) to represent the data. How many PCs should be used? Explain how you arrived at your conclusion. Provide graphs or other R outputs to support your choice.

#Apply PCA on ALL data set.transpose of the matrix of eprs values 
pca.ALL<-prcomp(t(exprs(ALL)), scale=TRUE) 
#print out summary of the PCA ALL
summary(pca.ALL)
PropVar <- summary(pca.ALL)$importance[2,] 
plot(1:length(PropVar), PropVar, xlab='number of principal components', ylab='proportion of variance explained',cex=0.3) 

#1(h) Do a SVM classifier of IsB using only the first five PCs. (The number K=5 is fixed so that we all use the same classifier. You do not need to choose this number in the previous part (g).) What is the sensitivity of this classifier?

#We now apply the support vector machine using svm() function 
install.packages('e1071');
library(e1071);
data.pca <- pca.ALL$x[,1:5]
#Do svm for classification of B.stage based on exprs.gene, use linear kernel.
smv <- svm(IsB~data.pca,type="C-classification",kernel="linear")
#use fit to get predicted classes
predicted.svm <- predict(smv,data.pca)
tpr.svm <- mean(predicted.svm==IsB&IsB==TRUE)/(mean(predicted.svm==IsB&IsB==TRUE)+mean(predicted.svm!=IsB&IsB==TRUE))
print("tpr.svm")
tpr.svm

#1(i) Use leave-one-out cross-validation to estimate misclassification rate (mcr) of the SVM classifier. Report your estimate.

#size of the whole sample
n <- length(IsB) 
# A vector to save misclassification rates (mcr) for each delete-one validation
mcrcvraw <-rep (NA,n) 
for (i in 1:n) {
  #delete each observation from training, then use it for testing.
  svm <- svm(data.pca[-i,],IsB[-i],type="C-classification",kernel="linear")
  #get prediction. transpose t() is used to make the vector into 1 by ncol matrix
  predicted.svm <- predict(svm,t(data.pca[i,])) 
  #misclassification proportion
  mcrcvraw[i] <-mean(predicted.svm!=IsB[i])
}
#average the mcr over all n rounds.

mcrcv<- mean(mcrcvraw)#average the mcr over all n rounds.
mcrcv



# Problem 2 

# K=1

install.packages(c('rpart','ROCR','VGAM'))
library(VGAM)
c("K = 1")
print("Classifier :          MCR, leave-one-out MCR")
pca.iris <- prcomp(iris[,1:4], scale=TRUE)
Species <- iris$Species
data.pca <- pca.iris$x[,1,drop = F]
n <- length(Species)
iris2 <- data.frame(Species, data.pca)
# Logistic Regression
log.reg<- vglm(Species~., family=multinomial, data=iris2) 
# get prediction
predicted.prob <- predict(log.reg, iris2[,-1,drop=F], type="response")
#Assign to the class with largest
predicted.lgr <- apply(predicted.prob, 1, which.max) 
#relabel 1,2,3 by species names
predicted.lgr <- factor(predicted.lgr, levels=c("1","2","3"), labels=levels(iris2$Species)) 
#misclassification rate
mcr.lgr <- mean(predicted.lgr!=iris2$Species)
### leave-one-out cross validation
#A vector to save mcr validation
mcr.cv.raw<-rep(NA, n) 
for (i in 1:n) {
  #fit logistic
  fit.log <- vglm(Species~., family=multinomial, data=iris2[-i,]) 
  #get prediction probability
  predicted.prob <- predict(fit.log, iris2[i,-1, drop=F], type="response") 
  #Assigning the class
  pred <- apply(predicted.prob, 1, which.max) 
  pred <- factor(pred, levels=c("1","2","3"), labels=levels(iris2$Species))
  #check misclassification
  mcr.cv.raw[i] <- mean(pred!=Species[i]) 
}
#average the mcr over all n rounds.
mcr.cv <- mean(mcr.cv.raw) 
c(mcr.lgr, mcr.cv)

# SVM
#train SVM
library(svm)
iris2.svm <- svm(data.pca, Species, type = "C-classification", kernel = "linear") 
#get SVM prediction.
svm.pred<- predict(iris2.svm , data.pca) 
#misclassification rate
mcr.svm<- mean(svm.pred!=Species) 
#A vector to save mcr validation
mcr.cv.raw<-rep(NA, n)
for (i in 1:n) {
  #train SVM without i-th observation
  svm.est <- svm(data.pca[-i,], Species[-i], type = "C-classification", kernel =
                  "linear") 
  #predict i-th observation.
  svm.pred<- predict(svm.est, t(data.pca[i,]))
  #misclassification rate
  mcr.cv.raw[i]<- mean(svm.pred!=Species[i]) 
}
#average the mcr over all n rounds.
mcr.cv<-mean(mcr.cv.raw) 
c(mcr.svm, mcr.cv)

# Classification Tree
fit <- rpart(Species ~ ., data = iris2, method = "class")
pred.tr<-predict(fit, iris2, type = "class")
mcr.tr <- mean(pred.tr!=Species) 
mcr.cv.raw <- rep(NA, n) 
for (i in 1:n) {
  fit.tr <- rpart(Species ~ ., data = iris2[-i,], method = "class") 
  pred <- predict(fit.tr, iris2[i,], type = "class")
  mcr.cv.raw[i] <- mean(pred!=Species[i]) 
}
mcr.cv<-mean(mcr.cv.raw)
c(mcr.tr, mcr.cv)

# K =2,3,4
data(iris)
pca.iris <- prcomp(iris[,1:4], scale=TRUE)
Species <- iris$Species
Analysis <- function(k){
  print(c("When K=",k))
  data.pca <- pca.iris$x[,1:k]
  n <- length(Species)
  iris2 <- data.frame(Species, data.pca)
  # Logistic Regression
  logis.reg<- vglm(Species~., family=multinomial, data=iris2) 
  predicted.prob <- predict(iris2.lgr, iris2[,-1], type="response") 
  predicted.lgr <- apply(predicted.prob, 1, which.max) 
  predicted.lgr <- factor(predicted.lgr, levels=c("1","2","3"), labels=levels(iris2$Species)) 
  mcr.lgr <- mean(predicted.lgr!=iris2$Species) 
  mcr.cv.raw<-rep(NA, n) 
  for (i in 1:n) {
    fit.log <- vglm(Species~., family=multinomial, data=iris2[-i,]) 
    predicted.prob <- predict(fit.log, iris2[i,-1], type="response") 
    pred <- apply(predicted.prob, 1, which.max) 
    pred <- factor(pred, levels=c("1","2","3"), labels=levels(iris2$Species))
    mcr.cv.raw[i] <- mean(pred!=Species[i]) 
  }
  mcr.cv <- mean(mcr.cv.raw) 
  print(c("For Logistic regression: emperical mcr, leave one out mcr", mcr.lgr, mcr.cv))
  
  # SVM
  iris2.svm <- svm(data.pca, Species, type = "C-classification", kernel = "linear") 
  svm.pred<- predict(iris2.svm , data.pca) 
  mcr.svm<- mean(svm.pred!=Species) 
  mcr.cv.raw<-rep(NA, n) 
  for (i in 1:n) {
    svm.est <- svm(data.pca[-i,], Species[-i], type = "C-classification", kernel =
                    "linear") 
    svm.pred<- predict(svm.est, t(data.pca[i,])) 
    mcr.cv.raw[i]<- mean(svm.pred!=Species[i]) 
  }
  mcr.cv<-mean(mcr.cv.raw) 
  print(c("For SVM: emperical mcr, leave one out mcr", c(mcr.svm, mcr.cv)))
  
  # Classification Tree
  fit <- rpart(Species ~ ., data = iris2, method = "class")
  pred.tr<-predict(fit, iris2, type = "class")
  mcr.tr <- mean(pred.tr!=Species) 
 
  mcr.cv.raw <- rep(NA, n) 
  for (i in 1:n) {
    fit.tr <- rpart(Species ~ ., data = iris2[-i,], method = "class") 
    pred <- predict(fit.tr, iris2[i,], type = "class")
    mcr.cv.raw[i] <- mean(pred!=Species[i]) 
  }
  mcr.cv<-mean(mcr.cv.raw) 
  print(c("For Classification Tree: emperical mcr, leave one out mcr", c(mcr.tr, mcr.cv)))
}

doClassifier(2)
doClassifier(3)
doClassifier(4) 

