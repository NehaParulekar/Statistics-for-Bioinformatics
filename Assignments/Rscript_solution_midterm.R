# Solution to Midterm 

# problem 1

rm(list=ls())
# a) Find E(X), E(Y), sd(X) and sd(Y).

# Finding E(X)
X.Range<- c(1,2,3)
f.x<- function(x) 2.469862*(x*exp(-x^2))
f_x<- function(x) f.x(x)*(x %in% X.Range)
E.X<- sum(X.Range*f_x(X.Range))
print("The E(X) is ")
print(E.X)

# Finding E(Y)
f.Y<- function(y) 2*(y*exp(-y^2))*(0<=y & y<=Inf)
E.Y<- integrate(function(y) y*f.Y(y), lower=0, upper=Inf)$value
print("The E(Y) is")
print(E.Y)

# Finding sd(X)
X.Range<-c(1,2,3)
f.x<- function(x) 2.469862*(x*exp(-x^2))
f_x<- function(x) f.x(x)*(x %in% X.Range)
E.X<- sum(X.Range*f_x(X.Range))
# findind variance
Var.X<- sum((X.Range-E.X)^2 * f_x(X.Range)) 
# finding standard deviatioin
sd.X<- sqrt(Var.X)
print("the sd(x) is")
print (sd.X)

# Finding sd(Y)
f.Y<- function(y) 2*(y*exp(-y^2))*(0<=y & y<=Inf)
E.Y<- integrate(function(y) y*f.Y(y), lower=0, upper=Inf)$value
Var.Y<-integrate(function(y) (y-E.Y)^2*f.Y(y), lower=0, upper=Inf)$value
sd.Y<- sqrt(Var.Y)
print("The sd(Y) is")
print(sd.Y)

# b) If X and Y are independent, find E(2X-3Y) and sd(2X-3Y).

# Finding E(2X-3Y)
(2*(E.X) - 3*(E.Y))

# Finding sd(2X-3Y)
Var<- (4*Var.X)+(9*Var.Y)
sd<- sqrt(Var)
print("The sd(2X-3Y) is")
print(sd)


# problem 2 

rm(list=ls())

# Estimate accurately to two decimal places

# X follows standard normal distribution
X<- rnorm(10000, mean=0, sd=1)
# Y follows chi-square distribution
Y<- rchisq(10000, df=4)
# solving the equation
eqn <- (X^2)/(X^2+Y)
E.eqn<-mean(eqn)
print(E.eqn)


# problem 3

rm(list=ls())

# What is the probability that our empirical coverage is greater than 94%?

1-pbinom(940,1000,0.92)


# problem 4 

rm(list=ls())

# Find the value of MLE ^??? on this data set

# download the file, put it in the working directory of your R session.
getwd()
# load it using command
y<-as.numeric(t(read.table(file = "normalData.txt", header=T)))

c(sqrt(sum((y-mean(y))^2)/length(y)))


# problem 5 

rm(list=ls())

# a) Use the t-test to test how many genes have mean expression values greater than 0.6. Use a FDR of 10%.

#loading the data 
data(golub, package = "multtest");
golub.fac <- factor(golub.cl, levels = 0:1, labels = c("ALL", "AML"))
#applying t-test
p.values <- apply(golub, 1, function(x){
  t.test(x, mu=0.6, alternative = c("greater"))$p.value
})
# finding the genes with mean expression greater than 0.6
mean <- p.values<0.05
print(" The genes with mean expression greater than 0.6")
sum(p.values<0.05)
# using FDR of 10%
p.fdr<-p.adjust(p = p.values, method="fdr")
print("the no of genes have mean expression values greater than 0.6 after 10% FDR are")
sum(p.fdr<0.10)

#b) Find the gene names of the top five genes with mean expression values greater than 0.6.
genes <- order(p.fdr, decreasing = FALSE)[1:5]
print("gene names of the top five genes with mean expression values greater than 0.6.")
golub.gnames[genes, 2]


# problem 6

# On the Golub et al. (1999) data set, compare the "GRO3 GRO3 oncogene" (at row 2715) with the "MYC V-myc avian myelocytomatosis viral oncogene homolog" (at row 2302).

rm(list=ls())

# Load the data for problem 6 
data(golub, package = "multtest")
gol.fac <- factor(golub.cl,levels=0:1, labels= c("ALL","AML"))

# a) Draw a histogram of the GRO3 gene expression values.

# row index of GRO3 gene is 2715. getting the data for GRO3
GRO3.data = golub[2715,] 
# to draw a histogram 
hist(GRO3.data)

# b) Draw a scatterplot of the GRO3 gene expression values versus MYC gene expression values, labeled with different colors for ALL and AML patients.

# row index of MYC gene is 2302. getting the data for MYC
MYC.data = golub[2302,] 
# row index of GRO3 gene is 2715. getting the data for GRO3
GRO3.data<- golub[2715,] 

ALL.x<- golub[2715, gol.fac=="ALL"]
AML.x<- golub[2715, gol.fac=="AML"]
ALL.y<- golub[2302, gol.fac=="ALL"]
AML.y<- golub[2302, gol.fac=="AML"]

# to draw a scatterplot
plot(ALL.x,ALL.y,xlab=golub.gnames[2715,2],ylab=golub.gnames[2302,2], col="blue")
points(AML.x,AML.y, col="red")
legend("topright",c("ALL","AML"),col=c("blue","red"),lty=c(1,1))

# c) Use a parametric t-test to check (the alternative hypothesis) if the mean expression value of GRO3 gene is less than the mean expression value of MYC gene.
# applying t-test
t.test(GRO3.data, MYC.data, paired=T, alternative = "less" )

# d) Use a formal diagnostic test to check the parametric assumptions of the t-test. Is the usage of the t-test appropriate here?
# applying shapiro on GRO3 and MYC with their row index 2715 and 2302 respectively
shapiro.test(golub[2715,])
shapiro.test(golub[2302,])

# e)Use a nonparametric test to check (the alternative hypothesis) if the median difference between the expression values of GRO3 gene and the expression values of MYC gene is less than zero
# apply wilcox test
wilcox.test (x= GRO3.data, y= MYC.data, paired=T, alternative="less")

# f) Calculate a nonparametric 95% one-sided upper confidence interval for the median difference between the expression values of GRO3 gene and of MYC gene.

data.diff<-abs(golub[2715,]-golub[2302,])
n<-length(data.diff)
nboot<-1000
boot.xbar<-rep(NA, nboot)
for(i in 1:nboot) {
  data.star <- data.diff[sample(1:n,replace=TRUE)]
  boot.xbar[i]<-median(data.star)
}
quantile(boot.xbar,c(0.95))


# g) Calculate a nonparametric bootstrap 95% one-sided upper confidence interval for the mean difference between the expression values of GRO3 gene and of MYC gene.
data.diff<-abs(golub[2715,]-golub[2302,])
n<-length(data.diff)
nboot<-1000
t.boot = rep(NA,nboot)
for(i in 1:nboot){
  data.star<-data.diff[sample(1:n,replace=TRUE)]
  t.boot[i]= mean(data.star)
}
quantile(t.boot,c(0.95))


# problem 7

# On the Golub et al. (1999) data set, complete the following:

rm(list=ls())

# a) Find the row number of the "HPCA Hippocalcin" gene.

# load the data 
data(golub, package = "multtest")
gol.fac <- factor(golub.cl,levels=0:1, labels= c("ALL","AML"))

# a) Find the row number of the "HPCA Hippocalcin" gene.
# Finding the roe index through grep

grep("HPCA Hippocalcin",golub.gnames[,2])

# b) Find the proportion among ALL patients that the "HPCA Hippocalcin" gene is negatively expressed (expression value<0).
HPCA.data = golub[118, gol.fac=="ALL"] 
mean(HPCA.data<0)

# c) We want to show that "HPCA Hippocalcin" gene is negatively expressed in at least half of the population of the ALL patients.
# using binomial test 
p.all<-sum(golub[118,gol.fac=="ALL"]<0)
binom.test(x=p.all,n=27,p= 0.5, alternative="greater")

# d) Find a 95% confidence interval for the difference of proportions in the ALL group versus in the AML group of patients with negatively expressed "HPCA Hippocalcin" gene.
# getting the row index for HPCA
grep("HPCA Hippocalcin",golub.gnames[,2])
# getting the data 
x.ALL <- golub[118,gol.fac=="ALL"]
nALL <- length(x.ALL)
x.AML<- golub[118,gol.fac=="AML"]
nAML <- length(x.AML)
nboot<-1000
boot.diff <- rep(NA,nboot)
for (i in 1:nboot) {
  data.ALL <- x.ALL[sample(1:nALL,replace=TRUE)]
  data.AML <- x.AML[sample(1:nAML,replace=TRUE)]
  data.diff <- data.ALL-data.AML
  boot.diff[i] <- mean(data.diff)
}
quantile(boot.diff,c(0.025,0.975))
