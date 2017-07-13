#solutions to homework module 7

# prolem 1

# a) Use FDR adjustments at the 0.05 level. How many genes are expressed
# higher in the ALL group?
rm(list=ls())
# load the data 
data(golub, package='multtest')
gol.fac<-factor(golub.cl, level=0:1, labels=c("ALL","AML"))

#for Wilcoxon two-sample tests
wilcox.data=NULL
for (i in 1:3051){
  wilcox.data[i]<-wilcox.test (golub[i,] ~gol.fac, paired=F, alternative="greater")$p.value
}

# to find the genes expressed higher in ALL group

gene.exp<-wilcox.data<0.05
sum(gene.exp)

fdr.wilcox<-p.adjust(p=wilcox.data, method="fdr")
sum(fdr.wilcox<0.05)

# order non-fdr and fdr

non.fdr<- order(wilcox.data, decreasing=FALSE)
print("Gene names for genes with top 3 smallest p-values before FDR adjustment")
golub.gnames[non.fdr[1:3],2]

fdr.AML<- order(fdr.wilcox, decreasing=FALSE)
print("Gene names for genes with top 3 smallest p-values after FDR adjustment")
golub.gnames[fdr.AML[1:3],2]

# finding difference between means 
ALL.mean = apply(golub[, gol.fac=="ALL"], 1, mean)
AML.mean = apply(golub[, gol.fac=="AML"], 1, mean)
Difference<- ALL.mean - AML.mean
diff.order<- order(Difference, decreasing=TRUE)
print("Largest difference")
golub.gnames[diff.order[1:3],2]


# problem 2

# For the Golub et al. (1999) data set, apply the Shapiro-Wilk test of normality to
#every gene's expression values in the AML group. How many genes do not pass
# the test at 0.05 level with FDR adjustment? Please submit your R script with the
# answer.

rm(list=ls())

# load the data 
data(golub, package='multtest')
gol.fac<-factor(golub.cl, level=0:1, labels=c("ALL","AML"))

#applying the test
 shapiro.test<- apply (golub[, gol.fac=="AML"], 1, function(x) {
   shapiro.test(x)$p.value })

# get fdr p values
fdr<- p.adjust(p=shapiro.test, method="fdr")

# calculating and printing the number of genes that failed test

print("the genes do not pass the test at 0.05 level with FDR adjustment")
print(sum(fdr<0.05))


# problem 3

# Gene "HOXA9 Homeo box A9" can cause leukemia (Golub et al., 1999). Use
# appropriate Wilcoxon two-sample tests to test if, for the ALL patients, the gene
# "HOXA9 Homeo box A9" expresses at the same level as the "CD33" gene. Please
# submit your R script with the answer.

rm(list=ls())

# load the data 
data(golub, package='multtest')
gol.fac<-factor(golub.cl, level=0:1, labels=c("ALL","AML"))

# getting the row index of genes
HOX<- grep("HOXA9 Homeo box A9",golub.gnames[,2])
print("the row index of HOXA9 Homeo box A9  ")
print(HOX)
CD33<- grep("CD33",golub.gnames[,2])
print(" the row index of CD33 is")
print(CD33)


# applying the test 

wilcox.test<-wilcox.test (x= golub[1391, gol.fac=="ALL"], y= golub[808, gol.fac=="ALL"], paired=T, alternative="two.sided")
print(wilcox.test)


# problem 4 

#The data set "UCBAdmissions" in R contains admission decisions by gender at six
#departments of UC Berkeley. For this data set, carry out appropriate test for
#independence between the admission decision and gender for each of the
#departments.


rm(list=ls())

#loading the source
source("http://www.bioconductor.org/biocLite.R")
abiocLite()
library(datasets); 

str(UCBAdmissions)

Dept <- c("Dept = A","Dept = B", "Dept = C", "Dept = D", "Dept = E", "Dept = F")

for (i in 1:6 ){
  print(Dept[i])
  Dept.Data <- matrix(c(UCBAdmissions[1,1,i], UCBAdmissions[2,1,i], UCBAdmissions[1,2,i], 
                       UCBAdmissions[2,2,i]), nrow=2, dimnames=list("Admit"=c("Admitted","Rejected"), 
                                                                    "Gender"=c("Male","Female")))
# applying the chi-square test and fisher test and printing the data   
  print(Dept.Data)
  print(chisq.test(Dept.Data))
  print(fisher.test(Dept.Data))
}

# problem 5

# Please program this permutation test in R. Use this nonparametric test on the
# "CD33" gene of the Golub et al. (1999) data set. Test whether the variance in the
# ALL group is smaller than the variance in the AML group. Please submit your R
# code with the answer


rm(list=ls())

# Formula used to perform the problem.

#install the package
install.packages('gtools')
library(gtools)
#load the data
data(golub, package = "multtest")
gol.fac <- factor(golub.cl,levels=0:1, labels= c("ALL","AML"))
#find the row index
gene <- grep("gene_name",golub.gnames[,2])

data <- golub[gene,]
n <- length(data)

T.obs <- var(data[gol.fac=="ALL"]) / var(data[gol.fac=="AML"])
#observe statistic
n.perm = 2000
T.perm = NULL
for(i in 1:n.perm) {
  data.perm = sample(data, n, replace=F) 
  T.perm[i] = var(data.perm[gol.fac=="ALL"]) / var(data.perm[gol.fac=="AML"])
}

mean(T.perm <= T.obs)

# applying the above formula to CD33 gene 

# installing ackage
install.packages('gtools')
library(gtools)

# loading data
data(golub, package = "multtest")
gol.fac <- factor(golub.cl,levels=0:1, labels= c("ALL","AML"))

# finding row index
CD33 <- grep("CD33",golub.gnames[,2])

data <- golub[CD33,]
n <- length(data)
# observe statistic
T.obs <- var(data[gol.fac=="ALL"]) / var(data[gol.fac=="AML"])

n.perm = 2000
T.perm = NULL
for(i in 1:n.perm) {
  data.perm = sample(data, n, replace=F) 
  T.perm[i] = var(data.perm[gol.fac=="ALL"]) / var(data.perm[gol.fac=="AML"])
}
# p-value
mean(T.perm <= T.obs)

