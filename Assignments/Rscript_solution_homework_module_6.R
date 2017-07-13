# solution to Module 6

# Problem 1

rm(list=ls())

# Setting Golub data set for Problem 1a-1f
# laoding the golub data set 
data(golub, package="multtest")
H4j<-grep("H4/j",golub.gnames[,2])
APS<-grep("APS Prostate specific antigen", golub.gnames[,2])
# creating factor
gol.fac<- factor(golub.cl, levels=0:1, labels=c("ALL","AML"))

# 1(a) The mean "H4/j gene" gene expression value in the ALL group is greater than -1.
print(t.test(golub[H4j, gol.fac=="ALL"], mu= -1, alternative="greater"))

# 1(b) The mean "H4/j gene" gene expression value in ALL group differs from the mean "H4/j gene" gene expression value in the AML group.
print(t.test(golub[H4j, gol.fac=="ALL"], golub[H4j, gol.fac=="AML"]))

# 1(c) In the ALL group, the mean' expression value for the "H4/j gene" gene is lower than the mean expression value for the "APS Prostate specific antigen" gene.
print(t.test(golub[H4j, gol.fac=="ALL"], golub[APS, gol.fac=="ALL"], alternative="less", paired=T))

# 1(d) Let plow denote the proportion of patients for whom the "H4/j gene" expression is lower than the "APS Prostate specific antigen" expression. We wish to show that plow in the ALL group is greater than half. Does this test conclusion agree with the conclusion in part (c)?
Plow<-sum(golub[H4j, gol.fac=="ALL"] < golub[APS, gol.fac=="ALL"])
H4jlength<-length(golub[H4j, gol.fac=="ALL"])
print(binom.test(x=Plow, n=H4jlength, p=0.5, alternative="greater"))

# 1(e) Let pH4j denotes the proportion of patients for whom the "H4/j gene" expression values is greater than -0.5. We wish to show that pH4j in the ALL group is less than 0.5.
PH4j<-sum(golub[H4j, gol.fac=="ALL"] > -0.5)
PH4jlen<-length(golub[H4j, gol.fac=="ALL"])
print(binom.test(x=PH4j, n=PH4jlen, p=0.5, alternative="less"))

# 1(f) pH4j in the ALL group differs from pH4j in the AML group.
PH4j.ALL<-sum(golub[H4j, gol.fac=="ALL"] > -0.5)
PH4j.AML<-sum(golub[H4j, gol.fac=="AML"] > -0.5)
ALL.length<- length(golub[H4j, gol.fac=="ALL"])
AML.length<- length(golub[H4j, gol.fac=="AML"])
print(prop.test(x= c(PH4j.ALL,PH4j.AML), n=c(ALL.length, AML.length), alternative="two.sided"))


# problem 2(b)

rm(list=ls())

# using dbinom to find the probability of fewer than 20 rejects.
sum(dbinom(x=0:19, size= 1000, prob =0.05))
# using pbinom to find the probability of fewer than 20 rejects.
pbinom(q=19, size=1000, prob=0.05)

# problem 3(a)

rm(list=ls())

#creating the data set
x.simul<- matrix(rnorm(10000*20, mean=3, sd=4), ncol=20)
# t-test 
tstat<- function(x) (mean(x)-3)/(sd(x)/sqrt(length(x)))
tstat.simul<-apply(x.simul,1,tstat)
#calculating the rejection rate
power.simul<- mean(tstat.simul > qt(0.3, df=19) & tstat.simul < qt(0.4, df=19) )
# type I error rate with its 95% CI
print(power.simul+ c(-1,0,1)* qnorm(0.975)*sqrt(power.simul*(1-power.simul)/10000))


# problem 4

rm(list=ls())

# 4(a) Use Bonferroni and FDR adjustments both at 0.05 level. How many genes are differentially expressed according to these two criteria?

# load the golub data set and create factor
data(golub, package="multtest")
gol.fac<- factor(golub.cl, levels=0:1, labels=c("ALL","AML"))

# to get the no. of genes and apply welch two-sample test 

length<- length(golub.gnames[,2])
pvalues <- NULL
for (i in 1:length){
  pvalue <- t.test(golub[i,gol.fac=="ALL"], golub[i,gol.fac=="AML"])$p.value
  pvalues <- c(pvalues, pvalue)
}

# performing Bonferroni and FDR adjustment
p.bon<-p.adjust(p=pvalues, method="bonferroni")
p.fdr<-p.adjust(p=pvalues, method="fdr")

# printing answer
print("total number of genes differentially expressed at 0.05 level  no adjustment")
sum(pvalues < 0.05 )
print("total number of genes differentially expressed at 0.05 level  with bonferroni ")
sum(p.bon<0.05)
print("total number of genes differentially expressed at 0.05 level  with FDR")
sum(p.fdr<0.05)

# 4(b) Find the gene names for the top three strongest differentially expressed genes (i.e., minimum p-values). Hint: the gene names are stored in golub.gnames

# load the golub data set and create factor
data(golub, package="multtest")
gol.fac<- factor(golub.cl, levels=0:1, labels=c("ALL","AML"))

# to get the no. of genes and apply welch two-sample test 

length<- length(golub.gnames[,2])
pvalues <- NULL
for (i in 1:length){
  pvalue <- t.test(golub[i,gol.fac=="ALL"], golub[i,gol.fac=="AML"])$p.value
  pvalues <- c(pvalues, pvalue)
}

# performing Bonferroni and FDR adjustment
p.bon<-p.adjust(p=pvalues, method="bonferroni")
p.fdr<-p.adjust(p=pvalues, method="fdr")

#printing results

print("top three strongest differentially espressed genes for FDR")
p.fdr<-p.adjust(p=pvalues,method="fdr")
orderAML<-order(p.fdr, decreasing=FALSE)
golub.gnames[orderAML[1:3],2]

print("top three strongest differentially espressed genes for Bonferroni")
p.bon<-p.adjust(p=pvalues,method="bon")
orderAL<-order(p.bon, decreasing=FALSE)
golub.gnames[orderAML[1:3],2]


# problem 5 

rm(list=ls())

# 5 (a) Program R functions to calculate the Wald CI, the Wilson CI and the Agresti-Coull CI for binomial proportion. 
 
wald.CI<- function(X,n,alpha=0.05){
  p<- X/n
  z<- qnorm(1-alpha/2)
  c(c(p,(p + c(-1,1) * z *sqrt((p*(1-p))/n))))
}


wilson.CI<- function(X,n,alpha=0.05){
  p<- X/n
   z<- qnorm(1-alpha/2)
  c(p,((1/(1+z^2/n))* (p+(z^2/(2*n))+ c(-1,1)*z*sqrt((p*(1-p))/n+z^2/(4*n^2)))))
}


AgC.CI <- function(X,n,alpha=0.05){
  z<- qnorm(1-alpha/2)
  N<- n + z^2
  p<- (X + z^2/2)/N
   c(p,(p + c(-1,1)*z*sqrt((p*(1-p)))))
}


# 5(b) Run a Monte Carlo simulation to check the coverage of the Wald CI, the Wilson CI and the Agresti-Coull CI for n=40 and p=0.2 at the nominal confidence level of 95%. Do 10,000 simulation runs for calculating the empirical coverages

n.40<- rbinom(n=1, size=40, p=0.2)
n.40.wald<-wald.CI(n.40,40)
n.40.wilson<-wilson.CI(n.40,40)
n.40.AgC<-AgC.CI(n.40,40)

print(" the proportion of success for n=40 & p=0.2 ")
print(n.40)
print("The 95% CI for n=40 & p=0.2 ")
print(" for Wald CI")
print(rbind(n.40.wald))
print(" for wilson CI")
print(rbind(n.40.wilson))
print(" for AgC CI")
print(rbind(n.40.AgC))


# run 10000 simulations to calculate the emprical changes.

simul<- rbinom(n=10000, size=40, p=0.2)
simul.wald<- NULL
simul.wilson<- NULL
simul.AgC<- NULL
for(i in simul){
  simul.wald<- rbind(simul.wald, wald.CI(i,40))
  simul.wilson<- rbind(simul.wilson, wilson.CI(i,40))
  simul.AgC<- rbind(simul.AgC, AgC.CI(i,40))
}

# calculating the coverage of CI intervals and printing the results

print("The estimated coverage after 10000 simulations of n=40, p=0.2")
print("for wald CI coverage")
wald.coverage<- mean(0.2 > simul.wald[,2] & 0.2 < simul.wald[,3])
print(wald.coverage)
print("for wilson CI coverage")
wilson.coverage<- mean(0.2 > simul.wilson[,2] & 0.2 < simul.wilson[,3])
print(wilson.coverage)
print("for AgC CI coverage")
AgC.coverage<- mean(0.2 > simul.AgC[,2] & 0.2 < simul.AgC[,3])
print(AgC.coverage)
