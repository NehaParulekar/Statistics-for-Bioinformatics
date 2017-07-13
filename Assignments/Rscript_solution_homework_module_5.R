# Solutions to module 5

# Problem 1
rm(list=ls())
# 1(b) A random sample of size 6 from the exp(??) distribution results in observations:
# 1.636, 0.374, 0.534, 3.015, 0.932, 0.179. Find the MLE on this data set in two
# ways: by numerical optimization of the likelihood and by the analytic formula.

# ANALYTIC MLE FORMULA
# to get mean 
(1.636 + 0.374 + 0.534 + 3.015 + 0.932 + 0.179)/6
# solving for lanbda 
1/1.111667

# NUMERICAL MLE FORMULA
nmle<- function(x) - sum(log(dexp(c(1.636, 0.374, 0.534, 3.015, 0.932, 0.179),x)))
nmle.results<-optim(1, nmle)
print(nmle.results$par)


# Problem 2

# 2(b) Finda one-sided 90% lower confidence interval of m. 

# mean = 100.8, sd= 12.4
LC<- 100.8 + (qt(0.1, df=53-1)*(12.4/sqrt(53)))
LC


# Problem 3

rm(list=ls())
# Using the Golub et al. (1999) data set, analyze the Zyxin gene expression
# data separately for the ALL and AML groups.

# 3(a)(a)Find the bootstrap 95% CIs for the mean and for the variance of the gene
# expression in each group separately

# bring golub daa set
data(golub, package="multtest")
gol.fac<-factor(golub.cl, levels=0:1, labels=c("ALL","AML"))
#finding out the Zyxin gene data row
Zrow<-grep("Zyxin",golub.gnames[,2])
Zrow


# ALL :Bootstrap 95% CIs for mean and variance 
ZALL<- golub[Zrow,gol.fac=="ALL"]
nALL<- length(ZALL)
nboot<-1000
ALL.xbar<- rep(NA, nboot)
ALL.xvar<- rep(NA, nboot)
for (i in 1:nboot){
  ALLdata.star<- ZALL[sample(1:nALL, replace=TRUE)]
   ALL.xbar[i]<-mean(ALLdata.star)
   ALL.xvar[i]<-var(ALLdata.star)
}
CI.ALL.mean<-quantile(ALL.xbar,c(0.025,0.975))
CI.ALL.var<-quantile(ALL.xvar,c(0.025,0.975))
print("Mean expression of Zyxin for ALL group")
print(mean(ZALL))
print("95% Bootstrap CI for ALL group Zyxin mean expression")
print(CI.ALL.mean)
print("variance expression of zyxin for ALL group")
print(var(ZALL))
print("95% Bootstrap CI for ALL group zyxin variance")
print(CI.ALL.var)




# bring golub daa set
data(golub, package="multtest")
gol.fac<-factor(golub.cl, levels=0:1, labels=c("ALL","AML"))
#finding out the Zyxin gene data row
Zrow<-grep("Zyxin",golub.gnames[,2])
Zrow

# AML :Bootstrap 95% CIs for mean and variance 
ZAML<- golub[Zrow,gol.fac=="AML"]
nAML<- length(ZAML)
nboot<-1000
AML.xbar<- rep(NA, nboot)
AML.xvar<- rep(NA, nboot)
for (i in 1:nboot){
  AMLdata.star<- ZAML[sample(1:nAML, replace=TRUE)]
  AML.xbar[i]<-mean(AMLdata.star)
  AML.xvar[i]<-var(AMLdata.star)
}
CI.AML.mean<-quantile(AML.xbar,c(0.025,0.975))
CI.AML.var<-quantile(AML.xvar,c(0.025,0.975))
print("Mean expression of Zyxin for AML group")
print(mean(ZAML))
print("95% Bootstrap CI for AML group Zyxin mean expression")
print(CI.AML.mean)
print("variance expression of zyxin for AML group")
print(var(ZAML))
print("95% Bootstrap CI for AML group zyxin variance")
print(CI.AML.var)


# 3(b) Find the parametric 95% CIs for the mean and for the variance of the gene
# expression in each group separately. (You need to choose the appropriate
#  approximate formula to use: z-interval, t-interval or chi-square interval.)

# t-interval for ALL Zyxin mean and chi-square for ALL zyxin varianve expression
Zrow<-grep("Zyxin",golub.gnames[,2])
ZALL<- golub[Zrow,gol.fac=="ALL"]
nALL<- length(ZALL)
ci.mean.ALL<-mean(ZALL)+qt(c(0.025,0.975), df=nALL-1)*sd(ZALL)/sqrt(nALL)
print("95% CI's (t-interval) for the mean for All" )
print(ci.mean.ALL)
ci.var.ALL<-((nALL-1)*var(ZALL))/qchisq(c(0.975,0.025), df=nALL-1)
print("95% CI's (chi-square) for variance for ALL")
print(ci.var.ALL)


# t-interval for AML Zyxin mean and chi-square for AML zyxin varianve expression
Zrow<-grep("Zyxin",golub.gnames[,2])
ZAML<- golub[Zrow,gol.fac=="AML"]
nAML<- length(ZAML)
ci.mean.AML<-mean(ZAML)+qt(c(0.025,0.975), df=nAML-1)*sd(ZAML)/sqrt(nAML)
print("95% CI's (t-interval) for the mean for AML" )
print(ci.mean.AML)
ci.var.AML<-((nAML-1)*var(ZAML))/qchisq(c(0.975,0.025), df=nAML-1)
print("95% CI's (chi-square) for variance for AML")
print(ci.var.AML)



# 3(c) Find the bootstrap 95% CI for the median gene expression in both groups
# separately.

# Bootstrap for 95% CI for median for ALL 

data(golub, package="multtest")
gol.fac<-factor(golub.cl, levels=0:1, labels=c("ALL","AML"))
Zrow<-grep("Zyxin",golub.gnames[,2])
Zrow
ZALL<- golub[Zrow,gol.fac=="ALL"]
nALL<- length(ZALL)
nboot<-1000
ALL.median<- rep(NA, nboot)
for (i in 1:nboot){
  ALLdata.star<- ZALL[sample(1:nALL, replace=TRUE)]
  ALL.median[i]<-median(ALLdata.star)
}
CI.ALL.median<- quantile(ALL.median, c(0.025, 0.975))
print("Median expression of zyxin for ALL")
print(CI.ALL.median)

# Bootstrap for 95% CI for median for AML 

data(golub, package="multtest")
gol.fac<-factor(golub.cl, levels=0:1, labels=c("ALL","AML"))
Zrow<-grep("Zyxin",golub.gnames[,2])
Zrow
ZAML<- golub[Zrow,gol.fac=="AML"]
nAML<- length(ZAML)
nboot<-1000
AML.median<- rep(NA, nboot)
for (i in 1:nboot){
  AMLdata.star<- ZALL[sample(1:nAML, replace=TRUE)]
  AML.median[i]<-median(AMLdata.star)
}
CI.AML.median<- quantile(AML.median, c(0.025, 0.975))
print("Median expression of zyxin for AML")
print(CI.AML.median)

# Problem 4

rm(list=ls())
# 4(a) Write a R-script to conduct a Monte Carlo study for the coverage probabilities
# of the two CIs. That is, to generate nsim=1000 such data sets from the Poisson
# distribution. Check the proportion of the CIs that contains the true parameter ??. 

# finding no of simulations and mean for formula 1

#Number of simulations and generate dataset
nsim <- 1000
lambda <-
  
getdata <- matrix(rpois(50*nsim,10),nrow=nsim)
new.lambda <- (apply(getdata,1,mean))
tdist <- qt(.05,49) * sqrt(new.lambda/50)

#90% CI for sample mean
method1.low = new.lambda+tdist
method1.High = new.lambda-tdist

#Check coverage probabilities
sum(method1.low<lambda & lambda<method1.High)/1000

# finding no of simulations and mean for formula 2

#Number of simulations and generate dataset
nsim <- 1000
lambda <-
getdata <- matrix(rpois(50*nsim,10),nrow=nsim)
new.lambda <- (apply(getdata,1,mean))

#90% CI for sample mean
method2.low = 49 *(new.lambda)/qchisq(.95,49)
method2.High = 49 *(new.lambda)/qchisq(.05,49)

#Check coverage probabilities
coverage<-sum(method2.low<lambda & lambda<method2.High)/1000
coverage

# 4(b) (b) Run the Monte Carlo simulation for nsim=1000 runs, at three different
# parameter values: ??=0.1, ??=1 and ??=10. Report the coverage probabilities of these
# two CIs at each of the three parameter values.

# finding no of simulations and mean for formula 1 when ??=0.1 

#Number of simulations and generate dataset

nsim<-1000
lambda<-0.1
getdata<-matrix(rpois(50*nsim,lambda),nrow=nsim)
new.lambda<-(apply(getdata,1,mean))
tdist<- qt(.05,49) * sqrt(new.lambda/50)

#90% CI for sample mean
method1.low=new.lambda+tdist
method1.High=new.lambda-tdist

#Check if calculated CI contains lambda
coverage1<-sum(method1.low<lambda & lambda<method1.High)/1000

# finding no of simulations and mean for formula 2 when ??=0.1

#Number of simulations and generate dataset
nsim <- 1000
lambda <- 0.1
getdata <- matrix(rpois(50*nsim,lambda),nrow=nsim)
new.lambda <- (apply(getdata,1,mean))

#90% CI for sample mean
method2.low = 49 *(new.lambda)/qchisq(.95,49)
method2.High = 49 *(new.lambda)/qchisq(.05,49)

#Check if calculated CI contains lambda
coverage2<-sum(method2.low<lambda & lambda<method2.High)/1000

# finding no of simulations and mean for formula 1 when ??=1

#Number of simulations and generate dataset
nsim<-1000
lambda<-1
getdata<-matrix(rpois(50*nsim,lambda),nrow=nsim)
new.lambda<-(apply(getdata,1,mean))
tdist<- qt(.05,49) * sqrt(new.lambda/50)

#90% CI for sample mean
method1.low=new.lambda+tdist
method1.High=new.lambda-tdist

#Check if calculated CI contains lambda
coverage3<-sum(method1.low<lambda & lambda<method1.High)/1000

# finding no of simulations and mean for formula 2 when ??=1

#Number of simulations and generate dataset
nsim <- 1000
lambda <- 1
getdata <- matrix(rpois(50*nsim,lambda),nrow=nsim)
new.lambda <- (apply(getdata,1,mean))

#90% CI for sample mean
method2.low = 49 *(new.lambda)/qchisq(.95,49)
method2.High = 49 *(new.lambda)/qchisq(.05,49)

#Check if calculated CI contains lambda
coverage4<-sum(method2.low<lambda & lambda<method2.High)/1000

# finding no of simulations and mean for formula 1 when ??=10

#Number of simulations and generate dataset
nsim<-1000
lambda<-10
getdata<-matrix(rpois(50*nsim,lambda),nrow=nsim)
new.lambda<-(apply(getdata,1,mean))
tdist<- qt(.05,49) * sqrt(new.lambda/50)

#90% CI for sample mean
method1.low=new.lambda+tdist
method1.High=new.lambda-tdist

#Check if calculated CI contains lambda
coverage5<-sum(method2.low<lambda & lambda<method2.High)/1000

# finding no of simulations and mean for formula 2 when ??=10

#Number of simulations and generate data
nsim <- 1000
lambda <- 10
getdata <- matrix(rpois(50*nsim,lambda),nrow=nsim)
new.lambda <- (apply(getdata,1,mean))

#90% CI for sample mean
method2.low = 49 *(new.lambda)/qchisq(.95,49)
method2.High = 49 *(new.lambda)/qchisq(.05,49)

#Check if calculated CI contains lambda
coverage6<-sum(method2.low<lambda & lambda<method2.High)/1000

print("For nsim = 1000 & lambda = 0.1 coverage prabability of lambda using method 1 is")
print (coverage1)
print("For nsim = 1000 & lambda = 0.1 coverage prabability of lambda using method 2 is")
print (coverage2)
print("For nsim = 1000 & lambda = 1 coverage prabability of lambda using method 1 is")
print (coverage3)
print("For nsim = 1000 & lambda = 1 coverage prabability of lambda using method 2 is")
print (coverage4)
print("For nsim = 1000 & lambda = 10 cverage prabability of lambda using method 1 is")
print (coverage5)
print("For nsim = 1000 & lambda = 10 cverage prabability of lambda using method 2 is")
print (coverage6)

