# Solutions to Module 3 Homework.

rm(list=ls())
# Problem 1

integrate (function(x) 2*exp(-2*(x-1)), lower=1, upper=4)

# Problem 2
# find P(X=1)

2*exp(-2)

# Find P(-2<X<4)
f_x<- function(x) (2^x/factorial(x))*exp(-2)
sum(f_x(0:3))

rm(list=ls())

# Problem 4

# For Y following a binomial (n = 3, p = 0.25) distribution, compute the
# following:

# P(Y ??? 2) =
 set.x<- c(0,1,2)
sum(dbinom(set.x, size=3, p=0.25))

# E(Y) = np = 3*0.25

3*0.25

# Var(Y) = np(1-p)

3*0.25*(1-0.25)

rm(list=ls())

# Problem 5

# For X following a Chi-square distribution with degree of freedom m = 3,
# compute the following:
#  P(1 < X < 4) = 
#  Method 1: cumulative Distribution Function
pchisq(4,3)-pchisq(1,3)
# Method 2: Integrating the pdf
integrate( function(x) dchisq(x,3), lower=1, upper=4)

# b) use a Monte Carlo simulation with sample size n=100,000 to estimate
# P(1 < X < 4).
X<-rchisq( n=100000, df=3)
mean((1<X)&(X<4))


rm(list=ls())

# Problem 7
# a) What is the probability that a randomly chosen patient have the Zyxin
# gene expression values between 1 and 1.6?
rm(list=ls())
integrate (function(x)dnorm(x, mean=1.6, sd=0.4), lower=1, upper=1.6)$value
# or 
pnorm(1.6, mean=1.6,sd=0.4)- pnorm(1, mean=1.6, sd=0.4)

# (b) Use a Monte Carlo simulation of sample size n=500,000 to estimate the
# probability in part (a). Give your R code, and show the value of your
# estimate

x<-rnorm(n=500000, mean=1.6, sd=0.4)
mean((x<1.6)&(1<x))

# (c) What is the probability that exactly 2 out of 5 patients have the Zyxin
# gene expression values between 1 and 1.6? 

dbinom(2, size=5, prob=0.4331928)


rm(list=ls())

# Problem 8

# (a) Hand in a R script that calculates the mean and variance of two random
# variables X~F(m=2,n=5) and Y~F(m=10,n=5) from their density functions. 

EX<- integrate(function(x) x*df(x, df1=2, df2=5), lower=0, upper=Inf)$value
EX
EY<- integrate(function(y) y*df(y, df1=10, df2=5), lower=0, upper=Inf)$value
EY
VarX<- integrate(function(x) (x-EX)^2*df(x, df1=2, df2=5), lower=0, upper=Inf)$value
VarX
VarY<- integrate(function(y) (y-EY)^2*df(y, df1=10, df2=5), lower=0, upper=Inf)$value
VarY

# (b) Use the formula in Table 3.4.1 to calculate the means and variances
# directly. 

# mean= n/(n-2)
mean<-5/(5-2)
mean
#variance= (2*n^2*(m+n-2))/(m*(n-2)^2*(n-4))
variance<- function(m,n) (2*n^2*(m+n-2))/(m*(n-2)^2*(n-4))
variance(2,5)
variance(10,5)
