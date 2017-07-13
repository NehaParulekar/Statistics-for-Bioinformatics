# rough self assesment

f_x<-function (x)(x^2)/338350
sum(f_x(10:49))

X_range<-c(1:100)
f.x<-function(x) (x^2)/338350
f_x<-function(x)f.x(x)*(x%in%X_range)
EX<-sum(X_range*f_x(X_range))
VarX<- sum((X_range-EX)^2*f_x(X_range))
VarX
stdX<-sqrt(VarX)
stdX


rm(list=ls())
X<-rchisq( n=100000, df=3)
mean((1<X)&(X<4))

integrate (function(x)rnorm(n=500000, mean=1.6, sd=0.4), lower=1, upper=1.6) 




x<-rnorm(n=500000, mean=1.6, sd=0.4)
mean((x<1.6)&(x<1))




x<-rnorm(n=1000)
x

x<-rf(100000,)

mean((0.4<x) & (x<1.5))