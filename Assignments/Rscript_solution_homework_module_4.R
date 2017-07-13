# Solutions to Module 4 

# Problem 2
rm(list=ls())
1-pnorm(15, mean=14,sd=2.04393/sqrt(100)) 

# Problem 3 
rm(list=ls())
require(mvtnorm)
data<-rmvnorm(50,mean=c(9,10),sigma=matrix(c(3,2,2,5),nrow=2))
meanxy<-apply(data,1,mean)
varxy<-apply(data,1,var)
mean(meanxy) + c(-1,1)*1.96*sqrt(var(meanxy)/50) #95%CI for mean
mean(varxy) + c(-1,1)*1.96*sqrt(var(varxy)/50) #95%CI for variance
sqrt(1.912846)
sqrt(4.173590)
rnorm(10000, mean=c(9.27048, 10.40593), sd=c(1.383057, 4.173590))


# Problem 4 

rm(list=ls())

#chi-square distribution
X1 <- rchisq(10000, df=10)
# gamma distribution
X2 <- rgamma(10000, shape=1, scale=2)
# t distribution
X3 <- rt(10000, df=3)
# calculate mean E(Y)
Y <- sqrt(X1)*(X2)+4*(X3^2)
meanY<- mean(Y)
# Print mean E(Y)
print (meanY)



# Problem 5

#take a sample of 1000 from standard normal distribuion
# Repeating it 1000 times
rm(list=ls())
size<- 1000
nsim<- 1000
my.data <- matrix(rnorm(size*nsim, mean=0, sd=1))
maxima <- apply(my.data,1, max)
#subtract an from the maxima value and then divide it by bn
n<- size
an <- sqrt(2*log(n)) - 0.5*(log(log(n))+log(4*pi))*(2*log(n))^(-1/2)
bn <- (2*log(n))^(-1/2)
plot_data<-(maxima - an/bn)
# now we plot the density from the normalised maxima
hist(plot_data, freq=F,xlim=c(-5,5), ylim=c(0,0.5))
# now add the density function 
f<-function(x){exp(-x)*exp(-exp(-x))}
curve(f, add=TRUE, col = "blue")
curve(dnorm, add=TRUE, col = "red")
legend("topright",c("maxima","f","dnorm"),col=c("yellow","red","blue"),lty=c(1,1,1))
