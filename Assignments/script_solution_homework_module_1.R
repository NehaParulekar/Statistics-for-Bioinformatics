#Solutions to Module 1 Homework

# Question 1

#(a) What is the class of the object defined be vec <-c(5,TRUE) ?
rm(list = ls())

vec<-c(5,TRUE)
class(vec)

# (b) Suppose I have vectors x <- 1:4 and y <- 1:2. What is the result of the
# expression x + y?

x<-1:4
y<-1:2
z<-(x+y)
z 

#(c) Suppose I define the following function in R:
#fsin<-function(x) sin(pi*x)
#What will be returned by fsin(1) ?

fsin<-function(x)sin(pi*x)
fsin(1)



#(d) What is returned by the R command c(1,2) %*% t(c(1,2)) ?

c(1,2)%*%t(c(1,2))

#(e) Suppose I define the following function in R:
f <- function(x) {
  g <- function(y) {
    y + z
  }
  z <- 4
  x + g(x)
}
# If I then run in R the following statements
z <- 15
f(3)
# What value is returned?
rm(list = ls())

f <- function(x) {
  g <- function(y) {
    y + z
  }
  z <- 4
  x + g(x)
}
z <- 15
f(3)

# Question 2

# Use R to calculate the given equation :

x=1:1000
y<-sum(x^2)
y

# Question 3

#Write an R script that does all of the following:

# a) Create a vector X of length 20, with the kth element in X = 2k, for
# k=1…20. Print out the values of X.

rm(list = ls())

v=(1:20)
k<-v
X<-(2*k)
print(X)

# b) Create a vector Y of length 20, with all elements in Y equal to 0. Print
# out the values of Y.

Y<-rep(0,20)
print (Y)

#c) Using a for loop, reassigns the value of the k-th element in Y, for k =
# 1…20. When k < 12, the kth element of Y is reassigned as the cosine
# of k. When the k = 12, the kth element of Y is reassigned as the value
# integrate(function(x) sqrt(t), lower=0, upper=k) 


integrand<-function(t)sqrt(t)
for(K in 1:20){
  if (k < 12) {
    Y[k]<-cos(k)
  }
  else  {
    Y[k]<-integrate(integrand,lower=0, upper=k)
  }
}
print(Y)

