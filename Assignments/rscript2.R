source("http://www.bioconductor.org/biocLite.R")
biocLite("multtest")
library(multtest)
data(golub)  

golub.gnames[1042,]
str(golub)
grep('Zyxin',golub.gnames[,2])
grep('Zyxin',golub.gnames[1042,3])
gol.fac<-factor(golub.cl,levels=0:1,labels = c("ALL","AML")) 
mean(golub[2124,gol.fac=="AML"])
labels=c("lab1","lab2","lab3","lab4")
x=sample(labels,10,replace=T)
x
table(x) 
pie(table(x))
barplot(table(x),main="Barplot")
x <- c((1:10)/10, (1:9)/3)
x
hist(x)
boxplot(golub[1042,] ~ gol.fac)
sex=sample(c("Male","Female"),10,replace=T)
smoke=sample(c("Smoke","Unsmoke"),10,replace=T)
table(sex,smoke)
plot(x,y,xlab="x-label",ylab="Y-label", main="This is the title", type="o", col="blue",pch="+",cex=2)
z <- y/2
plot(y~x,col="blue",type="o",ylab="y and z")
lines(x,z,col="red",type="o")
legend("topleft",c("this is y","this is z"),col=c("blue","red"),lty=c(1,1))



data(golub)
library(multtest)
sd.all<-apply(golub, 1,sd)
print(sum(sd.all>1))

source("http://www.bioconductor.org/biocLite.R")
biocLite("ALL")




str(trees)
plot(trees$Height~trees$Girth,col="blue",pch="+",xlim=c(0,30),ylim=c(05,25),xlab="Girth",ylab="Height and Volume")
points(trees$Girth,trees$Volume,col="red",pch="o",type="o")
legend("bottomright", c("Height","volume"), col=c("blue","red"), pch=c("+","O"))