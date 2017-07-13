# solutions to homework module 9

# Problem 1 

rm(list=ls())
#load the data from multtest package

data(golub,package="multtest") 
grep("GRO2 GRO2",golub.gnames[,2])
grep("GRO3 GRO3",golub.gnames[,2])
GRO2.data <- golub[2714,]
GRO3.data <- golub[2715,]

# 1(a)Find the correlation between the expression values of these two genes.

# apply Correlation test to check if true correlation = 0

cor.test(GRO2.data,GRO3.data)

# 1(b)Find the parametric 90% confident interval for the correlation with cor.test(). 

cor.test(GRO2.data, GRO3.data, conf.level = 0.90)

# 1(c) Find the bootstrap 90% confident interval for the correlation

# resampling 2000 times
nboot <- 2000
# giving a vector to save the resampled statistics
boot.cor <- matrix(0, nrow=nboot, ncol = 1) 
#data set with GRO2 and GRO3 in two columns 
data <- cbind(GRO2.data, GRO3.data)
for (i in 1:nboot){
# Resample the pairs
  dat.star <- data[sample(1:nrow(data),replace=TRUE),] 
# correlation on the resampled data 
  boot.cor[i,] <- cor(dat.star[,1], dat.star[,2])
}
# find the quantiles for resampled statistics
quantile(boot.cor[,1],c(0.05,0.95))

# 1(d) Test the null hypothesis that correlation = 0.64 against the one-sided
# alternative that correlation > 0.64 at the ?? = 0.05 level. What is your
# conclusion? Explain you reasoning supported by the appropriate R outputs.

# resampling 2000 times
nboot <- 2000
# giving a vector to save the resampled statistics
boot.cor <- matrix(0, nrow=nboot, ncol = 1) 
#data set with GRO2 and GRO3 in two columns 
data <- cbind(GRO2.data, GRO3.data)
for (i in 1:nboot){
  # Resample the pairs
  dat.star <- data[sample(1:nrow(data),replace=TRUE),] 
  # correlation on the resampled data 
  boot.cor[i,] <- cor(dat.star[,1], dat.star[,2])
}
# find the quantiles for resampled statistics
quantile(boot.cor[,1],c(0.025,0.975))


# Problem 2

rm(list=ls())
#load the data from multtest package
data(golub,package="multtest") 
grep("Zyxin",golub.gnames[,2])

# 2(a)How many of the genes have correlation values less than negative 0.5? 

Zyxin.data <- golub[2124,]
cor.data<-apply(golub, 1, function(x) cor.test(x, Zyxin.data)$estimate)
# to get  no of genes that have correlation values less than negative 0.5
print("no of genes that have correlation values less than negative 0.5")
print(sum(cor.data < -0.5))

# 2(b)Find the gene names for the top five genes that are most negatively
# correlated with Zyxin gene.

order.cor <- order(cor.data, decreasing=FALSE)
print("gene names for the top five genes that are most negatively correlated with Zyxin gene.")
golub.gnames[order.cor[1:5],2]

# 2(c) Using the t-test, how many genes are negatively correlated with the Zyxin
# gene? Use a false discovery rate of 0.05.

# using t-test
cor.ttest <- apply(golub, 1, function(x) cor.test(x, Zyxin.data, alternative = "l")$p.value)
print("genes are negatively correlated with the Zyxin gene.")
sum(cor.ttest < 0.05)
# after adjusting for FDR
cor.fdr <- p.adjust(p = cor.ttest, method = "fdr")
print("After FDR adjustment")
sum(cor.fdr < 0.05)


# Problem 3 

rm(list=ls())
#load the data from multtest package

data(golub,package="multtest") 
grep("GRO2 GRO2",golub.gnames[,2])
grep("GRO3 GRO3",golub.gnames[,2])
GRO2.data <- golub[2714,]
GRO3.data <- golub[2715,]

# 3(a)Is there a statistically significant linear relationship between the two genes'
# expression? Use appropriate statistical analysis to make the conclusion.
# What proportion of the GRO3 GRO3 oncogene expression's variation can
# be explained by the regression on GRO2 GRO2 oncogene expression?

reg.fit <- lm(GRO3.data ~ GRO2.data)
reg.fit  
summary(reg.fit)

# 3(b)Test if the slope parameter is less than 0.5 at the ?? = 0.05 level

#Show 95% 2-sided CIs from regression fit
confint(reg.fit, level = 0.95)

# 3(c) Find an 80% prediction interval for the GRO3 GRO3 oncogene expression
# when GRO2 GRO2 oncogene is not expressed (zero expression value).

reg.fit <- lm(GRO3.data ~ GRO2.data)
newdata <- data.frame(GRO2.data = 0)
predict(reg.fit, newdata, interval="prediction", level = 0.80)

# 3(d)Check the regression model assumptions. Can we trust the statistical
# inferences from the regression fit? 

# For Normal Assumption we plot/draw scatter plot
plot(reg.fit, which = 2)


# For Confirmation we apply Shapiro test to reg.fit and draw the plot 
shapiro.test(resid(reg.fit))

plot(reg.fit,which = 1)


# Problem 4 

# 4(a)Regress stack.loss on the other three variables. What is the fitted regression
# equation?

data(stackloss)
str(stackloss)

stack.loss<- as.data.frame(stackloss[,c('Air.Flow', 'Water.Temp', 'Acid.Conc.', 'stack.loss')])
lin.reg <- lm(stack.loss~Air.Flow+Water.Temp+Acid.Conc., data=stack.loss)
summary(lin.reg)

#(c) Find a 90% confidence interval and 90% prediction interval for stack.loss
# when Air.Flow=60, Water.Temp=20 and Acid.Conc.=90.

given.data <- data.frame(Air.Flow = 60, Water.Temp = 20, Acid.Conc. = 90)
predict(lin.reg, given.data, interval="confidence", level = 0.90)
predict(lin.reg, given.data, interval="prediction", level = 0.90)
