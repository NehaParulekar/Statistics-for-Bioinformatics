#Solutions to homework module 8

# Problem 1
rm(list=ls())

#load the package
install.packages("lmtest");library(ALL);data(ALL);library(lmtest);library(Biobase);

# 1(a)Conduct the one-way ANOVA. Do the disease stages affect the mean gene
# expression value?

# patients in stages 
ALLB1234<- ALL[,ALL$BT%in%c("B","B1","B2","B3","B4")]
# exprs gives gene expression values. here we take gene 109_at values
y <- exprs(ALLB1234)["109_at",]
#anova
anova(lm(y ~ ALLB1234$BT))

# 1(b)From the linear model fits, find the mean gene expression value among B3
# patients.

summary(lm(y~ ALLB1234$BT -1))

# 1(c) Which group's mean gene expression value is different from that of group B?

pairwise.t.test(y, ALLB1234$BT)

# 1(d)Use the pairwise comparisons at FDR=0.05 to find which group means are
# different. What is your conclusion?

# FDR-adjusted pairwise
pairwise.t.test(y,ALLB1234$BT,p.adjust.method='fdr') 

# (e) Check the ANOVA model assumptions with diagnostic tests? Do we need to
# apply robust ANOVA tests here? If yes, apply the appropriate tests and state
# your conclusion.

shapiro.test(residuals(lm(y ~ ALLB1234$BT)))

bptest(lm(y ~ ALLB1234$BT), studentize = FALSE)


# Problem 2

# 2(a) Use FDR adjustments at 0.05 level. How many genes are expressed
# different in some of the groups? 

ALLB1234 <- ALL[,which(ALL$BT %in% c("B","B1","B2","B3","B4"))]
data.ALL <- exprs(ALLB1234)[,]
ALLData <- apply(data.ALL,1,function(x) kruskal.test(x ~ ALLB1234$BT)$p.value)
ALLData.fdr <- p.adjust(p=ALLData,method="fdr")
sum(ALLData.fdr<0.05)

# 2(b) Find the probe names for the top five genes with smallest p-values. 

order.ALLfdr <- order(ALLData.fdr, decreasing=FALSE)
k=1
names = NULL
for (i in order.ALLfdr[1:5]){names[k] <- names(ALLData.fdr[i])
                            k=k+1}
print('Top five genes with smallest p-values =')
print(names)


# Problem 3

# 3(a)Conduct the appropriate ANOVA analysis. Does any of the two factors
# affects the gene expression values? Are there interaction between the two
# factors?

ALL.datasex <- ALL[,which(ALL$BT %in% c("B1","B2","B3","B4") & ALL$sex %in% c("M", "F"))]
ALL.dataSex <- exprs(ALL.datasex)["38555_at",]
B.cell<-ALL.datasex$BT
ALL.sex<-ALL.datasex$sex
anova(lm(ALL.dataSex~ B.cell*ALL.sex))

# 3(b)Check the ANOVA model assumption with diagnostic tests? Are any of the
# assumptions violated?

shapiro.test(residuals(lm(ALL.dataSex ~ B.cell+ALL.sex)))

bptest(lm(ALL.dataSex ~ B.cell+ALL.sex), studentize = FALSE)


# Problem 4

# 1(a)Program this permutation test in R. 

stages <- #Example: c("B1","B2","B3")
  gene <-  #Example: "109_at"
  
Statistics <- function(stages,gene){
  
  ALL.B <- ALL[,which(ALL$BT %in% stages)]
  data <- exprs(ALL.B)[gene,]
  group <- ALL.B$BT[,drop=T]
  
  g <- length(stages)
  Means <- summary(lm(data ~ group-1))[["coefficients"]][1:g]
  Total.Mean <- (1/g)*sum(Means)
  
  MUj_MU <- NULL
  for (i in 1:g){
    
    MUj_MU[i] <- (Means[i]-Total.Mean)^2
  } 
  T.obs <- (1/(g-1))*sum(MUj_MU) #Observed statistic
  
  n <- length(data)
  n.perm = 2000
  T.perm = NULL
  for(i in 1:n.perm) {
    data.perm = sample(data, n, replace=F)
    Means.Perm <- summary(lm(data.perm ~ ALL.B$BT-1))[["coefficients"]][1:g]
    Total.MeanPerm <- (1/g)*sum(Means.Perm)
    MUj_MU1 <- NULL
    for (k in 1:g){
      MUj_MU1[k] <- (Means.Perm[k]-Total.MeanPerm)^2
    } 
    T.perm[i] = (1/(g-1))*sum(MUj_MU1) #Permuted statistic
  }
  mean(T.perm>=T.obs) #p-value
}

# 4(b)Run this permutation test on the Ets2 repressor gene 1242_at on the patients
# in stage B1, B2, and B3 from the ALL data set.

stages <- c("B1","B2","B3")
gene <- "1242_at"

Statistics(stages,gene)
# function is defined in 4a

