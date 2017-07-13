# solutions to Homework module 10 

# problem 1

# Loading the yeast microarray data 

source("http://bioconductor.org/biocLite.R")
biocLite("annotate")
library("annotate")
biocLite("ArrayExpress")
library(ArrayExpress)
biocLite("expresso")
library(expresso)
yeast.raw <- ArrayExpress('E-MEXP-1551')
a
# 1(a) Preprocess the raw data set into an expression data set using: the "mas"
# background correction method, the "quantiles" normalization method, "pmonly"
# pm correction method and "medianpolish" summary method. Give the R command
# here for doing this task
yeast.raw <- ArrayExpress('E-MEXP-1551')
eset<- expresso(yeast.raw,bgcorrect.method="mas",
                    normalize.method="quantiles",
                    pmcorrect.method="pmonly",
                    summary.method="medianpolish")
exprs.yeast <- exprs(eset)

# 1(b) Print out the mean expression values for the first five genes across all samples.

apply(exprs.yeast[1:5,], 1, mean)

# 1(c) How many genes and how many samples are in the preprocessed expression
# data set?

str(exprs.yeast)


# Problem 2 

# load the data
biocLite("ArrayExpress") 
biocLite("annotate");
library(ArrayExpress)
yeast.raw <- ArrayExpress('E-MEXP-1551')

# 2a) What is the annotation package for the yeast data set in question 1? 

annotation(yeast.raw)

# 2b) Search the 1769308_at gene GO numbers related to Molecular Function (MF).
# How many GO numbers do you get?
biocLite("yeast2.db")
library(yeast2.db)
go1769308_at <- get("1769308_at", env = yeast2GO)
library(annotate)
mf.go1769308_at <- getOntology(go1769308_at,"MF")
length(mf.go1769308_at)

# 2(c) Find the GO parents of the GO IDs in part (b). How many GO parents are there?

biocLite("GO.db")
library(GO.db)

x <- get("1769308_at", env = yeast2GO)
gonr <- getOntology(go1769308_at, "MF")
#parents for GO numbers
go.P <- getGOParents(gonr)
#extract GO numbers for parents
parents <- sapply(go.P,function(x) x$Parents) 
parents
length(parents)

# 2(d) Find the GO children of the GO IDs in part (b). How many GO children are there?

#children for GO numbers
go.Children <- getGOChildren(gonr)
#extract GO numbers for children
ch <- sapply(go.Children,function(x) x$Children) 
ch
length(unlist(ch))

# Problem 3 

# load the data 
biocLite("genefilter")
biocLite("genefilter")
library("genefilter")
library("ALL")
data(ALL)
library("limma")

# 3(a) We look for genes expressed differently in stages B2 and B3. Use genefilter to
# program the Wilcoxon test and the Welch t-test separately for each gene. For each
# test, we select the genes with p-value<0.001. To save computational time, we set
# exact=F in the Wilcoxon test function

patient.B <- exprs(ALL)[,(ALL$BT %in% c("B2","B3"))]
factor <- droplevels(ALL$BT[ALL$BT %in% c("B2","B3")])
f1 <- function(x) (wilcox.test(x ~ factor, exact = F)$p.value < 0.001)
f2 <- function(x) (t.test(x ~ factor)$p.value < 0.001) 
wilcox <- genefilter(patient.B, filterfun(f1))
t.test <- genefilter(patient.B, filterfun(f2))

# (b) Compute a Venn diagram for the Wilcoxon test and the t-test, and plot it.

x <- apply(cbind(wilcox,t.test), 2, as.integer)
vc <- vennCounts(x, include="both")
vennDiagram(vc)

# 3(d) What is the annotation package for the ALL data set? Find the GO numbers for "oncogene". 
annotation(ALL)
biocLite("hgu95av2.db")
library(hgu95av2.db)

GO.term <- function(term) {
  GTL <- eapply(GOTERM, function(x) {grep(term, x@Term, value=TRUE)}) 
  Gl <- sapply(GTL, length)         
  names(GTL[Gl>0])                  
}
oncogene.id<- GO.term("oncogene")
oncogene.id

# 3(e) How many genes passing the filters in (a) are oncogenes?
selected <- wilcox & t.test
sel<- patientB[selected,]
probes <- hgu95av2GO2ALLPROBES$"GO:0090402"
print(sum(probes %in% rownames(sel)))

# Problem 4 

# loading the data 

library("limma")
library("ALL")
data(ALL)
library("genefilter")

# 4a)Select the persons with B-cell leukemia which are in stage B1, B2, and B3.

all.B <- ALL[,which(ALL$BT %in% c("B1","B2","B3"))]

# 4(b) Use the linear model to test the hypothesis of all zero group means

fac.B123 <- factor(all.B$BT)
design.ma <- model.matrix(~ 0 + fac.B123)
colnames(design.ma) <- c("B1","B2","B3")
fit <- lmFit(all.B, design.ma); fit <- eBayes(fit)
# p is small. reject null that means = 0; they are different
print(topTable(fit, number=5,adjust.method="fdr"), digits=4)
# differently expressed
sum(topTable(fit, number=Inf,adjust.method="fdr")$adj.P.Val<0.05) 

# B3 group:
print( topTable(fit, coef=3, number=5, adjust.method="fdr"), digits=4) 

# 4(c) Use two contrasts to perform analysis of variance to test the null hypothesis of
# equal group means. Do this with a false discovery rate of 0.01. How many
# differentially expressed genes are found?

ma <- makeContrasts(B1-B2,B2-B3, levels=factor(allB$BT))
fit1 <- contrasts.fit(fit, cont.ma);fit1 <- eBayes(fit1)
fdr.p.data <- topTable(fit1,number=Inf,adjust.method="fdr")$adj.P.Val
# differently expressed
sum(fdr.p.data<0.01) 
print(topTable(fit1, number=5,adjust.method="fdr"), digits=4) 
