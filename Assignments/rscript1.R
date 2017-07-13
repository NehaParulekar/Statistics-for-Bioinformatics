array(1:8,c(2,2,2))
array (0, c=(3,4,5))
array (0,dim=c(3,4,5))
mynum=c(2,4)
mynum
mystr=c("I","love","stats")
mystr
myboo=matrix(c(TRUE,FALSE,FALSE,TRUE), ncol=2)
myboo
mylist= list(mynum, mystr, myboo, 100)
mylist
mystr= "hello world!"
myvec= c(4:40)
mylist=list(mystr, myvec)
mylist
mylist[2]
mylist[c(2,4)]
mylist[[2]]
mylist[[2]][1]
mylist[[2]][1]= "you"
mylist[[2]]
mystr
rm(list=ls())
mynum=c(2,4,6)
mystr = c("I","love","stat")
myboo = c(TRUE,FALSE,TRUE)
mydf= data.frame(mynum, mystr, myboo)
mydf
mydf$mystr
mean(mydf$mynum)


gene1<- c(1.00,1.50,1.25)
gene2<- c(1.35,1.55,1.00)
cbind(gene1,gene2)
rbind(gene1,gene2)
gene12<-rbind(gene1,gene2)
gene12
gene3<-c(-1.1,-1.50,-1.25)
gene4<-c(-1.20,-1.30,-1.00)
genedat<-rbind(gene12,gene3,gene4)
colnames(genedat)<-c("Eric","Peter","Anna")
genedat
is.vector(genedat)
is.matrix(genedat)
is.data.frame(genedat)
is.numeric(genedat)
is.character(genedat)
str(genedat)
y<-as.matrix(c(1,2,3,4))
y
as.matrix(c(1,2,3,4),nrow=2)
as.matrix(c(1,2,3,4),ncol=2)
matrix(c(1,2,3,4),nrow=2)
z<-as.character(gene1)
z
mystr
as.numeric(mystr)


missing