## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite("siggenes")
library(siggenes)
#data(golub)
# sam.out <- sam(golub, golub.cl, rand = 123, gene.names = golub.gnames[,3])
# sam.out
setwd("/home/sepideh/Documents/thesis/database/breast")
breast <- as.matrix(read.delim("/home/sepideh/Documents/thesis/database/breast/breast.txt", header=FALSE, stringsAsFactors=FALSE))
net=as.matrix(breast[,1])
breast<-breast[,2:ncol(breast)]
class <- read.csv("/home/sepideh/Documents/thesis/database/breast/class.txt", sep="", stringsAsFactors=FALSE)

y<-as.vector(matrix(0,nrow(class),1))
y[which(grepl("tumor",class[,1]))]=1
y[which(grepl("normal",class[,1]))]=2
y=as.numeric(y)


sam.out <- sam(breast, y, method=d.stat, rand = 123, gene.names = net[,1])
sam.out
#Gene-specific information about the genes called differentially expressed using a specific value of ???
sum.sam.out <- summary(sam.out,1.4)
sum.sam.out

#To obtain just the names of the genes called significant using ??? = 0.3
x=list.siggenes(sam.out, 1.4)
#write.table(x, "../Data//DEG99.txt")

DEGs=matrix(0,length(x),1)
for(i in 1:length(x)){
  g=(which(netid[,1]==x[i]))
  if(length(g>0)){
    DEGs[i,1]=netid[g,2]}
}
DEGs<-unique(DEGs)
DEGs=DEGs[DEGs!=0]




print(sum.sam.out, varNames = "Proteins")

#The rows of golub that contain the values of the differentially expressed genes can also be obtained by:
sum.sam.out@row.sig.genes

#the general information about the set of significant genes by
sum.sam.out@mat.fdr

#and the gene-specific information by
sum.sam.out@mat.sig

#sam.out <- sam(golub, golub.cl, method = wilc.stat, rand = 123)
summary(sam.out)

#and the gene-specific information by
sum.sam.out@mat.sig
