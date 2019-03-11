source("https://bioconductor.org/biocLite.R")
biocLite("org.Hs.eg.db")
library(biomaRt)
library(varhandle)
library(AnnotationDbi)
library(annotate)
library(org.Hs.eg.db)
library(XML)
library(readxl)
##################################################
#targetscan
setwd("E:\\thesis\\database\\TG\\targetscan")
files <- list.files(pattern=".txt")
targetscan<-read.delim(files)
col1<-unfactor(Conserved_Site_Context_Scores[,5])
col2<-unfactor(Conserved_Site_Context_Scores[,2])
tsize=length(col1)
targetscan=matrix(0,tsize,2)
targetscan[,1]=col1
targetscan[,2]=col2
targetscan2=matrix(0,tsize,2)
ci=1
for(i in 1:nrow(targetscan)){
  if(grepl("hsa", targetscan[i,1])){
    z<-mget(targetscan[i,2], revmap(org.Hs.egSYMBOL), ifnotfound=NA)
    targetscan2[ci,2]<-matrix(unlist(z[1][[1]][1]), byrow=TRUE, nrow=length(z) )
    targetscan2[ci,1]=targetscan[i,1]
    ci=ci+1}
}
write.table(targetscan, file="targetscan_2v.txt", row.names=FALSE, col.names=FALSE)
###########################################################
setwd("E:\\thesis\\database\\TG")
mira=matrix(0,nrow(miranda),2)
##miranda
for(i in 1:nrow(miranda)){
  x=unlist(gregexpr(pattern ='-',miranda[i,2]))
  if(length(x)>=3){
    miranda[i,2]=substr(miranda[i,2], 1, x[3]-1)}}
mira1=unfactor(miranda[,2])
mira2=unfactor(miranda[,4])
mira[,1]=mira1
mira[,2]=mira2
write.table(mira, file="miranda.txt", row.names=FALSE, col.names=FALSE)
for(i in 1:nrow(mira)){
  mira[i,2]<-mget(mira[i,2], revmap(org.Hs.egSYMBOL), ifnotfound=NA)
}
##################################
##mirtarbase
setwd("E:\\thesis\\database\\TG\\mirtarbase")
file <- list.files(pattern=".xls")
tmp1<-read_excel(file[1])
col11<-tmp1[,2]
col22<-tmp1[,5]
ttsize=length(col11)
mirtarbase=matrix(0,ttsize,2)
mirtarbase[,1]=col11
mirtarbase[,2]=col22
write.table(mirtarbase, file="mirtarbase_v2.txt", row.names=FALSE, col.names=FALSE)
ci=1
mirtarb=matrix(0,nrow(mirtarbase),2)
for(i in 1:nrow(mirtarbase)){
  if(grepl("hsa", mirtarbase[i,1])){
    mirtarb[ci,]<-mirtarbase[i,]
    ci=ci+1}
}
######################################################
##pictar-dorina
size=0
setwd("E:\\thesis\\database\\TG\\HomoSapiens_mRNA_ENSEMBL78")
file_name=list.files(pattern=".txt")
nfile=length(file_name)
for(i in 1:nfile){
  setwd("E:\\thesis\\database\\TG\\HomoSapiens_mRNA_ENSEMBL78")
  df1 <- read.delim(file_name[i])
  size=nrow(df1)
  v<-gsub("_","-",file_name[i])
  if(grepl("hsa", mir)){
   x=unlist(gregexpr(pattern ='-',v))
   if(length(x)>=3){
     mir=substr(v, 1, x[4]-1)}else{
       mir=v
     }
  x=matrix(0,size,1)
  x[,1]<- sapply(unfactor(df1[,2]), function(k) strsplit(k, "_")[[1]][1])
  z<-mget(x[,1], revmap(org.Hs.egENSEMBL), ifnotfound=NA)
  setwd("E:\\thesis\\database\\TG\\pic")
  out = file(paste(mir,".txt", sep=""),'a')
  z<-matrix(unlist(z), byrow=TRUE, nrow=length(z) )
  write(unique(z), file=out, append=T)
  close(out)}
}
#####################################################
#mirbase
miRDB <- read.delim("E:/thesis/database/TG/miRBase/miRDB.txt", header=FALSE)
mirba=matrix(0,nrow(miRDB),2)
mirba[,1]<-unfactor(miRDB[,1])
mirba[,2]<-unfactor(miRDB[,2])
x=which(grepl("hsa", mirba[,1]))
mirbase=matrix(0,length(x),2)
col111<-mirba[x,1]
col222<-mirba[x,2]
mirbase[,1]=col111
mirbase[,2]=col222
z<-mget(mirbase[,2], revmap(org.Hs.egACCNUM), ifnotfound=NA)
xz<-matrix(unlist(z), byrow=TRUE, nrow=length(z) )
nz<-names(z)
nam<-mirbase[which(xz!="Na"),1]
n=length(unique(nam))
m<-unique(nam)
setwd("E:\\thesis\\database\\TG\\miRBase\\TG_mir")
for(i in 1:n){
  mir<-gsub("_","-",m[i])
  x=which(mirbase[,1]==m[i])
  if(length( xz[x])>0){
    x<-na.omit(xz[x])
    out = file(paste(mir,".txt", sep=""),'a')
    write(unique(x), file=out, append=T)
    close(out)
  }
}
#####################################################
## find TG
setwd("E:\\thesis\\database\\TG\\pic")
file_name=list.files(pattern=".txt")
nfile=length(file_name)-1
pictar=matrix(0,nfile,1)
pictar[,1]=file_name[1:nfile]
mirna=matrix(0,nfile,1)
pictar[,1]<-gsub("-vs","",pictar[,1])
mirna[,1]<- sapply(pictar[,1], function(k) strsplit(k, "\\.")[[1]][1])
setwd("E:\\thesis\\database\\TG\\miRBase\\TG_mir")
file_name1=list.files(pattern=".txt")
nfile1=length(file_name1)
mirbase=matrix(0,nfile1,1)
mirbase[,1]=file_name1
mirb=matrix(0,nfile1,1)
mirb[,1]<- sapply(mirbase[,1], function(k) strsplit(k, "\\.")[[1]][1])
tgs<-unique(targetscan[,1])
mtg<-mirna[,1]
mi<-unique(unfactor(mira[,1]))
mtb<-unique(mirtarbase[,1])
mrb<-unique(mirb[,1])
mir<-unique(c(tgs,mtg,mtb,mi,mrb))
mir<-mir[2:length(mir)]
n=length(mir)
for(i in 1:n){
  ind=0
  setwd("E:\\thesis\\database\\TG\\pic")
  x=which(file_name[]==paste(mir[i], ".txt", sep=""))
  gpictar=NULL
  if(length(x)>0){
    df1 <- read.delim(file_name[x])
    gpictar<-unique(df1[,1])
    gpictar<-na.omit(gpictar)
  }
  setwd("E:\\thesis\\database\\TG\\miRBase\\TG_mir")
  x=which(file_name1[]==paste(mir[i], ".txt", sep=""))
  gmirbase=NULL
  if(length(x)>0){
    df1 <- read.delim(file_name1[x])
    gmirbase<-unique(df1[,1])
    gmirbase<-na.omit(gmirbase)
  }
  x=which(targetscan[,1]==mir[i])
  gtarget<-targetscan[x,2]
  gtarget<-na.omit(gtarget)
  x=which(mira[,1]==mir[i])
  gmira<-mira[x,2]
  gmira<-na.omit(gmira)
  x=which(mirtarbase[,1]==mir[i])
  gmirtarbase<-mirtarbase[x,2]
  gmirtarbase<-na.omit(gmirtarbase)
  gene<-unique(c(gpictar,gtarget,gmira,gmirtarbase,gmirbase))
  size=length(gene)
  genes=matrix(0,size,1)
  x=matrix(0,size,1)
  x[which(gene%in%gtarget),1]=1
  y=matrix(0,size,1)
  y[which(gene%in%gmirtarbase),1]=1
  w=matrix(0,size,1)
  w[which(gene%in%gmira),1]=1
  z=matrix(0,size,1)
  z[which(gene%in%gpictar),1]=1
  q=matrix(0,size,1)
  q[which(gene%in%gmirbase),1]=1
  cg<-x[,1]+y[,1]+z[,1]+w[,1]+q[,1]
  genes<-gene[which(cg>=3)]
  if(length(genes)>=5){
    setwd("E:\\thesis\\database\\TG\\3_targets")
    out = file(paste(mir[i], ".txt", sep=""),'a')
    write(unique(genes), file=out, append=T)
    close(out)}
}

###########################
#mapping mirna
setwd("E:\\thesis\\database\\TG\\targets")
file_name=list.files(pattern=".txt")
n=length(file_name)
mirs=NULL
for(i in 1:n){
  setwd("E:\\thesis\\database\\TG\\targets")
  print(i)
  df1 <- read.delim(file_name[i])
  mir<-mapping(unfactor(df1[,1]),ppi)
  setwd("E:\\thesis\\database\\TG\\targets_map")
  out = file(paste(file_name[i], sep=""),'a')
  write(unique(mir), file=out, append=T)
  close(out)
  mirs=c(mirs,mir)
}
mirs=unique(mirs)

##################################
setwd("E:\\thesis\\database\\TG\\targets_map")
file_name=list.files(pattern=".txt")
n=length(file_name)
mirs=NULL
for(i in 1:n){
  setwd("E:\\thesis\\database\\TG\\targets_map")
  if(file.size(file_name[i])>0){
    df1 <- read.delim(file_name[232])
    mir<-unfactor(df1[,1])
  mirs=c(mirs,mir)}
}
mirs=unique(mirs)

setwd("E:\\thesis\\database\\TG")
out = file("mirs.txt",'a')
write(mirs, file=out, append=T)
close(out)
#############################################
setwd("E:\\thesis\\database\\TG\\3_targets")
file_name=list.files(pattern=".txt")
n=length(file_name)
ci=0
r=NULL
x=which(grepl("-3p", file_name))
for(i in 1:length(x)){
  y=gsub("-3p","", file_name[x[i]])
  z=gsub("-5p","", file_name[x[i]+1])
  if(z==y){
    ci=ci+1
    r=c(r,x[i],x[i]+1)
    setwd("E:\\thesis\\database\\TG\\3_targets")
    df1 <- read.delim(file_name[x[i]])
    df2 <- read.delim(file_name[x[i]+1])
  df1<-matrix(unlist(df1), byrow=TRUE, nrow=length(df1) )
  df2<-matrix(unlist(df2), byrow=TRUE, nrow=length(df2) )
  z<-gsub(".txt","",file_name[x[i]] )
  write.csv(df1,paste("E:\\thesis\\database\\TG\\2_targets\\",z,".csv",sep=""))
  z<-gsub(".txt","",file_name[x[i]+1] )
  write.csv(df2,paste("E:\\thesis\\database\\TG\\2_targets\\",z,".csv",sep=""))
  }
}
for(i in 1:n){
  if(!(i %in% r)){
    setwd("E:\\thesis\\database\\TG\\3_targets")
    df1 <- read.delim(file_name[i])
    y=gsub("-3p","", file_name[i])
    y=gsub("-5p","", y)
    y<-gsub(".txt","", y)
    df1<-matrix(unlist(df1), byrow=TRUE, nrow=length(df1) )
    write.csv(df1,paste("E:\\thesis\\database\\TG\\1_targets\\",y,".csv",sep=""))
  }
}
