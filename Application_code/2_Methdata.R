memory.limit(10^10)
library(data.table)

# Working directory
setwd('D:/Personal working files/Lab/Application_final_v')

############# Meth data reading #################
LUAD.meth<-fread('TCGA-LUAD.methylation450.tsv',header=T)
LUSC.meth<-fread('TCGA-LUSC.methylation450.tsv',header=T)
#LUAD.meth<-(readr::read_delim('TCGA-LUAD.methylation450.tsv\t", escape_double = FALSE, trim_ws = TRUE))
#LUSC.meth<-(readr::read_delim('TCGA-LUSC.methylation450.tsv','\t', escape_double = FALSE, trim_ws = TRUE))
LUSC.meth=LUSC.meth[,-1]
Meth<-cbind(LUAD.meth,LUSC.meth)
rm(LUAD.meth,LUSC.meth)


sample=colnames(Meth)
Meth=na.omit(Meth)
site=Meth[,1]
Meth=Meth[,-1]
Meth=as.data.frame(lapply(Meth,as.numeric))
Meth=t(Meth)
colnames(Meth)=t(site)
rownames(Meth)=sample[-1]
length(rownames(Meth))
ID=data.frame(sample[-1])
colnames(ID)='ID'


#Meth=Meth[,-which(colSums(Meth<0.1)>915*0.90)]

#Meth=Meth[,-which(colSums(Meth>0.9)>915*0.90)]
Meth=cbind(ID,Meth)

#which(colnames(Meth)=='cg04767756')

