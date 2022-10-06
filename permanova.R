setwd("~/Code/EDPcG/Microbiota/")
# setwd("M:/Code/IPRS/PERMANOVA/")
rm(list = ls())
library(vegan)
library(dplyr)
meta<-read.csv("Result/meta.tsv",header = T,
               sep = "\t") 
meta$Stage2<-factor(meta$Stage2,levels = c("YL","PM","KK",
                                           "LC","ZJ"))
meta<-meta[order(meta$Stage2),]

data<-read.csv("Result/Beta/distance-matrix.tsv",
               sep = "\t",row.names = 1,header = T,check.names = F)
data<-data[meta$sample.id,meta$sample.id]

#########
for(i in c(3)){
  tmeta<-meta
  tdata<-data[tmeta$sample.id,tmeta$sample.id]
  sink(paste(colnames(tmeta)[i],".txt",sep = ""))  
  print(adonis2(tdata ~ tmeta[,i],data = meta,permutations = 999))
  sink()
}


