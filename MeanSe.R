rm(list = ls())

library(stringr)
library(dplyr)

setwd("~/Code/EDPcG/Microbiota/")
# setwd("M://Code/IPRS/")

meta<-read.csv("Result/meta.tsv",header = T,
               sep = "\t") 
meta$Stage2<-factor(meta$Stage2,levels = c("YL","PM","KK",
                                           "LC","ZJ"))
meta<-meta[order(meta$Stage2),]


##############Phylumn
data<-read.table("Result/Barplot/L2.tsv",row.names = 1,
                 header = T,check.names = F)
data<-data[,meta$sample.id]
data<-data[rowSums(data) >0,]

stde <- function(x) sd(x)/sqrt(length(x))
i<-3
# for(i in 2:ncol(metadata)){
new<-cbind(t(data),meta[,i,drop=F])
colnames(new)[ncol(new)]<-c("Type")

#mean
meann<-new %>% group_by(Type) %>% summarise_all(list(mean = mean)) %>%
  as.data.frame()
rownames(meann)<-meann$Type
meann<-meann[,2:ncol(meann)] %>%
  t() %>% as.data.frame()
rownames(meann)<-str_replace_all(rownames(meann),"_mean","")

#se

se<-new %>% group_by(Type) %>% summarise_all(list(se = stde)) %>%
  as.data.frame()
rownames(se)<-se$Type
se<-se[,2:ncol(se)] %>% 
  t() %>% as.data.frame()
rownames(se)<-str_replace_all(rownames(se),"_se","")


x<-as.data.frame(t(data)) %>%
  cbind(meta[,i,drop=F])
nr<-nrow(x)
pv<-data.frame()
for (i in 1:(ncol(x)-1)) {
  tmp<-x[1:nr,c(i,ncol(x))]
  colnames(tmp)<-c("v","t")
  kt<-kruskal.test(v ~ t, data=tmp)
  pv[i,1]<-kt$p.value
}
rownames(pv)<-colnames(x)[1:(ncol(x)-1)]

###### merge
all<-merge(meann,se,by=0,all=T) %>%
  merge(pv,by.x="Row.names",by.y = 0,all = T)

n<-ncol(all)
for (k in 1:ncol(meann)) {
  all[,(n+k)]<-paste(100*signif(all[,(k+1)],2),"+",
                     100*signif(all[,(ncol(meann)+k+1)],2),
                     "%",sep="")
}
n<-unique(meta$Stage2)
colnames(all)<-c("Taxonomy",
                 paste(n,"_mean",sep = ""),
                 paste(n,"_se",sep = ""),"p",
                 paste(n,"_mean+se",sep = ""))
write.csv(all,"Phy_meanse.csv")

###########Genus
data<-read.table("Result/Barplot/L6.tsv",row.names = 1,
                 header = T,check.names = F)
data<-data[,meta$sample.id]
data<-data[rowSums(data) >0,]

stde <- function(x) sd(x)/sqrt(length(x))
i<-3
# for(i in 2:ncol(metadata)){
new<-cbind(t(data),meta[,i,drop=F])
colnames(new)[ncol(new)]<-c("Type")

#mean
meann<-new %>% group_by(Type) %>% summarise_all(list(mean = mean)) %>%
  as.data.frame()
rownames(meann)<-meann$Type
meann<-meann[,2:ncol(meann)] %>%
  t() %>% as.data.frame()
rownames(meann)<-str_replace_all(rownames(meann),"_mean","")

#se

se<-new %>% group_by(Type) %>% summarise_all(list(se = stde)) %>%
  as.data.frame()
rownames(se)<-se$Type
se<-se[,2:ncol(se)] %>% 
  t() %>% as.data.frame()
rownames(se)<-str_replace_all(rownames(se),"_se","")


x<-as.data.frame(t(data)) %>%
  cbind(meta[,i,drop=F])
nr<-nrow(x)
pv<-data.frame()
for (i in 1:(ncol(x)-1)) {
  tmp<-x[1:nr,c(i,ncol(x))]
  colnames(tmp)<-c("v","t")
  kt<-kruskal.test(v ~ t, data=tmp)
  pv[i,1]<-kt$p.value
}
rownames(pv)<-colnames(x)[1:(ncol(x)-1)]

###### merge
all<-merge(meann,se,by=0,all=T) %>%
  merge(pv,by.x="Row.names",by.y = 0,all = T)

n<-ncol(all)
for (k in 1:ncol(meann)) {
  all[,(n+k)]<-paste(100*signif(all[,(k+1)],2),"+",
                     100*signif(all[,(ncol(meann)+k+1)],2),
                     "%",sep="")
}
n<-unique(meta$Stage2)
colnames(all)<-c("Taxonomy",
                 paste(n,"_mean",sep = ""),
                 paste(n,"_se",sep = ""),"p",
                 paste(n,"_mean+se",sep = ""))
write.csv(all,"Genus_meanse.csv")
