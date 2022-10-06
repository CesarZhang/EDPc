rm(list = ls())

setwd("~/Code/EDPcG/Microbiota/")
# setwd("C://Users/Cesar/Documents/Code/YZFP/Barplot/")
library(vegan)
library(ggplot2)

meta<-read.csv("Result/meta.tsv",header = T,
               sep = "\t") 

meta$Stage2<-factor(meta$Stage2,levels = c("YL","PM","KK",
                                           "LC","ZJ"))
meta$Stage1<-factor(meta$Stage1,levels = c("KKQ","KKH"))
meta<-meta[order(meta$Stage2),]

data<-read.csv("Result/Beta/distance-matrix.tsv",
                sep = "\t",row.names = 1,header = T,check.names = F)
data<-data[meta$sample.id,meta$sample.id]

ord<-metaMDS(data)

NMDS<-list()
NMDS[[1]]<-data.frame(x=ord$points[,1],y=ord$points[,2], 
                      meta1= as.factor(meta[,"Stage2"]))

#NMDS[[i]][,3]<-factor(NMDS[[i]][,3],levels = c("Chuff","Oral","Fecal","Rectal","Skin","Freshwater"))
colnames(NMDS[[1]])<-c("x","y",colnames(meta)[3])
color_munual<-c("#FFB415","#3983AD","#64638E","#A43E55","#419583","#78767D")

tiff(filename = paste(colnames(meta)[3], ".tiff",sep = ""),
    width = 2048,height = 2048)
print(ggplot(data = NMDS[[1]],aes(x,y,col=NMDS[[1]][,3]))+
        geom_point(size=40)+
        scale_color_manual(values=color_munual)+
        theme_light()+
        guides(color = guide_legend(""))+
        theme(title= element_text(size=80),
              axis.text= element_text(size=80),
              legend.text = element_text(size = 100))+
        
        labs(x="NMDS1",y="NMDS2"))
dev.off()

