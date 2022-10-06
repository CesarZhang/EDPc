############
# install.packages("ggpubr")
# library(ggplot2)
library(ggpubr)
rm(list = ls())
setwd("~/Code/EDPcG/Microbiota/")
############
meta<-read.csv("Result/meta.tsv",header = T,
               sep = "\t") 
meta$Stage2<-factor(meta$Stage2,levels = c("YL","PM","KK",
                                           "LC","ZJ"))
meta$Stage1<-factor(meta$Stage1,levels = c("KKQ","KKH"))
meta<-meta[order(meta$Stage2),]


aldn<-c("observed_features","shannon","ace","chao1")
ald<-data.frame()
for (i in 1:length(aldn)) {
  path<-"M:/Code/EDPcG/Microbiota/Result/Alpha/"
  if (i == 1) {
   ald<-read.csv(paste(path,aldn[i],"/alpha-diversity.tsv",sep = ""),
                   sep = "\t")
   } else {
   tmp<-read.csv(paste(path,aldn[i],"/alpha-diversity.tsv",sep = ""),
                sep = "\t")
   ald<-merge(ald,tmp,by.x = "X")
   }
  print(aldn[i])
}

######## 1
df<-merge(subset(ald,X %in% meta$sample.id),meta[,c("sample.id","Stage2")],
          by.x="X",by.y = "sample.id")

colnames(df)<-c("ID","observed_features","shannon","ace","chao","Type")
  
  
df1<-NULL
for (i in (3:dim(df)[2])-1){
  tmp<-data.frame(row.names=NULL,Sample=df[,"ID"],
                  Index=rep(colnames(df)[i],dim(df)[1]),
                  Value=df[,i],Type=df[,"Type"])
  print(colnames(df)[i])
  #it can't be used? [order(tmp$Type,tmp$Value),]
  if(i==1){df1<-tmp} else {df1<-rbind(df1,tmp)}
}

my_comparisons_pd<-list( c("YL","PM"),c("PM", "KK"),c("KK", "LC"),
                         c("LC", "ZJ"))
color_munual<-c("#FFB415","#3983AD","#64638E","#A43E55","#419583")

df1$Index<-factor(df1$Index,levels = c("observed_features",
                                       "ace","chao",
                                       "shannon"))

tiff(filename = "alpha3.tiff",
     width = 1240,height = 1240)
ggplot(df1,aes(x=Type,y=Value,color=Type))+
  geom_boxplot(lwd=1.5)+
  theme_bw()+
  ylab("")+
  xlab("")+
  theme(legend.position = "none",plot.title = element_text(size = 8))+
  stat_compare_means(
                     # comparisons = my_comparisons_pd,
                     size=8,
                     # method = "t.test",
                     # symnum.args =list(cutpoints = c(0.0001, 0.001, 0.01, 0.05, 1),
                     #                     symbols = c("***", "**", "*", "ns")),
                     # symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
                     #                     symbols = c("****", "***", "**", "*", "ns")),
                     bracket.size = 5)+

  scale_color_manual(values=color_munual)+
  facet_wrap(. ~ Index, drop=F,scales = "free",nrow = 2)+
  geom_jitter(shape=16,size=7 ,position=position_jitter(0.2))+
  theme(title= element_text(size=40),
        strip.text = element_text(size = 50),
        
        axis.text.y = element_text(size = 40),
        legend.key.height = unit(2.5, 'cm'), 
        legend.key.width  = unit(3.2, 'cm'), 
        legend.text = element_text(size = 40),
        legend.title = element_blank(),
        axis.text.x = element_text(size = 40))
  
dev.off()


######### 2
df<-merge(subset(ald,X %in% meta$sample.id),meta[,c("sample.id","Stage1")],
          by.x="X",by.y = "sample.id")

colnames(df)<-c("ID","observed_features","shannon","faith_pd","Type")


df1<-NULL
for (i in (3:dim(df)[2])-1){
  tmp<-data.frame(row.names=NULL,Sample=df[,"ID"],
                  Index=rep(colnames(df)[i],dim(df)[1]),
                  Value=df[,i],Type=df[,"Type"])
  print(colnames(df)[i])
  #it can't be used? [order(tmp$Type,tmp$Value),]
  if(i==1){df1<-tmp} else {df1<-rbind(df1,tmp)}
}

my_comparisons_pd<-list( c("KKQ","KKH"))
color_munual<-c("#FFB415","#3983AD","#64638E","#A43E55","#419583")

df1$Index<-factor(df1$Index,levels = c("observed_features",
                                       "shannon",
                                       "faith_pd"))

tiff(filename = "alpha2.tiff",
     width = 2200,height = 860)
ggplot(df1,aes(x=Type,y=Value,color=Type))+
  geom_boxplot(lwd=1.5)+
  theme_bw()+
  ylab("")+
  xlab("")+
  theme(legend.position = "bottom",plot.title = element_text(size = 8))+
  stat_compare_means(
                     comparisons = my_comparisons_pd,
                     size=8,
                     # method = "t.test",
                     # symnum.args =list(cutpoints = c(0.0001, 0.001, 0.01, 0.05, 1),
                     #                     symbols = c("***", "**", "*", "ns")),
                     # symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
                     #                     symbols = c("****", "***", "**", "*", "ns")),
                     bracket.size = 2)+
  
  scale_color_manual(values=color_munual)+
  facet_wrap(. ~ Index, drop=F,scales = "free_y",nrow = 1)+
  geom_jitter(shape=16,size=7 ,position=position_jitter(0.2))+
  theme(title= element_text(size=40),
        strip.text = element_text(size = 50),
        
        axis.text.y = element_text(size = 40),
        legend.key.height = unit(2.5, 'cm'), 
        legend.key.width  = unit(3.2, 'cm'), 
        legend.text = element_text(size = 40),
        legend.title = element_blank(),
        axis.text.x = element_blank())

dev.off()

