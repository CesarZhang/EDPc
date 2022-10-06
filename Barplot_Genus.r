rm(list = ls())

library(stringr)
library(ggplot2)

setwd("~/Code/EDPcG/Microbiota/")

library(paletteer) 
# paletteer_d("basetheme::clean")
cols<-c(rep(as.character(paletteer_d("basetheme::clean")),1),"#B1B2C8FF")

cols1<-c(rep(as.character(paletteer_d("ggthemes::Tableau_20"))[1:19],1),"#B1B2C8FF")
  


##############
meta<-read.csv("Result/meta.tsv",header = T,
               sep = "\t") 
meta$Stage2<-factor(meta$Stage2,levels = c("YL","PM","KK",
                                           "LC","ZJ"))
meta<-meta[order(meta$Stage2),]

data<-read.table("Result/Barplot/L6.tsv",row.names = 1,
                 header = T,check.names = F)
data<-data[,meta$sample.id]

x<-data
x<-x[order(rowSums(x),decreasing = T),]
N<-19
xn<-rownames(x[rownames(x) != c("Other"),])[1:N]
new_x<-rbind(x[rownames(x)%in% xn,],Others=colSums(x[!rownames(x) %in% xn,]))
# new_x<-t(new_x)


######
tax<-str_split_fixed(rownames(new_x),";",6)
tax[tax[,6] == "__",6]<-paste(tax[tax[,6] == "__",6],"Other_in_",
                              tax[tax[,6] == "__",5],sep = "")
tax[tax[,6] == "g__uncultured",6]<-paste(tax[tax[,6] == "g__uncultured",6],"_in_",
                                         tax[tax[,6] == "g__uncultured",5],sep = "")
tax[tax[,6] == "__Other_in___",6]<-paste(tax[tax[,6] == "__Other_in___",6],
                                         tax[tax[,6] == "__Other_in___",3],sep = "")
tax[,6]<-str_replace_all(tax[,6],"g__","")
tax[,6]<-str_replace_all(tax[,6],"f__","")
tax[,6]<-str_replace_all(tax[,6],"c__","")
tax[,6]<-str_replace_all(tax[,6],"__","")
tax[N+1,6]<-tax[N+1,1]
rownames(new_x)<-tax[,6]
new_x<-t(new_x)

x_p<-cbind(as.data.frame(new_x),meta)
x_p<-x_p[order(x_p$Stage2,decreasing = F),]

df<-NULL
#### get the purpose column
p<-3
for (i in 1:(N+1)){
  tmp<-data.frame(row.names=NULL,Sample=rownames(x_p),Taxa=rep(colnames(x_p)[i],dim(x_p)[1]),
                  Value=x_p[,i],Type=x_p[,(N+1+p)])
  #it can't be used? [order(tmp$Type,tmp$Value),]
  if(i==1){df<-tmp} else {df<-rbind(df,tmp)}
}
df$Taxa<-factor(df$Taxa,levels = colnames(x_p)[1:(N+1)])
df$Sample<-factor(df$Sample,levels =rownames(x_p))
df$Type<-factor(df$Type,levels = c("YL","PM","KK",
                                           "LC","ZJ"))

tiff(filename = "barplot_gen1.tiff",
     width = 2800,height = 2200)
ggplot(df,aes(Sample,Value,fill=Taxa))+
  geom_bar(stat="identity",position = "stack")+
  facet_grid(. ~ Type, drop=F,scale="free",space="free_x")+
  scale_fill_manual(values = cols1)+
  theme_classic()+ylab("")+
  guides(fill=guide_legend(title = "Taxonomy",ncol=1))+
  theme(axis.text.x=element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        axis.text= element_text(size=50),
        axis.title.y = element_text(size=60),
        axis.title.x = element_text(size=60),
        axis.line.y = element_blank(),
        legend.key.height = unit(4,"line"),
        legend.text = element_text(size = 60),
        legend.title = element_text(size = 60),
        strip.text = element_text(size=50))
dev.off()
