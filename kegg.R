#mean se 
###########################################

rm(list = ls())
# install.packages("ggplotify")
library(ggplotify)

library(stringr)
library(dplyr)
library(pheatmap)

setwd("~/Code/EDPcG/Microbiota/")
# setwd("M://Code/IPRS/")

meta<-read.csv("Result/meta.tsv",header = T,
               sep = "\t") 
meta$Stage2<-factor(meta$Stage2,levels = c("YL","PM","KK",
                                           "LC","ZJ"))
meta<-meta[order(meta$Stage2),]


##############Phylumn
data<-read.table("Result/KEGG/path_abun_unstrat.tsv",row.names = 1,
                 header = T,check.names = F)
data<-data[,meta$sample.id]
data<-data[rowSums(data) >0,]

###############
x<-t(data)
x<-x/rowSums(x) 
data<-t(x)
########

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
write.csv(all,"kegg_meanse.csv")

###################
kegginf<-read.csv("Result/KEGG/KEGG_pathways_info.tsv",sep = "\t",header = T)
dfk<-subset(all,p<0.05)
dfk<-merge(dfk,kegginf,by.x = "Taxonomy",by.y = "ID")
write.csv(dfk,"kegg_p.csv")
####################
dfk2<-read.csv("kegg_p_r.csv",sep = ",",header = T,
               row.names =1)
dfk2<-dfk2[,c(1:5,17:19)]
dfk2<-dfk2[order(dfk2$Level1,dfk2$Level2,dfk2$Name),]
colnames(dfk2)[1:5]<-c("YL","PM","KK","LC","ZJ")

dfk2m<-subset(dfk2,Level1 %in% c("Metabolism"))


for (i in unique(dfk2m$Level2)) {
  tmp<-subset(dfk2m,Level2 %in% i)
  # rownames(tmp)<-tmp$Name
  tiff(filename = paste(i,".tiff",sep = ""),
       width = 1600,height = 900)
    pheatmap(tmp[,1:5],
             # colorRampPalette(c("white","grey50","grey30","grey10","grey1","black"))(250) ,
             show_rownames = T,show_colnames = T,
             cluster_rows = T,cluster_cols = F,
             scale = "row",
             main = c(i),
             # clustering_distance_rows = "binary",
             # clustering_method = "average",
             treeheight_row = 80,
             treeheight_col = 40,
             legend = T,
             border_color = "white",
             #kmeans_k = 30,
             legend_breaks = c(1.5,0.5,-0.5,-1.5),
             # treeheight_col = 0,
             # annotation_row = dfk2m[,c("Level2"),drop=F],
             # annotation_names_col = T,
             labels_row = tmp$Name,
             cellwidth=50, cellheight=50,
             fontsize_row = 30,
             fontsize_col = 50,
             fontsize = 30)
    dev.off()
}

############################################################

dfk2o<-subset(dfk2,!Level1 %in% c("Metabolism"))
# rownames(dfk2o)<-paste(dfk2o$Level2,"___",dfk2o$Name,sep = "")

tiff("kegg-o2.tiff",width = 2400,height = 1800)
pheatmap(dfk2o[,1:5],
         # colorRampPalette(c("white","grey50","grey30","grey10","grey1","black"))(250) ,
         show_rownames = T,show_colnames = T,
         cluster_rows = T,cluster_cols = F,
         scale = "row",
         # cutree_rows = 7,
         # clustering_distance_rows = "binary",
         # clustering_method = "average",
         treeheight_row = 80,
         treeheight_col = 40,
         legend = T,
         border_color = "white",
         #kmeans_k = 30,
         legend_breaks = c(1.5,0.5,-0.5,-1.5),
         # treeheight_col = 0,
         annotation_row = dfk2o[,c("Level1"),drop=F],
         # annotation_names_col = T,
         labels_row = paste(dfk2o$Level2,"   (",dfk2o$Name,")",sep = ""),
         cellwidth=50, cellheight=50,
         fontsize_row = 30,
         fontsize_col = 50,
         fontsize = 30)
dev.off()
