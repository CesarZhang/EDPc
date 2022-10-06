##############
library(phyloseq)
library(ggClusterNet)
library(stringr)
library(tidyr)
library(igraph)
library(dplyr)
library(ggplot2)
library(sna)
library(ggrepel)
library(paletteer) 


rm(list = ls())
setwd("~/Code/EDPcG/Microbiota/")

##############
kegg<-read.csv("KEGG-Pathway.csv")
kegg<-kegg[order(kegg$Level.1,kegg$num,decreasing = T),]
kegg$Level.2<-factor(kegg$Level.2,levels = rev(kegg$Level.2))

col1<-c(rep(as.character(paletteer_d("ggthemes::Tableau_20"))[1:19],1),"#B1B2C8FF")


tiff("Kegg.tiff",width = 2140,height = 1980)
ggplot()+
  geom_bar(aes(y=Level.2,x=num,fill=Level.1),data = kegg,
           stat="identity")+
  scale_color_manual(aesthetics = c("fill"),values = col1)+
  labs(x="Number of Genes",y="")+
  guides(fill=guide_legend(title=""))+
  theme(legend.text = element_text(size = 30),
        legend.key = element_blank(),
        axis.title.x = element_text(size = 40),
        axis.text = element_text(size = 35),
        # axis.ticks = element_line(size = 2),
        # axis.ticks.length = unit(.5, "cm"),
        legend.title = element_text(size = 30),
        legend.key.size = unit(2, 'cm'),
        legend.key.height = unit(2, 'cm'),
        legend.key.width  = unit(2, 'cm'))

dev.off()  

####################
metadata <- read.csv("Result/meta.tsv",row.names = 1,header = T,sep = "\t")
otutab <- read.csv("Result/Relative-table.tax.tsv",row.names = 1,sep = "\t",
                   header = T)
taxo<-separate(otutab,taxonomy,into = c("Kingdom","Phylum","Class","Order","Family","Genus","Species"),sep = "; ")
taxo[is.na(taxo)] <- NA

taxo[,16]<-str_replace_all(taxo[,16],"d__","")
taxo[,17]<-str_replace_all(taxo[,17],"p__","")
taxo[,18]<-str_replace_all(taxo[,18],"c__","")
taxo[,19]<-str_replace_all(taxo[,19],"o__","")
taxo[,20]<-str_replace_all(taxo[,20],"f__","")
taxo[,21]<-str_replace_all(taxo[,21],"g__","")
taxo[,22]<-str_replace_all(taxo[,22],"s__","")
taxo<-taxo[,16:22]
sp<-c("Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium")
taxo[taxo[,6] %in% sp ,6]<-paste("Allorhizobium-Neorhizobium-","\n",
                                  "Pararhizobium-Rhizobium",sep = "")


otutree<-read_tree("Result/tree.nwk")

ps = phyloseq(sample_data(metadata),
              otu_table(as.matrix(otutab[,1:(ncol(otutab)-1)]), taxa_are_rows=TRUE),
              tax_table(as.matrix(taxo)),
              phy_tree(otutree)
)
tdata<-read.csv("../Transcriptome/Result/Gene_Exp1.csv",row.names = 1,header = T)
tdata<-tdata[order(rowSums(tdata),decreasing = T),][1:2208,]

tat<-read.csv("../Transcriptome/Result/Gene_Exp_tax.csv",row.names = 1,header = T)



#################
##Immune system
kt<-read.csv("Immune-system.csv")
tdata1<-tdata[kt$GeneID,]
tat1<-tat[kt$GeneID,]
ts = phyloseq(sample_data(metadata),
              otu_table(tdata1, taxa_are_rows=T),
              tax_table(as.matrix(tat1))
              # phy_tree(otutree)
)

ps.merge <- ggClusterNet::merge16S_ITS(ps16s = ps,
                                       psITS = ts,
                                       N16s = 370,
                                       NITS = 2088
)

map =  phyloseq::sample_data(ps.merge)
map$Group = "one"
phyloseq::sample_data(ps.merge) <- map

set.seed(925)
result <- corBionetwork(ps = ps.merge,
                        N = 0,
                        # lab = data,
                        r.threshold = 0.8, # 相关阈值
                        p.threshold = 0.05,
                        group = "Group",
                        # env = data1, # 环境指标表格
                        # envGroup = Gru,# 环境因子分组文件表格
                        # layout = "fruchtermanreingold",
                        # path = Envnetplot,# 结果文件存储路径
                        # fill = "Phylum", # 出图点填充颜色用什么值
                        size = "igraph.degree", # 出图点大小用什么数据
                        scale = F, # 是否要进行相对丰度标准化
                        bio = TRUE, # 是否做二分网络
                        zipi = F, # 是否计算ZIPI
                        step = 100, # 随机网络抽样的次数
                        width = 12,
                        label = TRUE,
                        height = 10,
                        big = TRUE,
                        select_layout = TRUE,
                        layout_net = "model_maptree2",
                        clu_method = "cluster_fast_greedy"
                        
                        
)
dev.off()

col<-c("#FFB415","#3983AD")
edges<-result[[3]]
nodeG<-result[[4]] %>%
  subset(igraph.degree>0)
# nodeG<-subset(nodeG,(group %in% c("bac") & igraph.degree>10) | group %in% c("fun"))

nodeB<-subset(nodeG,group %in% "bac")
nodeB[is.na(nodeB$Genus),"Genus"]<-paste("Other_in_",
                                         nodeB[is.na(nodeB$Genus),"Family"],
                                         sep = "")
p<-ggplot() + geom_segment(aes(x = X1, y = Y1, xend = X2, yend = Y2,color = as.factor(cor)),
                           data = edges, size = 1.5,alpha = 0.5) +
  scale_colour_brewer(aesthetics = c("color"),name="Correlation",palette = "Set1") +
  scale_size_continuous(range = c(10, 50),name = "Igraph.degree")+
  geom_point(aes(x = X1, y = X2,size = igraph.degree,fill = group),pch = 21, data =  nodeG) +
  scale_color_manual(aesthetics = c("fill"),name="Group",values = col,
                     labels = c("Bacteria","Gene"))+
  scale_x_continuous(breaks = NULL) + scale_y_continuous(breaks = NULL) +
  geom_text_repel(aes(x = X1, y = X2,label=Genus),data = nodeB,size=8,segment.size=1,
                  box.padding = 1.5, max.overlaps = Inf)+
  # labs( title = paste(layout,"network",sep = "_")) + 
  guides(fill=guide_legend(nrow=6,byrow=F,override.aes = list(size=10)),
         color=guide_legend(override.aes = list(size=1.5)))+
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
  theme(legend.box = "vertical")+
  theme(legend.text = element_text(size = 30),
        legend.key = element_blank(),
        legend.key.size = unit(2, 'cm'),
        legend.key.height = unit(2, 'cm'),
        legend.key.width  = unit(2, 'cm'),
        legend.title = element_text(size = 30))+
  theme(legend.background = element_rect(colour = NA)) +
  theme(panel.background = element_rect(fill = "white",  colour = NA)) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())

tiff(filename = "network_is.tiff",
     width = 1980,height = 1460)
print(p)
dev.off()
###########
write.csv(result[[4]],"nodeG-is.csv")
write.csv(result[[3]],"edge-is.csv")
write.csv(result[[5]],"cor-is.csv")

###########
nod<-subset(nodeG,filed == "bac" & igraph.degree >0)
nod$label<- nod$Genus
nod[is.na(nod$label),"label"] <- paste("Other_in_",
                                       nod[is.na(nod$label),"Family"],sep = "")

gn<-nod |>
  group_by(label) |>
  summarise(
    count=n())

gn<-merge(gn,nod[,c("label","Phylum")],by.x = "label") %>% unique()
gn<-gn[order(gn$count,decreasing = T),]
gn$label<-factor(gn$label,levels = rev(gn$label))

tiff("ASV.count.is.tiff",width = 1080,height = 980)
ggplot()+
  geom_bar(aes(y=label,x=count,fill=Phylum),data = gn,
           stat="identity")+
  scale_x_continuous(limits = c(0,15))+
  scale_color_manual(aesthetics = c("fill"),values = col)+
  labs(x="ASV Count",y="")+
  theme(legend.text = element_text(size = 30),
        legend.key = element_blank(),
        axis.title.x = element_text(size = 40),
        axis.text = element_text(size = 30),
        # axis.ticks = element_line(size = 2),
        # axis.ticks.length = unit(.5, "cm"),
        legend.title = element_text(size = 30),
        legend.key.size = unit(2, 'cm'),
        legend.key.height = unit(2, 'cm'),
        legend.key.width  = unit(2, 'cm'))

dev.off()  
###############################
#####Venn
library(VennDiagram)

set1<-read.csv("Carbohydrate_metabolism.csv")[,1]
set2<-read.csv("Amino_acid_metabolism.csv")[,1]
set3<-read.csv("Lipid_metabolism.csv")[,1]

color_m<-c("#FFB415","#3983AD","#78767D")
venn.diagram(
  x = list(set1, set2, set3),
  category.names = c("Carbohydrate\nmetabolism" , "Amino acid\nmetabolism" , "Lipid\nmetabolism"),
  filename = 'metabolism.png',
  output=TRUE,
  
  # Output features
  imagetype="png" ,
  height = 980 , 
  width = 980 , 
  resolution = 600,
  compression = "lzw",
  
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = color_m,
  
  # Numbers
  cex = .6,
  # fontface = "bold",
  fontfamily = "sans",
  
  # Set names
  cat.cex = .5,
  # cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135),
  cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans",
  rotation = 1
)
##############################
##Metabolism
kt<-read.csv("Metabolism.csv")
tdata1<-tdata[kt$Gene.ID,]
tat1<-tat[kt$Gene.ID,]
ts = phyloseq(sample_data(metadata),
              otu_table(tdata1, taxa_are_rows=T),
              tax_table(as.matrix(tat1))
              # phy_tree(otutree)
)

ps.merge <- ggClusterNet::merge16S_ITS(ps16s = ps,
                                       psITS = ts,
                                       N16s = 370,
                                       NITS = 46
)

map =  phyloseq::sample_data(ps.merge)
map$Group = "one"
phyloseq::sample_data(ps.merge) <- map


result <- corBionetwork(ps = ps.merge,
                        N = 0,
                        # lab = data,
                        r.threshold = 0.8, # 相关阈值
                        p.threshold = 0.05,
                        group = "Group",
                        # env = data1, # 环境指标表格
                        # envGroup = Gru,# 环境因子分组文件表格
                        # layout = "fruchtermanreingold",
                        # path = Envnetplot,# 结果文件存储路径
                        # fill = "Phylum", # 出图点填充颜色用什么值
                        size = "igraph.degree", # 出图点大小用什么数据
                        scale = F, # 是否要进行相对丰度标准化
                        bio = TRUE, # 是否做二分网络
                        zipi = F, # 是否计算ZIPI
                        step = 100, # 随机网络抽样的次数
                        width = 12,
                        label = TRUE,
                        height = 10,
                        big = TRUE,
                        select_layout = TRUE,
                        layout_net = "model_maptree2",
                        clu_method = "cluster_fast_greedy"
                        
                        
)
dev.off()

col<-c("#FFB415","#3983AD")
edges<-result[[3]]
nodeG<-result[[4]] %>%
  subset(igraph.degree>0)
nodeG[is.na(nodeG$Genus),"Genus"]<-paste("Other_in_",
                                         nodeB[is.na(nodeG$Genus),"Family"],
                                         sep = "")
write.csv(nodeG,"nodeG-mb.csv")
nodeG<-read.csv("nodeG-mb.csv",row.names = 1)

p<-ggplot() + geom_segment(aes(x = X1, y = Y1, xend = X2, yend = Y2,color = as.factor(cor)),
                           data = edges, size = 1.5,alpha = 0.5) +
  scale_colour_brewer(aesthetics = c("color"),name="Correlation",palette = "Set1") +
  scale_size_continuous(range = c(10, 50),name = "Igraph.degree")+
  geom_point(aes(x = X1, y = X2,fill = group),size=35,pch = 21, data =  nodeG) +
  scale_color_manual(aesthetics = c("fill"),name="Group",values = col,
                     labels = c("Bacteria","Gene"))+
  scale_x_continuous(breaks = NULL) + scale_y_continuous(breaks = NULL) +
  geom_text_repel(aes(x = X1, y = X2,label=Genus),data = nodeG,size=8,segment.size=1,
                  box.padding = 0, max.overlaps = Inf)+
  # labs( title = paste(layout,"network",sep = "_")) + 
  guides(fill=guide_legend(nrow=6,byrow=F,override.aes = list(size=10)),
         color=guide_legend(override.aes = list(size=1.5)))+
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
  theme(legend.box = "vertical")+
  theme(legend.text = element_text(size = 30),
        legend.key = element_blank(),
        legend.key.size = unit(2, 'cm'),
        legend.key.height = unit(2, 'cm'),
        legend.key.width  = unit(2, 'cm'),
        legend.title = element_text(size = 30))+
  theme(legend.background = element_rect(colour = NA)) +
  theme(panel.background = element_rect(fill = "white",  colour = NA)) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())

tiff(filename = "network_Mb.tiff",
     width = 1080,height = 980)
print(p)
dev.off()
################
##Metabolism-all
kt<-read.csv("Metabolism-all.csv")
tdata1<-tdata[kt$Gene.ID,]
tat1<-tat[kt$Gene.ID,]
ts = phyloseq(sample_data(metadata),
              otu_table(tdata1, taxa_are_rows=T),
              tax_table(as.matrix(tat1))
              # phy_tree(otutree)
)

ps.merge <- ggClusterNet::merge16S_ITS(ps16s = ps,
                                       psITS = ts,
                                       N16s = 370,
                                       NITS = 1528
)

map =  phyloseq::sample_data(ps.merge)
map$Group = "one"
phyloseq::sample_data(ps.merge) <- map

set.seed(9251)
result <- corBionetwork(ps = ps.merge,
                        N = 0,
                        # lab = data,
                        r.threshold = 0.8, # 相关阈值
                        p.threshold = 0.05,
                        group = "Group",
                        # env = data1, # 环境指标表格
                        # envGroup = Gru,# 环境因子分组文件表格
                        # layout = "fruchtermanreingold",
                        # path = Envnetplot,# 结果文件存储路径
                        # fill = "Phylum", # 出图点填充颜色用什么值
                        size = "igraph.degree", # 出图点大小用什么数据
                        scale = F, # 是否要进行相对丰度标准化
                        bio = TRUE, # 是否做二分网络
                        zipi = F, # 是否计算ZIPI
                        step = 100, # 随机网络抽样的次数
                        width = 12,
                        label = TRUE,
                        height = 10,
                        big = TRUE,
                        select_layout = TRUE,
                        layout_net = "model_maptree2",
                        clu_method = "cluster_fast_greedy"
                        
                        
)
dev.off()

col<-c("#FFB415","#3983AD")
edges<-result[[3]]
nodeG<-result[[4]] %>%
  subset(igraph.degree>0)
# nodeG<-subset(nodeG,(group %in% c("bac") & igraph.degree>10) | group %in% c("fun"))
write.csv(result[[4]],"nodeG-mta.csv")
write.csv(result[[3]],"edge-mta.csv")
write.csv(result[[5]],"cor-mta.csv")

nodeB<-read.csv("nodeB-mta.csv",row.names = 1)
p<-ggplot() + geom_segment(aes(x = X1, y = Y1, xend = X2, yend = Y2,color = as.factor(cor)),
                           data = edges, size = 1.5,alpha = 0.5) +
  scale_colour_brewer(aesthetics = c("color"),name="Correlation",palette = "Set1") +
  scale_size_continuous(range = c(10, 50),name = "Igraph.degree")+
  geom_point(aes(x = X1, y = X2,size = igraph.degree,fill = group),pch = 21, data =  nodeG) +
  scale_color_manual(aesthetics = c("fill"),name="Group",values = col,
                     labels = c("Bacteria","Gene"))+
  scale_x_continuous(breaks = NULL) + scale_y_continuous(breaks = NULL) +
  geom_text_repel(aes(x = X1, y = X2,label=Genus),data = nodeB,size=8,segment.size=1,
                  box.padding = 1.5, max.overlaps = Inf)+
  # labs( title = paste(layout,"network",sep = "_")) + 
  guides(fill=guide_legend(nrow=6,byrow=F,override.aes = list(size=10)),
         color=guide_legend(override.aes = list(size=1.5)))+
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
  theme(legend.box = "vertical")+
  theme(legend.text = element_text(size = 30),
        legend.key = element_blank(),
        legend.key.size = unit(2, 'cm'),
        legend.key.height = unit(2, 'cm'),
        legend.key.width  = unit(2, 'cm'),
        legend.title = element_text(size = 30))+
  theme(legend.background = element_rect(colour = NA)) +
  theme(panel.background = element_rect(fill = "white",  colour = NA)) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())

tiff(filename = "network_metabolism_all.tiff",
     width = 1980,height = 1460)
print(p)
dev.off()
################
##############
#heatmap
library(pheatmap)
lab<-nodeB
lab["bac_08c1a1f473cbc677485d3c3ad30ee762","Genus"]<-c("ANPR")
rownames(lab)<-str_replace(rownames(lab),"fun_","") %>%
                str_replace("bac_","")
ba<- rownames(subset(lab,group == "bac"))
dfba<- otutab[ba,rownames(metadata)]
lab1<-lab[ba,]

tiff("Metab-b-13.tiff",width = 2340,height = 860,res = 200)
pheatmap(dfba,
         show_rownames = T,show_colnames = T,
         cluster_cols = F,cluster_rows = T,
         fontsize = 15,
         scale = "row",
         angle_col = c("45"),
         # display_numbers = T,
         # border_color = "black",
         labels_row = lab1$Genus,
         cellwidth	= 30,
         cellheight = 15)
dev.off() 
################

ge<- rownames(subset(lab,group == "fun"))
dfge<- tdata[ge,rownames(metadata)]
lab1<-lab[ge,]

tiff("Metab-g-13.tiff",width = 2340,height = 860,res = 200)
pheatmap(dfge,
         show_rownames = T,show_colnames = T,
         cluster_cols = F,cluster_rows = T,
         fontsize = 15,
         scale = "row",
         angle_col = c("45"),
         # display_numbers = T,
         # border_color = "black",
         labels_row = lab1$Genus,
         cellwidth	= 30,
         cellheight = 15)
dev.off() 


################
####Line plot
metadata$Stage2<-factor(metadata$Stage2,levels = c("YL","PM","KK",
                                           "LC","ZJ"))
metadata<-metadata[order(metadata$Stage2),]

bal<-c("8c3ccefdcddf204d68de94534500528f")
gel<-c("LS_GLEAN_10011059")
bal1<-otutab[bal,rownames(metadata)]
gel1<-tdata[gel,rownames(metadata)]
df<-data.frame(Value=paste(c(bal1[1,],gel1[1,])),Stage=rep(rownames(metadata),2),
                Name=c(rep("Flavobacterium",15),rep("ALDH7A1",15)))
ggplot(data=df, aes(x=Stage, y=Value, group=Name)) +
  geom_line()+
  geom_point()
