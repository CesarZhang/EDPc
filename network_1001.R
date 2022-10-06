##############
library(phyloseq)
library(ggClusterNet)
library(stringr)
library(tidyr)
library(igraph)
library(dplyr)
library(ggplot2)
library(sna)
library(paletteer)


rm(list = ls())
setwd("~/Code/EDPcG/Microbiota/")
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

n<-floor(nrow(otutab)/10)
result = corMicro (ps = ps,
                   N = n,
                   method.scale = "TMM",
                   r.threshold=0.8,
                   p.threshold=0.05,
                   method = "spearman"
)

cor = result[[1]]

ps_net = result[[3]]

otu_table = ps_net %>% 
  vegan_otu() %>%
  t() %>%
  as.data.frame()
tax = ps_net %>% vegan_tax() %>%
  as.data.frame()

cols1<-c(rep(as.character(paletteer_d("ggthemes::Tableau_20"))[1:19],1),"#B1B2C8FF")


###################################
set.seed(13)
netClu  = modulGroup( cor = cor,cut = 10,method = "cluster_fast_greedy" )
result2 = model_maptree2(cor = cor,
                         method = "cluster_fast_greedy"
)

# result2 = PolygonRrClusterG (cor = cor,nodeGroup =group2 )
node = result2[[1]]

# ---node节点注释#-----------
nodes = nodeadd(plotcord =node,otu_table = otu_table,tax_table = tax)

nodes2 = nodes %>% inner_join(netClu,by = c("elements" = "ID"))
nodes2$group = paste("Module_",nodes2$group,sep = "")

#-----计算边#--------
edge = edgeBuild(cor = cor,node = node)

### 出图
pnet <- ggplot() + 
  geom_segment(aes(x = X1, y = Y1, xend = X2, yend = Y2,color = as.factor(cor)),
               data = edge, size = 1.5, alpha=0.7) +
  scale_colour_brewer(aesthetics = c("color"),name="Correlation",palette = "Set1") +
  geom_point(aes(X1, X2,fill = group,size = mean),pch = 21, data = nodes2) +
  scale_size_continuous(range = c(0, 50),name = "Mean")+
  scale_color_manual(aesthetics = c("fill"),name="Module",values = cols1)+
  scale_x_continuous(breaks = NULL) + scale_y_continuous(breaks = NULL) +
  theme(panel.background = element_blank()) +
  guides(fill=guide_legend(nrow=6,byrow=F,override.aes = list(size=10)),
         color=guide_legend(override.aes = list(size=1.5)))+
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
  theme(legend.box = "vertical")+
  theme(legend.text = element_text(size = 30),
        legend.key = element_blank(),
        legend.key.size = unit(1, 'cm'),
        legend.key.height = unit(2, 'cm'),
        legend.key.width  = unit(3, 'cm'),
        legend.title = element_text(size = 30))+
  theme(legend.background = element_rect(colour = NA)) +
  theme(panel.background = element_rect(fill = "white",  colour = NA)) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())
# pnet

tiff(filename = "network_1001.tiff",
     width = 1840,height = 1240)
print(pnet)
dev.off()
write.csv(edge,"edge_1001.csv")
write.csv(nodes2,"nodes_1001.csv")


####################################
result4 = nodeEdge(cor = cor)
#提取变文件
edge = result4[[1]]
#--提取节点文件
node = result4[[2]]
igraph  = igraph::graph_from_data_frame(edge, directed = FALSE, vertices = node)

hub = hub_score(igraph)$vector %>%
  sort(decreasing = TRUE) %>%
  as.data.frame()
colnames(hub) = "hub_sca"

hub_2<-merge(hub,nodes2,by.x = 0,by.y = "elements")
hub_2<-hub_2[order(hub_2$hub_sca,decreasing = T),][1:20,]
hub_2$Label<-hub_2$Genus
hub_2[is.na(hub_2$Label),"Label"]<-paste("Other_in_",
                                         hub_2[is.na(hub_2$Label),"Family"],
                                         sep = "")
hub_2$Label<-paste(hub_2$Label,hub_2$group,sep = " ")
# hub_2<-subset(hub_2,hub_sca>0.5)

# tiff("hub.tiff",width = 780,height = 1080)
p1<-ggplot(hub_2) +
  geom_bar(aes(x = hub_sca,y =reorder(rownames(hub_2),hub_sca)),stat = "identity",fill = "#4DAF4A")+
  scale_y_discrete(breaks=reorder(rownames(hub_2),hub_2$hub_sca),
                   labels=reorder(hub_2$Label,hub_2$hub_sca)) +
  labs(x="Hub Sca",y="")+
  theme(axis.text = element_text(size=20),
        axis.title.x = element_text(size = 25))
# dev.off()

p2<-ggplot(hub_2)+
  geom_bar(aes(x = mean,y =reorder(rownames(hub_2),hub_sca)),stat = "identity",fill = "#5DA5DAFF")+
  # scale_y_discrete(breaks=reorder(rownames(hub_2),hub_2$hub_sca),
  #                labels=reorder(hub_2$Label,hub_2$hub_sca)) +
  labs(x="Mean",y="")+
  theme(axis.text = element_text(size=20),
        axis.title.x = element_text(size = 25))

library(ggpubr)

tiff("hub_mean1.tiff",width = 1080,height = 960)
ggarrange(p1,p2+ rremove("y.text"),
          # labels = c("A", "B"),
          ncol = 2, nrow = 1,
          widths = c(1.5,0.7))
dev.off()
