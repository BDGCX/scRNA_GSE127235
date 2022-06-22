library(Seurat)
library(ggplot2)
library(scales)
library(tidyverse)
library(dplyr)
#载入seurat处理后的rds文件
sce <- readRDS(file = "sce_tutorial.rds")

#UMAP
set.seed(1234)#设置随机数种子，这样就不会每次图都变化
sce <- RunUMAP(sce, dims = 1:15, label = T)
head(sce@reductions$umap@cell.embeddings) # 提取UMAP坐标值。
#为了调用ggplot2我们把UMAP的坐标放到metadata中
sce<-AddMetaData(sce,sce@reductions$umap@cell.embeddings,col.name = colnames(sce@reductions$umap@cell.embeddings))
sce@meta.data$new_clusters <- sce@active.ident
#颜色设置
allcolours <- c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00",
            "#808000","#FF00FF","#FA8072","#7B68EE","#9400D3","#800080","#A0522D","#D2B48C","#D2691E","#87CEEB","#40E0D0","#5F9EA0",
            "#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4","#32CD32","#F0E68C","#FFFFE0","#EE82EE",
            "#FF6347","#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887")
allcolour <- hue_pal()(5)
#用ggplot画一个带有标签的umap图

class_avg <- sce@meta.data %>%
  group_by(new_clusters) %>%
  summarise(
    UMAP_1 = median(UMAP_1),
    UMAP_2 = median(UMAP_2)
  )


ggplot(sce@meta.data ,aes(x=UMAP_1,y=UMAP_2))+
  geom_point(aes(color=new_clusters))+
  scale_color_manual(values = allcolour)+
  geom_text(aes(label = new_clusters), data = class_avg)+
  theme(text=element_text(family="Arial",size=18)) +
  theme(panel.background = element_rect(fill='white', colour='black'), 
        panel.grid=element_blank(), 
        axis.title = element_text(color='black',family="Arial",size=18),axis.ticks.length = unit(0.4,"lines"), 
        axis.ticks = element_line(color='black'), 
        axis.ticks.margin = unit(0.6,"lines"),
        axis.line = element_line(colour = "black"), 
        axis.title.x=element_text(colour='black', size=18),
        axis.title.y=element_text(colour='black', size=18),
        axis.text=element_text(colour='black',size=18),
        legend.title=element_blank(),
        legend.text=element_text(family="Arial", size=18),
        legend.key=element_blank())+
  theme(plot.title = element_text(size=22,colour = "black",face = "bold"))  + 
  guides(colour = guide_legend(override.aes = list(size=5)))


                       