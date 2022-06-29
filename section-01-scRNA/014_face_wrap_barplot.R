rm(list = ls())
options(stringsAsFactors = F)
#载入包
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyr)
library(tibble)

#读入seurat处理后的rds文件
sce <- readRDS(file = "sce_tutorial.rds")

markers <- c("Lrg1","Vegf","Ctgf","Tnf")
#使用CaseMatch去除不存在的基因
library(limma)
markers <- CaseMatch(markers,rownames(sce))

# 提取表达矩阵
gotData<-function(seurat.obj,features,groups){
  mat<-GetAssayData(seurat.obj,assay = "RNA",slot = "data")[features,]#注意这里是从data中提取，而不是counts
  plotData<-mat%>%
    as.data.frame()%>%
    tibble::rownames_to_column(var="gene")%>%
    as_tibble()%>%
    tidyr::pivot_longer(names_to = "cell",values_to="exp",cols=2:(ncol(mat)+1))
  cellmeta<-seurat.obj@meta.data%>%
    tibble::rownames_to_column(var="cell")%>%
    as_tibble()%>%
    select(cell,sym(groups))
  plotData<-plotData%>%
    left_join(cellmeta,by="cell")%>%
    setNames(c("gene","cell","value","cellID"))
  plotData
}

plotData <- gotData(sce,markers,"new_clusters")

#进一步使用处理分组treatment进行图像绘制
#提取ac里的处理分组信息，Control或Diabetes
trt <- unlist(lapply(plotData$cell, function(x){
  strsplit(as.character(x),'_')[[1]][4]
}))
plotData$treatment <- trt
plotData1 <- plotData[order(plotData$gene,plotData$cellID,plotData$treatment),]


ggplot(plotData1,aes(x=cell,y=value,fill=cellID))+
  geom_bar(stat = "identity")+
  theme_test()+
  facet_grid(gene~cellID+treatment,
             scales = "free",
             space = "free",
             switch = "y"
             )+#多图panel组合，第一个参数按行分，第二个按列分
  #xlab(xlab)+
  ylab("Expression(log-scale)")+
  #scale_fill_manual(values = cols)+
  theme(panel.spacing.y =unit(0.13,"cm"),#各个面板之间的距离panel.spacing.x/y可以不同方向分别设置
        panel.spacing.x =unit(0.13,"cm"),
        strip.background.x = element_rect(fill=c("#8787FFFF","#6EE2FFFF"),colour = c("#8787FFFF","#6EE2FFFF"),size =0.5 ),
        strip.background.y = element_rect(fill = "transparent",colour = "white"),
        axis.text.x = element_blank(),
        #axis.ticks.x = element_blank(),#轴须标签，刻度线
        axis.ticks.y = element_blank(),
        panel.border = element_rect(size=0.6,colour = "black"),
        #strip.text.x  = element_text(margin =margin(0,0,0,0) ),
        strip.text.x  =element_blank(),
        strip.text.y = element_text(angle = 0,size = 8),
        axis.text.y = element_blank(),
        axis.title.x.bottom = element_blank(),
        legend.key.size = unit(15,"pt") ,#调整图例大小
        axis.title.y = element_text(size = 10)) 
 
ggsave("barplot.pdf",width = 3000,height = 1500,dpi = 400,units = "px")

#上方分组使用热图的聚类分组用AI组合