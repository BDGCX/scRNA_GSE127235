rm(list = ls())
options(stringsAsFactors = F)
#载入包
library(Seurat)
library(ggplot2)
library(ggsignif)
library(ggpubr)
library(RColorBrewer)
library(dplyr)
library(tidyr)
library(tibble)

#读入seurat处理后的rds文件
sce <- readRDS(file = "sce_tutorial.rds")

markers <- c("Nphs2","Kdr","Cd74","Pdgfrb","Atp1b1")
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

dat <- plotData1[plotData1$gene=="Kdr"&plotData1$treatment=="Control",]
#绘图
p <- ggplot(data=dat)+ 
  geom_boxplot(mapping=aes(x=cellID,y=value,colour = cellID ), #箱线图
               alpha = 1,
               size=0.2,
               width = 0.4,
               outlier.size=0.6,
               outlier.shape = 16)+ 
  geom_jitter(mapping=aes(x=cellID,y=value,colour = cellID), #散点
              alpha = 0.3,size=0.6)+
  #scale_color_manual(limits=c("A","B","C","D"), 
                     #values=c("#85B22E","#5F80B4","#E29827","#922927"))+ #颜色
  geom_signif(mapping=aes(x=cellID,y=value), # 不同组别的显著性
              comparisons = list(c("EC", "MC"), # 哪些组进行比较
                                 c("EC", "TC"),
                                 c("MC", "IC")),
              map_signif_level=T, # T显示显著性，F显示p value
              tip_length=c(0,0,0,0,0,0), # 修改显著性线两端的长短
              y_position = c(4.9,5.5,3.4), # 设置显著性线的位置高度
              size=0.2, # 修改线的粗细
              textsize = 1.5, # 修改显著性标记的大小
              test = "t.test")+ # 检验的类型
  theme_classic(  # 主题设置，这个是无线条主题
    base_line_size = 0.4)+ # 坐标轴的粗细
  labs(title=" ",x="",y="Expression")+ # 添加标题，x轴，y轴内容
  theme(plot.title = element_text(size = 6,
                                  colour = "black",
                                  hjust = 0.5),
        axis.title.y = element_text(size = 6, 
                                    # family = "myFont", 
                                    color = "black",
                                    #face = "bold", 
                                    vjust = 1.9, 
                                    hjust = 0.5, 
                                    angle = 90),
        legend.title = element_text(color="black", # 修改图例的标题
                                    size=6),
        legend.key.size = unit(8,"pt"),#调整图例大小
        legend.text = element_text(color="black", # 设置图例标签文字
                                   size = 4), 
        axis.text.x = element_text(size = 6,  # 修改X轴上字体大小，
                                   color = "black", # 颜色
                                   #face = "bold", #  face取值：plain普通，bold加粗，italic斜体，bold.italic斜体加粗
                                   vjust = 0.5, # 位置
                                   hjust = 0.5, 
                                   angle = 0), #角度
        axis.text.y = element_text(size = 6,  
                                   color = "black",
                                   #face = "bold", 
                                   vjust = 0.5, 
                                   hjust = 0.5, 
                                   angle = 0) 
  )
ggsave("./figure/signif_barplot.png",p,width = 1500,height = 1000,
       dpi = 400,units = "px")
