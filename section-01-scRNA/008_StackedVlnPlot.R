#方法一，通过Seuart->scanpy来实现，现在问题就变成为Seurat对象到scanpy对象的转换
#读入seurat处理后的rds文件，以下是R代码
rm(list = ls())
options(stringsAsFactors = F)

library(Seurat)
#install.packages("hdf5r") loomR安装需要hdf5r
#devtools::install_github(repo = "mojaveazure/loomR", ref = "develop")
library(loomR)
library(SeuratDisk)
library(ggplot2)
library(scales)
library(dplyr) 
library(scater)
library(reticulate)
#设定python环境
use_condaenv("C:/Users/bdgcx/anaconda3")
py_config()#安装的python版本环境查看，显示anaconda和numpy的详细信息。
py_available()#[1] TRUE   #检查您的系统是否安装过Python
py_module_available("pandas")#检查“pandas”是否安装
py_module_available("loompy")

## 在R中载入python模块
scanpy <-  import("scanpy")
loompy <-  import("loompy")
pandas <- import("pandas")
numpy <-  import("numpy")
#载入seurat处理后的rds文件
sce <- readRDS(file = "sce_tutorial.rds")

#将sce@meta.data$new_clusters转换为字符串以便使用python进行绘图时分组
sce@meta.data$new_clusters <- as.character(sce@meta.data$new_clusters)
# seurat对象转换为loom文件,使用SeuratDisk包
sce.loom <- as.loom(x = sce, filename = "C:/Users/bdgcx/OneDrive/single cell/GSE127235_DKD1_mice/sce.loom", verbose = FALSE)
# Always remember to close loom files when done
sce.loom$close_all()

#直接在R中设定python环境
library(reticulate)
#设定python环境
use_condaenv("C:/Users/bdgcx/anaconda3")
py_config()#安装的python版本环境查看，显示anaconda和numpy的详细信息。
py_available()#[1] TRUE   #检查您的系统是否安装过Python
py_module_available("pandas")#检查“pandas”是否安装
py_module_available("loompy")

## 在R中载入python模块,未探索明白，遂弃
scanpy <-  import("scanpy")
loompy <-  import("loompy")

##接下来使用python
#Scanpy读取loom文件转换为能够操作的anndata对象，以下是python代码
#具体代码见python文件夹中的stackedvlnplot.py文件

#方法二，，利用使用基于scanpy包衍生的scanyuan包来说实现，数据处理及流程同方法一

#方法三，利用R函数实现单细胞StackedVlnPlot
library(Seurat)
library(ggplot2)
sce@meta.data$new_clusters <- sce@active.ident
markers <- c("Cd300lg","Egfl7","Rasgrp3","Ptprb","Cyyr1", "Gpihbp1", "Ramp3",
             "Tspan7","Adgrl4","Fabp4","Sfrp2","Agtr1a","Mgp","Ptn", 
             "Hopx","S1pr3","Cd248","P2rx1","Lhfp","Fhl2","Clic3","Cdkn1c",
             "Enpep","Dpp4","Plce1","Magi2","Rhpn1","Srgap1","Ildr2",
             "Robo2","Mgat5")
#1.直接使用VlnPlot函数的stack = T
#使用CaseMatch去除不存在的基因
library(limma)
markers <- CaseMatch(markers,rownames(sce))
markers <- as.character(markers)
idents <- c("TC","IC","POD","MC","EC")
colpalette <- c("#F8766D" ,"#A3A500", "#00BF7D" ,"#00B0F6", "#E76BF3")
VlnPlot(sce,features = markers,
        pt.size = 0,group.by = "new_clusters",
        cols = colpalette,
        fill.by="ident",#小提琴颜色按照分组进行划分
        sort = F,
        idents=idents,
        #flip = T,转换x轴和y轴的位置
        stack = T)+NoLegend()
ggsave("svp.png",width = 3000,height = 1700,dpi = 300,units = "px")

#2.别人写好的函数(推荐此方法，与文献中绘图一致)
###载入需要的R包
library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(tibble)
##load测试数据
sce <- readRDS(file = "sce_tutorial.rds")

markers <- c("Cd300lg","Egfl7","Rasgrp3","Ptprb","Cyyr1", "Gpihbp1", "Ramp3",
             "Tspan7","Adgrl4","Fabp4","Sfrp2","Agtr1a","Mgp","Ptn", 
             "Hopx","S1pr3","Cd248","P2rx1","Lhfp","Fhl2","Clic3","Cdkn1c",
             "Enpep","Dpp4","Plce1","Magi2","Rhpn1","Srgap1","Ildr2",
             "Robo2","Mgat5")

colpalette <- c("#F8766D" ,"#A3A500", "#00BF7D" ,"#00B0F6", "#E76BF3")

##第一个函数从Seurat对象中获取表达值并转化为需要的格式tidy
##tibble::rownames_to_column将行名转换为第一列
gotData<-function(seurat.obj,features,groups){
  mat<-GetAssayData(seurat.obj,assay = "RNA",slot = "data")[features,]
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

##第二个函数画图
plot_stacked_violin<-function(plotData,xlab,ylab,cols){
  ggplot(plotData,aes(y=xaxis,x=value,fill=cellID))+
    geom_violin()+
    theme_test()+
    facet_grid(vars(gene), vars(cellID),#多图panel组合，第一和参数行分，第二个按列分
    switch = "y")+#右边的标签放左边
    xlab(xlab)+
    ylab(ylab)+
    scale_fill_manual(values = cols)+
    theme(panel.spacing.y =unit(0.13,"cm"),#各个面板之间的距离panel.spacing.x/y可以不同方向分别设置
          panel.spacing.x =unit(0.3,"cm"),
          strip.background = element_rect(fill="transparent",color = "white"),
          axis.text.x = element_blank(),
          #axis.ticks.x = element_blank(),#轴须标签，刻度线
          axis.ticks.y = element_blank(),
          panel.border = element_rect(size=0.7,colour = "black"),
          strip.text.x = element_text(size=15,face = "bold"),
          strip.text.y = element_text(angle = 0,size = 13),
          axis.text.y = element_blank(),
          axis.title.x.bottom = element_text(size = 13),
          axis.title.y = element_text(size = 13))+
    NoLegend()
}

plot_stacked_violin_1<-function(plotData,xlab,ylab,cols){
  ggplot(plotData,aes(y=xaxis,x=value,fill=cellID))+
    geom_violin()+
    facet_wrap(gene~cellID)+#多图panel组合，只有一个因素facet_wrap(.~cellID)，多个因素facet_wrap(~cellID+gene),gene~cellID代表行列
    theme_test()+
    xlab(xlab)+
    ylab(ylab)+
    scale_fill_manual(values = cols)+
    theme(panel.spacing=unit(0.05,"cm"),
          strip.background = element_rect(fill="transparent",color = "white"),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.border = element_rect(size=0.7,colour = "black"),
          strip.text = element_text(size=10,face = "italic"),
          axis.text.y = element_text(size = 6,face="bold"),
          axis.title.x.bottom = element_text(size = 8),
          axis.title.y = element_text(size = 13))+
    NoLegend()
}

plotData <- gotData(sce,markers,"new_clusters")#"Cyyr1", "Gpihbp1", "Ramp3"),"new_clusters")
                    
plotData$xaxis <- c(rep(1,times=nrow(plotData)))

## 调整因子水平
plotData$gene <- factor(plotData$gene,levels=c("Cd300lg","Egfl7","Rasgrp3","Ptprb","Cyyr1", "Gpihbp1", "Ramp3",
                                               "Tspan7","Adgrl4","Fabp4","Sfrp2","Agtr1a","Mgp","Ptn", 
                                               "Hopx","S1pr3","Cd248","P2rx1","Lhfp","Fhl2","Clic3","Cdkn1c",
                                               "Enpep","Dpp4","Plce1","Magi2","Rhpn1","Srgap1","Ildr2",
                                               "Robo2","Mgat5"))

plot_stacked_violin(plotData,"Expression (log-scale)"," ",colpalette)#空白处是y轴标签

ggsave("Stackedviolin-ggplot2.png",width = 2500,height = 4000,dpi = 300,units = "px")

#3.使用别人写好的函数
#具体用法见https://divingintogeneticsandgenomics.rbind.io/post/stacked-violin-plot-for-visualizing-single-cell-data-in-seurat/
modify_vlnplot <- function(obj, feature, pt.size = 0, plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),...) {
  p <- VlnPlot(obj, feature = feature, pt.size = pt.size, ...
               ) +
    xlab("") + ylab(feature) + ggtitle("") +
    theme(legend.position = "none",
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_line(),
          axis.title.y = element_text(size = rel(1), angle = 0, vjust = 0.5),
          plot.margin = plot.margin )
  return(p)
}
library(patchwork)
## main function
StackedVlnPlot <- function(obj, features, pt.size = 0, plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"), ...) {
  plot_list <- purrr::map(features, function(x) modify_vlnplot(obj = obj,feature = x, ...))
  plot_list[[length(plot_list)]]<- plot_list[[length(plot_list)]] +
    theme(axis.text.x=element_text(), axis.ticks.x = element_line())
  p <- patchwork::wrap_plots(plotlist = plot_list, ncol = 1)
  return(p)
}

#配色方案
library(ggsci)
colpalettes<-unique(c(pal_npg("nrc")(10),pal_aaas("default")(10),pal_nejm("default")(8),pal_lancet("lanonc")(9),
                      pal_jama("default")(7),pal_jco("default")(10),pal_ucscgb("default")(26),pal_d3("category10")(10),
                      pal_locuszoom("default")(7),pal_igv("default")(51),
                      pal_uchicago("default")(9),pal_startrek("uniform")(7),
                      pal_tron("legacy")(7),pal_futurama("planetexpress")(12),pal_rickandmorty("schwifty")(12),
                      pal_simpsons("springfield")(16),pal_gsea("default")(12)))
#绘图
markers <- c("Cd300lg","Egfl7","Rasgrp3","Ptprb","Cyyr1", "Gpihbp1", "Ramp3",
             "Tspan7","Adgrl4","Fabp4","Sfrp2","Agtr1a","Mgp","Ptn", 
             "Hopx","S1pr3","Cd248","P2rx1","Lhfp","Fhl2","Clic3","Cdkn1c",
             "Enpep","Dpp4","Plce1","Magi2","Rhpn1","Srgap1","Ildr2",
             "Robo2","Mgat5")
markers <-  c("Cd300lg","Egfl7","Rasgrp3","Ptprb","Cyyr1")
StackedVlnPlot(sce, markers, 
               pt.size=0, 
               cols=colpalettes,idents=c("EC","MC","POD"))#idents选择需要绘图的组别
  +theme(panel.background = element_rect(fill='white', colour='black'), 
        panel.grid=element_blank()) 
