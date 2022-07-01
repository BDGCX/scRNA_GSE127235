rm(list = ls())
options(stringsAsFactors = F)
#载入包
library(Seurat)
library(ggplot2)
library(tidyverse)
library(stringr)
library(tibble)
#读入seurat处理后的rds文件
sce <- readRDS(file = "sce_tutorial.rds")
sce<- subset(sce,idents=c("EC","MC","POD"))

#提取表达矩阵
markers <- c("Plat","Tmeme108","Col18a1", "Micul", "Dag1",
             "Loxl2", "Fgfr1", "Lgr4" ,"Bmp7", "Ncam1", "Npnt","Lypd1",
             "Vegfa", "Adam9", "Cyr61", "Sdc1", "Itga4", "Mmp2", "Pdgfrb",
             "Art3","Angpt2", "Ngf","Nt5e", "Pdgfra", "Ptn", "Notch3",
             "Lrp1", "Mfge8", "Apoe", "Cacna2d1", "Epha3", 'Nrp1', "B2m",
             "Lgals3bp", "Ephb4", "Ramp2", "Fgfr3", "Plaur", "Efna1",
             "Egfl7", "Cdh5", "F8", "Sema6a", "Gas6", "Pdia5")

#使用CaseMatch去除不存在的基因
library(limma)
markers <- CaseMatch(markers,rownames(sce))

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
#将长数据转化为宽数据
cts<-reshape2::dcast(plotData,gene~cell+cellID,value.var="value")

dat01 <- cts
#计算z-scores
a <- dat01 %>% 
  rowwise() %>% #rowwise按行进行数据处理
  mutate(mean_value = mean(c_across(2:607)),
         sd_value = sd(c_across(2:607))) %>% 
  mutate(across(2:607,~(.x-mean_value)/sd_value)) %>% 
  select(-c(mean_value,sd_value)) -> dat01.2
#宽格式转换为长格式
c <- unlist(lapply(colnames(dat01.2)[2:607], function(x){
  strsplit(as.character(x),'_')[[1]][1:4]
}))
dat01.2 %>% 
  reshape2::melt(id.vars="gene") %>% 
  mutate(new_var=str_replace_all(variable,c("EXP[12]_"="",
                                        "COL[0123456789][0123456789]_"="",
                                        "ROW[0123456789][0123456789]_"="",
                                        "Control_"="",
                                        "Diabetes_"=""))) %>% 
  group_by(gene,new_var) %>% 
  summarise(mean_value=mean(value)) %>%
  ungroup() -> dat01.3

#调整因子水平
dat01.3$new_var<-factor(dat01.3$new_var,
                        levels = c("EC","MC","POD"))
dat01.3$gene <-factor(dat01.3$gene,
                      levels = as.character(rev(markers)))#rev将向量倒序
dat <- dat01.3
dat$mean_value[dat$mean_value>=1]=1
#绘制热图
#library(paletteer)
ggplot(data = dat,
       aes(x=new_var,y=gene))+
  geom_tile(aes(fill=mean_value),
            color="white")+
  #scale_fill_paletteer_c("basetheme::ink",
                         #direction = 1,
                         #name="Expression level (Z-score)",
                         #limits=c(-2,2))+
  scale_fill_gradient2(low = "blue",
                       mid = "white",
                       high ="red",
                       midpoint = 0,
                       guide = "colourbar"
                       )+
  scale_x_discrete(position = "top")+
  scale_y_discrete(position = "right")+
  labs(x=NULL,y=NULL)+
  theme_minimal()+
  ggtitle("Ligand")+
  theme(panel.grid = element_blank(),
        plot.title = element_text(size = 10,face = "bold",hjust = 0.5),#标题居中
        legend.position = "left",
        legend.key.width = unit(8,units = "pt"),
        legend.title = element_blank(),
        axis.text.x = element_text(angle = 0,
                                   hjust = 0.5,
                                   vjust=1),
        )+
  guides(color=guide_colorbar(direction ="vertical",
                              ticks = FALSE))
                        
ggsave("./figure/z-score.png",width = 1000,height = 2500,
       dpi = 400,units = "px")
