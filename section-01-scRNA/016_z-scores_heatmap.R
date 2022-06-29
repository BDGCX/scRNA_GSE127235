rm(list = ls())
options(stringsAsFactors = F)
#载入包
library(Seurat)
library(ggplot2)
library(tidyverse)
library(stringr)

#读入seurat处理后的rds文件
sce <- readRDS(file = "sce_tutorial.rds")
sce<- subset(sce,idents=c("EC","MC","POD"))
cts <- GetAssayData(sce,slot = "data")
dat01 <- as.matrix(cts)
#计算z-scores
dat01 %>% 
  rowwise() %>% #rowwise按行进行数据处理
  mutate(mean_value = mean(c_across(2:16)),
         sd_value = sd(c_across(2:16))) %>% 
  mutate(across(2:16,~(.x-mean_value)/sd_value)) %>% 
  select(-c(mean_value,sd_value)) -> dat01.2
#宽格式转换为长格式
dat01.2 %>% 
  reshape2::melt(id.vars="Gene") %>% 
  mutate(new_var=str_replace(variable,'-[123]','')) %>% 
  group_by(Gene,new_var) %>% 
  summarise(mean_value=mean(value)) %>%
  ungroup() -> dat01.3




#绘制热图
library(ggplot2)  
library(paletteer)

dat01.3$new_var<-factor(dat01.3$new_var,
                        levels = c("Ph","Sb","Xy","Pi","Le1"))

ggplot(data = dat01.3,
       aes(x=Gene,y=new_var))+
  geom_tile(aes(fill=mean_value),
            color="white")+
  scale_fill_paletteer_c("ggthemes::Classic Red-Green",
                         direction = -1,
                         name="Expression level (Z-score)",
                         limits=c(-2,2))+
  scale_y_discrete(position = "right")+
  labs(x=NULL,y=NULL)+
  theme_minimal()+
  theme(panel.grid = element_blank(),
        legend.position = "top",
        axis.text.x = element_text(angle = 60,
                                   hjust = 1,
                                   vjust=1),
        plot.margin = unit(c(0.2,0.2,0.2,1),'cm'))+
  guides(fill=guide_colorbar(title.position = "top",
                             title.hjust = 0.5,
                             barwidth = 10,
                             barheight = 0.5,
                             ticks = FALSE))

