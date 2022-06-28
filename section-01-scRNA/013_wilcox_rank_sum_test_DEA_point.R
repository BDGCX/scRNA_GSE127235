rm(list = ls())
options(stringsAsFactors = F)
#options(future.globals.maxSize = 1000 * 1024^2)

library(dplyr)
library(tibble)
library(data.table)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(Seurat)
#载入数据
#方法1.使用FindMarkers函数的wilcoxon rank sum test
sce <- readRDS("sce_tutorial.rds")
#从Seurat对象中选择细胞群cluster
deg <- FindMarkers(object = sce, 
                   ident.1 = "Diabetes",#病例组。ident.2=NULL时默认剩下所有细胞为对照组
                   ident.2 = "Control",
                   test.use='wilcox',
                   group.by = "treatment",
                   subset.ident="EC")
DT::datatable(deg)
write.csv(deg,file="deg.csv")

deg$gene <- rownames(deg)
#火山图使用avg_log2FC和pct.1-pct.2的绝对值绘制离群点比较不离谱
ggplot(data = deg,aes(x=pct.1,y=pct.2))+#这里可能需要用abs(avg_log2FC)和-log10(p_val)画
  geom_point(color="#888888",size=1.3)+
  xlab("Diabetes")+
  ylab("Control")+
  theme(axis.text.x = element_text(colour = "black"))+
  theme_light()+
  geom_text_repel(data=subset(deg, avg_log2FC> 1.5| avg_log2FC < -1), 
                  aes(label=gene),
                  col="black",
                  alpha = 0.8,
                  size=3)+
  geom_point(data=subset(deg,deg$p_val<0.05 & deg$avg_log2FC < -1),color="blue") +
  geom_point(data=subset(deg,deg$p_val<0.05 & deg$avg_log2FC > 1.5),color="red") 

ggsave("./figure/DEA_point.png",width = 2500,height = 1600,dpi = 400,units = "px")
#方法2.各个单细胞亚群独立在两个分组做差异分析
table(sce$new_clusters,sce$treatment )#均为meta.data中的分组信息
# 如下所示，两个分组，每个分组里面的都是有这些不同单细胞亚群
#Control Diabetes
#EC      163      216
#MC      106       41
#POD      56       24
#IC       10       58
#TC       11       19

#进行批量差异分析
Idents(sce) <-  paste0('c',sce$treatment )
table(Idents(sce))
degs <-  lapply(unique(sce$new_clusters), function(x){
  FindMarkers(sce[,sce$new_clusters==x],ident.1 = 'cDiabetes',
              ident.2 = 'cControl')
})
EC <- degs[[1]]
MC <- degs[[2]]
do.call(rbind,lapply(degs, function(x){
  table(x$avg_log2FC > 0 )
}))


#方法3.使用stats包的wilcoxon rank sum test
counts <- sce@assays$RNA@counts
counts[1:4,1:4]
exprset <- counts
exprset <- as.matrix(exprset)

rm(sce)

#判断样本归类
a <- unlist(lapply(colnames(exprset),function(x){
  strsplit(as.character(x),"_")[[1]][4]
})
)
group_list <- ifelse(a =="Control",'control','diabetes')
table(group_list)

conNum <- 346#control组织样品数目
treaNum <- 358#diabetes组织样品数目

grp <- c(rep(1,conNum),rep(2,treaNum))#将对照组织变为1，DM组织为2

outTab <- data.frame()

for (i in rownames(exprset)) {
  geneName=rownames(exprset)
  
  rt=rbind(expression=exprset[i,],group=grp)
  rt=as.matrix(t(rt))
  
  wilcoxTest<-wilcox.test(expression~group,data=rt)  
  
  conGeneMean=mean(exprset[i,1:conNum])
  treatGeneMean=mean(exprset[i,(conNum+1):ncol(exprset)])
  
  log2FC=log2(treatGeneMean)-log2(conGeneMean)#log2(treatGeneMean/conGeneMean)
  
  pValue=wilcoxTest$p.value
  
  conMed=median(exprset[i,1:conNum])
  treatMed=median(exprset[i,(conNum+1):ncol(exprset)])
  
  diffMed=treatMed-conMed
  
  if(((log2FC>0) & (diffMed>0)) | ((log2FC<0) & (diffMed<0)) ){
    outTab=rbind(outTab,cbind(gene=i,conMean=conGeneMean,treaMean=treatGeneMean,log2FC=log2FC,pValue=pValue
    )) }
  
}

outTab$adj.P.Val <- p.adjust(as.numeric(outTab$pValue), method = "BH")

outTab$FDR=p.adjust(as.numeric(as.vector(outTab$pValue)),method = "fdr")


for (i in rownames(exprset)) {
  geneName=rownames(exprset)
  
  rt=rbind(expression=exprset[i,],group=grp)
  rt=as.matrix(t(rt))
  
  wilcoxTest<-wilcox.test(expression~group,data=rt)  
  
  conGeneMean=mean(exprset[i,1:conNum])
  treatGeneMean=mean(exprset[i,(conNum+1):ncol(exprset)])
  
  log2FC=log2(treatGeneMean)-log2(conGeneMean)#log2(treatGeneMean/conGeneMean)
  
  pValue=wilcoxTest$p.value
  
  conMed=median(exprset[i,1:conNum])
  treatMed=median(exprset[i,(conNum+1):ncol(exprset)])
  
  diffMed=treatMed-conMed
  outTab=rbind(outTab,cbind(gene=i,conMean=conGeneMean,treaMean=treatGeneMean,log2FC=log2FC,pValue=pValue))
                              
}

DEG <- outTab[which(abs(as.numeric(outTab$log2FC))>=1),]

DEG$adj.P.Val <- p.adjust(as.numeric(DEG$pValue), method = "BH")

write.csv(DEG,file = "w_DEA.csv")

DEG <- read.csv("w_DEA.csv",header = T)
DEG <- DEG[,-1]
rownames(DEG) <- DEG[,1]
DEG <- DEG[,-1]

