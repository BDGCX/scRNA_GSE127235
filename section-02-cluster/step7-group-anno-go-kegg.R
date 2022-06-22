rm(list=ls())
options(stringsAsFactors = F)
#加载包
library(Seurat)
library(gplots)
library(ggplot2)
library(clusterProfiler)
library(org.Mm.eg.db)
#使用sce.markers文件
sce <- readRDS("sce_tutorial.rds")
sce.markers <- FindAllMarkers(object = sce, only.pos = TRUE, min.pct = 0.25, 
                              thresh.use = 0.25)
#id转换
ids=bitr(sce.markers$gene,'SYMBOL','ENTREZID','org.Mm.eg.db')
sce.markers <- merge(sce.markers,ids,by.x='gene',by.y='SYMBOL')
gcSample <- split(sce.markers$ENTREZID, sce.markers$cluster)
gcSample # entrez id , compareCluster 

xx <- compareCluster(gcSample, fun="enrichGO",
                     OrgDb="org.Mm.eg.db", pvalueCutoff=0.05)
#KEGG下载物种名报错
install.packages('R.utils')
R.utils::setOption( "clusterProfiler.download.method",'auto' )

xx <- compareCluster(gcSample, fun="enrichKEGG",
                     organism="mmu", pvalueCutoff=0.05)
#organism人为hsa,小鼠为mmu,大鼠为rno
p <- dotplot(xx) 
p+ theme(axis.text.x = element_text(angle = 45, 
                                    vjust = 0.5, hjust=0.5))
p

# msigdb， 单细胞相关基因集
source('code/group_kegg_go.R')
pro='all'
group_kegg_go(gcSample, pro )

