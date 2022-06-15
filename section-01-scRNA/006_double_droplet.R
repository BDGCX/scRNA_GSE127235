load("../sce_tutorial.Rdata")
library(Seurat)
##1、detect special cells----
#empty droplet
#BiocManager::install("DropletUtils")
library(scDblFinder)
library(DropletUtils)
e.out <- emptyDrops(GetAssayData(sce,slot="counts",assay="RNA"))
#Error in testEmptyDrops(m, lower = lower, ...) : 
#no counts available to estimate the ambient profile
##https://support.bioconductor.org/p/123554/#123562
#如上回答所说，empty droplet往往在第一步就已经过滤掉了，而一般上传到GEO的也都是过滤掉空液滴的。

#double droplet
#https://osca.bioconductor.org/doublet-detection.html
#BiocManager::install("scran")
head(sce@meta.data)
library(scran)
#GetAssayData(scRNA,slot="counts",assay="RNA")[1:8,1:4]
?doubletCluster #检查有无double droplet聚在一起的类
db.test <- doubletCluster(GetAssayData(sce,slot="counts",assay="RNA"),
                          clusters=sce@meta.data$seurat_clusters)
head(db.test)
table(sce@meta.data$seurat_clusters)
library(scater)
chosen.doublet <- rownames(db.test)[isOutlier(db.test$N, 
                                              type="lower", log=TRUE)]
chosen.doublet #结果显示没有
#还有其它多种方法