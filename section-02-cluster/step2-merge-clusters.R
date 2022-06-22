## 
### ---------------
###
### Create: Jianming Zeng
### Date: 2021-06-26 16:13:22
### Email: jmzeng1314@163.com
### Blog: http://www.bio-info-trainee.com/
### Forum:  http://www.biotrainee.com/thread-1376-1-1.html
### CAFS/SUSTC/Eli Lilly/University of Macau
### Update Log:  2021-06-26  First version
###
### ---------------



rm(list = ls())
library(Seurat)
# devtools::install_github('satijalab/seurat-data')
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
load(file = 'basic.sce.pbmc.Rdata')
levels(Idents(pbmc))

DimPlot(pbmc, reduction = 'umap', 
        label = TRUE, pt.size = 0.5) + NoLegend()

sce=pbmc

# 参考；https://mp.weixin.qq.com/s/9d4X3U38VuDvKmshF2OjHA 
genes_to_check = c('PTPRC', 'CD3D', 'CD3E',  'PRF1' , 'NKG7',
                   'CD19', 'CD79A', 'MS4A1' ,
                   'CD68', 'CD163', 'CD14',
                   "FCER1A", "PPBP")
DotPlot(sce, group.by = 'seurat_clusters',
        features = unique(genes_to_check)) + RotatedAxis()


# 参考：https://mp.weixin.qq.com/s/YYJWbluM86rp9y4CNHBKpQ 
# method 1
Idents(sce)
levels(sce)
head(sce@meta.data)
#  method : 1 
new.cluster.ids <- c("T", "Mono", "T", 
                     "B", "T", "Mono",
                     "T", "DC", "Platelet")#将想要合并的细胞群取相同名字
names(new.cluster.ids) <- levels(sce)
sce <- RenameIdents(sce, new.cluster.ids)
DimPlot(sce, reduction = 'umap', 
        label = TRUE, pt.size = 0.5) + NoLegend()

#  method : 2  
# 一个向量：
cluster2celltype <- c("0"="T",
                  "1"="Mono", 
                  "2"="T", 
                  "3"= "B", 
                  "4"= "T", 
                  "5"= "Mono",
                  "6"= "T", 
                  "7"= "DC", 
                  "8"= "Platelet")
sce[['cell_type']] = unname(cluster2celltype[sce@meta.data$seurat_clusters])
DimPlot(sce, reduction = 'umap', group.by = 'cell_type',
        label = TRUE, pt.size = 0.5) + NoLegend()

#method3
# 一个数据框 
(n=length(unique(sce@meta.data$seurat_clusters)))
celltype=data.frame(ClusterID=0:(n-1),
                    celltype='unkown')
celltype[celltype$ClusterID %in% c(0,2,4,6),2]='T'
celltype[celltype$ClusterID %in% c(3),2]='B'
celltype[celltype$ClusterID %in% c(1,5),2]='Mono' 
celltype[celltype$ClusterID %in% c(7),2]='DC' 
celltype[celltype$ClusterID %in% c(8),2]='Platelet'
sce@meta.data$celltype = "NA"
for(i in 1:nrow(celltype)){
  sce@meta.data[which(sce@meta.data$seurat_clusters == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}
table(sce@meta.data$celltype)
phe=sce@meta.data
DimPlot(sce, reduction = 'umap', group.by = 'celltype',
        label = TRUE, pt.size = 0.5) + NoLegend()

head(sce@meta.data)
table(sce@meta.data$cell_type,
      sce@meta.data$celltype)

p1=DimPlot(sce, reduction = 'umap', group.by = 'celltype',
           label = TRUE, pt.size = 0.5) + NoLegend()
p2=DotPlot(sce, group.by = 'celltype',
           features = unique(genes_to_check)) + RotatedAxis()
library(patchwork)
p1+p2


