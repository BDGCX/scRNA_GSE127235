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
# 首先提取T细胞子集
DimPlot(pbmc, reduction = 'umap', 
        label = TRUE, pt.size = 0.5) + NoLegend()
FeaturePlot(pbmc, features = c("CD3D","CD3E"))
sce=pbmc
table(Idents(sce))
t_sce = sce[, Idents(sce) %in% 
                 c( "Naive CD4 T" ,  "Memory CD4 T" ,
                    'CD8 T','NK')]
# 然后进行标准的降维聚类分群
# 代码不要变动
sce=t_sce
sce <- NormalizeData(sce, normalization.method = "LogNormalize",
                     scale.factor = 1e4) 
sce <- FindVariableFeatures(sce, selection.method = 'vst',
                            nfeatures = 2000)
sce <- ScaleData(sce, vars.to.regress = "percent.mt")
sce <- RunPCA(sce, features = VariableFeatures(object = sce)) 
sce <- FindNeighbors(sce, dims = 1:10)
sce <- FindClusters(sce, resolution = 1 )
# Look at cluster IDs of the first 5 cells
head(Idents(sce), 5)
table(sce$seurat_clusters) 
sce <- RunUMAP(sce, dims = 1:10)
DimPlot(sce, reduction = 'umap')

sce.markers <- FindAllMarkers(object = sce, only.pos = TRUE, 
                              min.pct = 0.25, 
                              thresh.use = 0.25)
DT::datatable(sce.markers)
# write.csv(sce.markers,file=paste0(pro,'_sce.markers.csv'))
 
library(dplyr) 
# 在seurat V4里面， 是 avg_log2FC ， 但是如果是V3版本包，仍然是 avg_logFC
top10 <- sce.markers %>% group_by(cluster) %>% top_n(5, avg_logFC)
DoHeatmap(sce,top10$gene,size=3)  
p <- DotPlot(sce, features = unique(top10$gene),
             assay='RNA' )  + coord_flip()

p
ggsave('check2-top5-for-all-Tcells.pdf',height = 18)

# 然后查看文献的标记基因
# 2021年5月的文章：《A single-cell map of intratumoral changes during anti-PD1 treatment of patients with breast cancer
# 参考：https://mp.weixin.qq.com/s/FhEASxwTF8e7lNoKxTuxrg 

marker_genes= c("LEF1","TCF7","SELL","IL7R","CD40LG","ANXA1","FOS",
                "JUN","FOXP3","SAT1","IL2RA","CTLA4","PDCD1","CXCL13","CD200",
                "TNFRSF18","CCR7","NELL2","CD55","KLF2","TOB1","ZNF683","CCL5",
                "GZMK","EOMES","ITM2C","CX3CR1","GNLY","GZMH","GZMB","LAG3","CCL4L2",
                "FCGR3A","FGFBP2","TYROBP","AREG","XCL1","KLRC1","TRDV2","TRGV9","MTRNR2L8",
                "KLRD1","TRDV1","KLRC3",
                "CTSW","CD7","MKI67","STMN1","TUBA1B","HIST1H4C" )


p <- DotPlot(sce, features = marker_genes,
             assay='RNA'  )  + coord_flip() +
  theme(axis.text.x = element_text(angle = 90))

p
ggsave('check_g1_markers_by_Tcell-SubType.pdf')


marker_genes =c("CD3D","CD4","CD8A","CCR7","LEF1","SELL","TCF7","GNLY","IFNG","NKG7","PRF1",
                "GZMA","GZMB","GZMH","GZMK","HAVCR2","LAG3","PDCD1","CTLA4","TIGIT","BTLA","KLRC1",
                "ANXA1","ANKRD28","CD69","CD40LG","ZNF683","FOXP3","IL2RA","IKZF2","NCR1","NCAM1",
                "TYROBP","KLRD1","KLRF1","KLRB1","CX3CR1","FCGR3A",
                "XCL1","XCL2","TRDV2","TRGV9","TRGC2","MKI67","TOP2A")

p <- DotPlot(sce, features = marker_genes,
             assay='RNA'  )  + coord_flip() +
  theme(axis.text.x = element_text(angle = 90))

p
ggsave('check_g2_markers_by_Tcell-SubType.pdf')

# naive (LEF1, SELL, TCF7),
# effector (IFNG), 
# cytotoxicity (GZMB, PRF1), 
# early and general exhaustion (PDCD1, CTLA4, ENTPD1 ) .
# antigen presentation (CD74, HLA-DRB1/5, HLA-DQA2)

genes_to_check = c('PTPRC', 'CD3D', 'CD3E', 'FOXP3',
                   'CD4','IL7R','NKG7','CD8A',
                   'LEF1', 'SELL', 'TCF7', # naive marker
                   'IFNG','GZMB', 'PRF1',
                   'PDCD1', 'CTLA4', 'ENTPD1'  )
p <- DotPlot(sce, features = genes_to_check,
             assay='RNA'  )  + coord_flip() +
  theme(axis.text.x = element_text(angle = 90))

p








