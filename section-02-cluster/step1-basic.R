##见code →section-01-scRNA

#以下为生信技能树演示代码
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

# Load the PBMC dataset
# https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz
pbmc.data <- Read10X(data.dir = "filtered_gene_bc_matrices/hg19/")
# Initialize the Seurat object with the raw (non-normalized data).
#barcode是细胞名字，gene是基因名字，matrix是表达量矩阵
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k",
                           min.cells = 3, min.features = 200)
pbmc

pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
#Visualize QC metrics as a violin plot
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# FeatureScatter is typically used to visualize feature-feature relationships, but can be used for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt") 
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") 
plot1 + plot2
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 1e4) 
pbmc <- FindVariableFeatures(pbmc, selection.method = 'vst', nfeatures = 2000)
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)
# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2#图太大显示不出来报错，可以直接储存到pdf或图片格式
pbmc <- ScaleData(pbmc, vars.to.regress = "percent.mt")
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))  
# 
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)
# Look at cluster IDs of the first 5 cells
head(Idents(pbmc), 5)
table(pbmc$seurat_clusters) 
pbmc <- RunUMAP(pbmc, dims = 1:10)
# note that you can set `label = TRUE` or use the LabelClusters function to help label individual clusters
DimPlot(pbmc, reduction = 'umap')

# 参考： https://mp.weixin.qq.com/s/enGx9_Sv5wKLdtygL7b4Jw 

VlnPlot(pbmc, features = c("MS4A1", "CD79A")) 
FeaturePlot(pbmc, features = c("MS4A1", "GNLY", "CD3E", 
                               "CD14", "FCER1A", "FCGR3A", 
                               "LYZ", "PPBP", "CD8A"))

new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T",
                     "B", "CD8 T", "FCGR3A+ Mono", "NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, reduction = 'umap', 
        label = TRUE, pt.size = 0.5) + NoLegend()


# https://satijalab.org/seurat/archive/v3.0/visualization_vignette.html
VlnPlot(pbmc, features = c("MS4A1", "CD79A")) 
FeaturePlot(pbmc, features = c("MS4A1", "CD79A"))
RidgePlot(pbmc, features = c("MS4A1", "CD79A"), ncol = 1)
features= c('IL7R', 'CCR7','CD14', 'LYZ',  'IL7R', 'S100A4',"MS4A1", "CD8A",'FOXP3',
            'FCGR3A', 'MS4A7', 'GNLY', 'NKG7',
            'FCER1A', 'CST3','PPBP')
DotPlot(pbmc, features = unique(features)) + RotatedAxis()
DoHeatmap(subset(pbmc ), features = features, size = 3)
DoHeatmap(subset(pbmc, downsample = 100), 
          features = features, size = 3)

save(pbmc,file = 'basic.sce.pbmc.Rdata')





