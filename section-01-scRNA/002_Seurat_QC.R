rm(list = ls()) 
Sys.setenv(R_MAX_NUM_DLLS=999)
options(stringsAsFactors = F)
## 首先载入数据,GSE127235处理后的数据
load(file='./input.Rdata')

library(pheatmap)
library(ggplot2)
library(loomR)
library(SeuratDisk)

fivenum(apply(counts,1,function(x) sum(x>0) ))
boxplot(apply(counts,1,function(x) sum(x>0) ))
fivenum(apply(counts,2,function(x) sum(x>0) ))
hist(apply(counts,2,function(x) sum(x>0) ))
# 上面检测了 counts 和 meta 两个变量，后面需要使用
# using raw counts is the easiest way to process data through Seurat.
counts[1:4,1:4];dim(counts)
library(stringr) 
table(str_split(colnames(counts),'_',simplify = T)[,4])
#Control Diabetes 
#407      426 

meta <- metadata
head(meta) 

#二.创建Seurat对象,根据Seurat官网的标准流程
library(Seurat)
###1.以sce为例创建Seurat对象####
#sce.meta <- data.frame(Patient_ID=group$Patient_ID,row.names = group$sample)
#head(sce.meta)
#table(sce.meta$Patient_ID)
#使用scanpy建立Seurat对象后直接调至细胞周期评分
sce <- CreateSeuratObject(counts = counts,
                          meta.data =meta,
                          project = "sce")
#不使用scanpy
#文献中(2) cells with <200 genes were excluded,(min.features) = 200
sce <- CreateSeuratObject(counts = counts,
                         meta.data =meta,
                         project = "sce", 
                         min.cells = 10,#少于10个细胞覆盖的基因被过滤 
                         min.features = 200)
sce
#An object of class Seurat 
#14162 features across 826 samples within 1 assay 
#Active assay: RNA (14162 features, 0 variable features)
#过滤掉了一万多个基因，但只过滤掉7个细胞，但文献中只得到318个糖尿病细胞

####2.数据质控####

#质控的参数主要有两个： (1)每个细胞测到的unique feature数目
#（unique feature代表一个细胞检测到的基因的数目，可以根据数据的质量进行调整） 
#(2).每个细胞检测到的线粒体基因的比例，理论上线粒体基因组与核基因组相比，
#只占很小一部分。所以线粒体基因表达比例过高的细胞会被过滤。

#提取并计算线粒体基因比例
sce[["percent.mt"]] <- PercentageFeatureSet(sce, pattern = "^mt-") #人为"^MT-"
VlnPlot(sce, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

VlnPlot(object = sce, 
        features = c("nFeature_RNA", "nCount_RNA" ), 
        group.by = 'treatment',
        ncol = 2)
#图中显示percent.mt为0，可能线粒体基因被过滤过？
table(grepl("ND1",rownames(sce)))
#经检索ND1等存在，但不带mt-,使用下述方法手动添加
mt.genes <- c("ND1","ND2","ND3","ND4","ND4L","ND5","ND6",
              "CYTB","COX1","COX2","COX3","ATP6","ATP8")
             # "TA","TN","TR","TD","TC","TE","TQ","TG","TH",
              #"TI","TL1","TL2","TK","TM","TF","TP","TS1",
              #"TS2","TT","TW","TY","TV","RNR1","RNR2")
C<-GetAssayData(object = sce, slot = "counts")
percent.mt <- Matrix::colSums(C[mt.genes,])/Matrix::colSums(C)*100
sce <- AddMetaData(sce, percent.mt, col.name = "percent.mt")

##nFeature_RNA代表每个细胞测到的基因数目
#nCount_RNA代表每个细胞测到所有基因的表达量之和,即基因的所有转录本，因为一个基因可能不只有一个转录本，所以n_Count通常大于n_Feature
#percent.mt代表测到的线粒体基因的比例

#提取计算核糖体基因比例，核糖体基因的命名规则通常人是RP，小鼠是Rp
rb.genes <- rownames(sce)[grep("^Rp[sl]",rownames(sce))]#人使用"^RP[SL]"
C<-GetAssayData(object = sce, slot = "counts")
percent.ribo <- Matrix::colSums(C[rb.genes,])/Matrix::colSums(C)*100
sce <- AddMetaData(sce, percent.ribo, col.name = "percent.ribo")

## 提取ERCC基因
#方法1
table(grepl("^ERCC-",rownames(sce)))#查看是否有ERCC基因
#FALSE 
#14162 说明没有ERCC基因
#External RNA Control Consortium，是常见的已知浓度的外源RNA分子spike-in的一种
#指标含义类似线粒体含量，ERCC含量大，则说明total sum变小，往往趋势和n_Feature相反
sce[["percent.ERCC"]] <- PercentageFeatureSet(sce, pattern = "^ERCC-")
#网上看了相关教程，一般ERCC占比不高于10%
sce <- subset(sce, subset = percent.ERCC < 10)
#方法2
ercc.genes <- grep(pattern = "^ERCC-", x = rownames(x = sce), value = TRUE)
percent.ercc <- Matrix::colSums(sce[ercc.genes, ]) / Matrix::colSums(sce)
# AddMetaData adds columns to object@meta.data, and is a great place to stash QC stats
sce <- AddMetaData(object = sce, metadata = percent.ercc,
                   col.name = "percent.ercc")

#计算G2/M Score和S Score
#用到Seurat的CellCycleScoring()函数，依据的基因集已内置在Seurat中
s.genes <- Seurat::cc.genes.updated.2019$s.genes
g2m.genes <- Seurat::cc.genes.updated.2019$g2m.genes
# 使用CellCycleScoring函数计算细胞周期评分
sce <- CellCycleScoring(sce, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
# 查看细胞周期评分及分期（phase assignments）
head(sce[[]])

# seurat对象转换为loom文件,使用SeuratDisk包
adata.loom <- as.loom(x = sce, filename = "C:/Users/bdgcx/OneDrive/python/scanpy/adata.loom", verbose = FALSE)
# Always remember to close loom files when done
adata.loom$close_all()
# Visualize the distribution of cell cycle markers across
#山峦图
RidgePlot(sce, features = c("PCNA", "TOP2A", "MCM6", "MKI67"), ncol = 2)

p <- VlnPlot(sce, features =c("nFeature_RNA"),group.by = "EXP",ncol = 2)+ylim(0,8000) +xlab(" ")
png(filename = "1.png",width = 1980,height = 1500,units = "px",res = 400)
p
dev.off()
p <- VlnPlot(sce, features = c("percent.mt","percent.ribo"), ncol = 2,group.by = "EXP")
png(filename = "2.png",width = 1980,height = 1500,units = "px",res = 400)
p
dev.off()

p <- VlnPlot(sce, features = c("G2M.Score","S.Score"), ncol = 2,group.by = "EXP")
png(filename = "3.png",width = 1980,height = 1500,units = "px",res = 400)
p
dev.off()
#还可以分析免疫相关基因，缺氧相关基因等有待进一步探索

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(sce, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(sce, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 <- FeatureScatter(sce, feature1 = "nCount_RNA", feature2 = "percent.ribo")
scatter <- plot1 + plot2+plot3

png(filename = "scatter.png",width = 4500,height = 1500,units = "px",res = 400)
scatter
dev.off()
#去除线粒体基因表达比例过高的细胞，和一些极值细胞

#文献中cells with ,200 genes were excluded
#(3) cells with >30% mitochondrial RNA reads were excluded.
sce <- subset(sce, subset = nFeature_RNA > 200 & percent.mt < 30)
#sce <- subset(sce, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 30)
#图中显示percent.mt为0，可能线粒体基因被过滤过？文献中没有给极值信息所以此步跳过
dim(sce)
#14162  704  文献中得到644个细胞，可能线粒体基因仍没添加全？ 

#还可以分析免疫相关基因，缺氧相关基因等有待进一步探索

###3.标准化####
#文献中The raw read counts were normalized per cell using the NormalizeData function
#by dividing the total number of reads in that cell, then multiplying by a scale factor of 10,000 and taking log2transformed values
sce <- NormalizeData(sce, normalization.method = "LogNormalize", scale.factor = 10000)
#鉴定细胞间表达量高变的基因（feature selection）
#这一步的目的是鉴定出细胞与细胞之间表达量相差很大的基因，用于后续鉴定细胞类型，
#We selected 2167 highly variable genes on the basis ofthe average expression and dispersion per gene using FindVariableGenes function with parameters x.low.cutoff=0.1 (lower bound 0.1 for average expression), x.high.cutoff=3 (upper bound 3 for average expression), and y.cutoff=1 (low bound 1 for dispersion)
sce <- FindVariableFeatures(sce, selection.method = "vst", nfeatures = 2167,
                           x.low.cutoff=0.1,x.high.cutoff=3,y.cutoff=1)
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(sce), 10)
# plot variable features with and without labels
plot1 <- VariableFeaturePlot(sce)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = T)
p <- plot1+plot2
png(filename = "HVG_scatte.png",width = 4900,height = 2100,units = "px",res = 400)
p
dev.off()

#对矩阵进行Scale,如果保留的线粒体基因较多，因此我这里选择矫正的参数为"percent.mt","nCount_RNA"和"percent.ribo"
#文献We regressed out cell-cell variation in gene expression driven by the number of reads, mitochondrial gene content, and ribosomal gene content using the ScaleData function.
sce <- ScaleData(sce,vars.to.regress=c("percent.mt","nCount_RNA","percent.ribo"))

#有需要的话可以标准化细胞周期，使用ScaleData函数进行数据标准化，并设置vars.to.regress参数指定对细胞周期评分进行回归处理，消除细胞周期异质性的影响
sce <- ScaleData(sce,vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(sce))



