rm(list = ls())
options(stringsAsFactors = F)
#加载包
library(Seurat)
library(ggplot2)
library(monocle)
## 首先创建对象
#monocle需要的用来构建 CellDataSet 对象的三个数据集
#1.表达量矩阵exprs:数值矩阵 行名是基因, 列名是细胞编号.
#2.细胞的表型信息phenoData(pd): 第一列是细胞编号，其他列是细胞的相关信息
#3.基因注释featureData(fd): 第一列是基因编号, 其他列是基因对应的信息
#并且这三个数据集要满足如下要求:
#表达量矩阵必须：
#保证它的列数等于phenoData的行数
#保证它的行数等于featureData的行数
#而且phenoData的行名需要和表达矩阵的列名匹配
#featureData和表达矩阵的行名要匹配
#featureData至少要有一列“gene_short_name”, 就是基因的symbol

# 将Seurat中的对象转换为monocle识别的对象
#cds <- importCDS(GetAssayData(sce))
#cds <- importCDS(sce)  不知道为啥，使用imporCDS总报错,所以手动载入对象

## 手动载入Seurat的数据
sce <- readRDS("sce_tutorial.rds")
#从Seurat对象中选择做拟时序的细胞群cluster
Mono_tj<-subset(sce, idents = c("EC"))
rm(sce)
#Extract data, phenotype data, and feature data from the SeuratObject
data <- as(as.matrix(Mono_tj@assays$RNA@counts), 'sparseMatrix')#转变为稀疏矩阵
pd <- new('AnnotatedDataFrame', data = Mono_tj@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)

#Construct monocle cds
monocle_cds <- newCellDataSet(data,
                              phenoData = pd,
                              featureData = fd,
                              expressionFamily = negbinomial.size())

HSMM<-monocle_cds
## 归一化 
HSMM <- estimateSizeFactors(HSMM)
HSMM <- estimateDispersions(HSMM)
dim(HSMM)
#Filtering low-quality cells
##基因层面的过滤
## 起初是 14162 features, 379 samples 
HSMM <- detectGenes(HSMM, min_expr = 3 )
print(head(fData(HSMM)))

expressed_genes <- row.names(subset(fData(HSMM),
                                    num_cells_expressed >= 5)) 

length(expressed_genes)
HSMM <- HSMM[expressed_genes,]
HSMM
# 过滤基因后是  10986 features, 379 samples 
print(head(pData(HSMM)))

#作图，查看标准化的基因表达值的分布，也可跳过此步
L <- log(exprs(HSMM[expressed_genes,]))
#将每个基因都标准化，melt方便作图
library(reshape)
melted_dens_df <- melt(Matrix::t(scale(Matrix::t(L))))
qplot(value, geom = "density", data = melted_dens_df) +
  stat_function(fun = dnorm, size = 0.5, color = 'red') +
  xlab("Standardized log(FPKM)") +
  ylab("Density")

#聚类
#挑选合适的基因进入下游分析，比如PCA或者tSNE
#第一步是决定用哪些基因来进行聚类（特征选择）我们可以根据平均表达水平筛选基因，我们还可以选择细胞间异常变异的基因。这些基因往往对细胞状态有很高的信息含量。
disp_table <- dispersionTable(HSMM)
#筛选平均表达量 > 0.1的基因用于聚类
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)
HSMM <- setOrderingFilter(HSMM, unsup_clustering_genes$gene_id)
plot_ordering_genes(HSMM)

plot_pc_variance_explained(HSMM, 
                           return_all = F,
                           max_components = 20)

# 其中 num_dim 参数选择基于上面的PCA图
HSMM <- reduceDimension(HSMM, 
                       max_components = 2, 
                       num_dim = 10,
                       reduction_method = 'tSNE', 
                       verbose = T)

HSMM <- clusterCells(HSMM, num_clusters = 2)

plot_cell_clusters(HSMM,color_by = 'Cluster')
plot_cell_clusters(HSMM,color_by = 'Cluster',
                   markers = c("Kdr"))
#Monocle允许我们减去“无趣的”变异源的影响，以减少它们对集群的影响。也可以不用去除
#可以使用到clusterCells和其他几个Monocle函数的residualmodelformula astr参数来实现这一点。此参数接受一个R模型公式字符串，该字符串指定要在cluster之前减去。
HSMM <- reduceDimension(HSMM, max_components = 2, num_dim = 2,
                        reduction_method = 'tSNE',
                        residualModelFormulaStr = "~Size_Factor + num_genes_expressed",
                        verbose = T)
HSMM <- clusterCells(HSMM, num_clusters = 5)
plot_cell_clusters(HSMM, 1, 2, color = "Cluster")  

plot_cell_clusters(HSMM, 1, 2, color = "Cluster") +
  facet_wrap(~CellType)

#构建轨迹分3步
#Step 1: choosing genes that define progress
#Step 2: reducing the dimensionality of the data
#Step 3: ordering the cells in pseudotime

#Step 1基因的选择
#差异分析
Sys.time()

diff_test_res <- differentialGeneTest(HSMM[expressed_genes],
                                      fullModelFormulaStr = "~treatment")
Sys.time()
#推断发育轨迹，挑选差异显著的基因（q<0.01）做降维，给细胞排序
ordering_genes <- row.names (subset(diff_test_res, qval < 0.01))
length(ordering_genes)

HSMM <- setOrderingFilter(HSMM, ordering_genes)
plot_ordering_genes(HSMM)
#降维
HSMM <- reduceDimension(HSMM, max_components = 2,
                        reduction_method  = 'DDRTree')
# 细胞排序
HSMM <- orderCells(HSMM)
# 可视化
plot_cell_trajectory(HSMM, color_by = "treatment")
# 展现marker基因
plot_genes_in_pseudotime(HSMM[cg,],
                         color_by = "treatment")



