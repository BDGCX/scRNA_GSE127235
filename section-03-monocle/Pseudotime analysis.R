rm(list = ls())
options(stringsAsFactors = F)
#加载包
library(Seurat)
library(ggplot2)
library(monocle)
library(dplyr)
library(patchwork)
## 首先创建对象
#monocle需要的用来构建 CellDataSet 对象的三个数据集
#1.表达量矩阵exprs:数值矩阵 行名是基因, 列名是细胞编号.
#2.细胞的表型信息phenoData(pd): 第一列是细胞编号，其他列是细胞的相关信息，即meta.data
#3.基因注释featureData(fd): 第一列是基因编号, 其他列是基因对应的信息,即基因名称
#并且这三个数据集要满足如下要求:
#表达量矩阵必须：
#保证它的列数等于phenoData的行数
#保证它的行数等于featureData的行数
#而且phenoData的行名需要和表达矩阵的列名匹配
#featureData和表达矩阵的行名要匹配
#featureData至少要有一列“gene_short_name”, 就是基因的symbol

####1.加载seurat数据集####
getOption('timeout')
options(timeout=10000)

sce <- readRDS("sce_tutorial.rds")
DimPlot(pbmc, reduction = 'umap', 
        label = TRUE, pt.size = 0.5) + NoLegend()
               
table(Idents(sce))
#从Seurat对象中选择做拟时序的细胞群cluster
sce<-subset(sce, idents = c("EC"))
#sce <-  sce[,Idents(sce)  %in% c("EC")]这步等同于subset
table(Idents(sce))
levels(Idents(sce))
####check monocytes(此步可省略)####
markers_df  <- FindMarkers(object = sce, 
                           ident.1 = "Diabetes",#病例组。ident.2=NULL时默认剩下所有细胞为对照组
                           ident.2 = "Control",
                           group.by = "treatment",
                           subset.ident="EC",
                           #logfc.threshold = 0,
                           min.pct = 0.25)
head(markers_df)
# 在seurat V4里面， 是 avg_log2FC ， 但是如果是V3版本包，仍然是 avg_logFC
cg_markers_df <- markers_df[abs(markers_df$avg_log2FC) >1,]
#cg_markers_df <- markers_df[abs(markers_df$avg_logFC) >1,]
dim(cg_markers_df) 
cg_markers_df <- cg_markers_df[order(cg_markers_df$avg_log2FC),]
DotPlot(sce,
        features = rownames(cg_markers_df)) + 
  theme(axis.text.x = element_text(angle = 45, 
                                   vjust = 0.5, 
                                   hjust=0.5))
DoHeatmap(sce,
          features = rownames(cg_markers_df)) 
####2.创建monocle对象####
#Extract data, phenotype data, and feature data from the SeuratObject
data <- as(as.matrix(sce@assays$RNA@counts), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = sce@meta.data)
fData <- data.frame(gene_short_name = row.names(data), 
                    row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)

#Construct monocle cds
sc_cds <- newCellDataSet(data,
                         phenoData = pd,
                         featureData = fd,
                         lowerDetectionLimit = 0.5,
                         expressionFamily = negbinomial.size())

sc_cds
#Estimate size factors and dispersions归一化
cds <- sc_cds
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds) 
#Filtering low-quality cells
cds <- detectGenes(cds, min_expr = 1)
print(head(fData(cds)))

expressed_genes <- row.names(subset(fData(cds),
                                    num_cells_expressed >= 10))
print(head(pData(cds)))
#以下步骤可省略
L <- log(exprs(cds[expressed_genes,]))
#将每个基因都标准化，melt方便作图
library(reshape)
melted_dens_df <- melt(Matrix::t(scale(Matrix::t(L))))
#作图，查看标准化的基因表达值的分布
qplot(value, geom = "density", data = melted_dens_df) +
  stat_function(fun = dnorm, size = 0.5, color = 'red') +
  xlab("Standardized log(FPKM)") +
  ylab("Density")

cds <- cds[fData(sc_cds)$num_cells_expressed > 10, ]
cds
####3.Classifying and Counting Cells####
#并不是所有的基因都有作用，所以先进行挑选，合适的基因用来进行聚类。
#官网提供四种细胞分类方法，此处以第2、3种为例
#1.Classifying cells by type
#2.Clustering cells without marker genes
#3.Clustering cells using marker genes
#4.Imputing cell type

##2.Clustering cells without marker genes
disp_table <- dispersionTable(cds)
unsup_clustering_genes <- subset(disp_table,
                                 mean_expression >= 0.1)
cds <- setOrderingFilter(cds, 
                         unsup_clustering_genes$gene_id)
plot_ordering_genes(cds)#红线表示单片基于这种关系对色散的期望。我们标记用于聚类的基因用黑点表示，其他的用灰点表示。 
plot_pc_variance_explained(cds, return_all = F) # norm_method='log'

#3.Clustering cells using marker genes (Recommended)
MYF5_id <- row.names(subset(fData(cds), gene_short_name == "MYF5"))
ANPEP_id <- row.names(subset(fData(cds),
                             gene_short_name == "ANPEP"))

cth <- newCellTypeHierarchy()
cth <- addCellType(cth, "Myoblast", classify_func =
                     function(x) { x[MYF5_id,] >= 1 })
cth <- addCellType(cth, "Fibroblast", classify_func = function(x)
{ x[MYF5_id,] < 1 & x[ANPEP_id,] > 1 })

marker_diff <- markerDiffTable(cds[expressed_genes,],
                               cth,
                               residualModelFormulaStr = "~treatment + num_genes_expressed",
                               cores = 1)
candidate_clustering_genes <-row.names(subset(marker_diff, 
                                              qval < 0.01))
marker_spec <-calculateMarkerSpecificity(cds[candidate_clustering_genes,], 
                                         cth)
head(selectTopMarkers(marker_spec, 3))
#To cluster the cells, we'll choose the top 500 markers for each of these cell types
semisup_clustering_genes <- unique(selectTopMarkers(marker_spec, 500)$gene_id)
HSMM <- setOrderingFilter(cds, semisup_clustering_genes)
plot_ordering_genes(cds)
plot_pc_variance_explained(HSMM, return_all = F)

#接2、3种方法
# 其中 num_dim 参数选择基于上面的PCA图
cds <- reduceDimension(cds, max_components = 2, 
                       num_dim = 6,
                       reduction_method = 'tSNE', 
                       verbose = T)
cds <- clusterCells(cds, num_clusters = 2) 
plot_cell_clusters(cds, 1, 2,color_by ="Cluster",
                   markers =c("Lyz2","Ces2e") )
plot_cell_clusters(cds, 1, 2,color_by ="treatment")
#减去“无趣的”变异源的影响，以减少它们对集群的影响
cds <- reduceDimension(cds, max_components = 2, num_dim = 2,
                       reduction_method = 'tSNE',
                       residualModelFormulaStr = "~Size_Factor + num_genes_expressed",
                       verbose = T)
cds <- clusterCells(cds, num_clusters = 2)
plot_cell_clusters(cds, 1, 2, color = "treatment")  # 图不放了

plot_cell_clusters(cds, 1, 2, color = "Cluster") +
  facet_wrap(~treatment)

####4.Constructing Single Cell Trajectories构建轨迹#####
#构建轨迹分3步
#Step 1: choosing genes that define progress
#Step 2: reducing the dimensionality of the data
#Step 3: ordering the cells in pseudotime

# 接下来很重要，到底是看哪个性状的轨迹
colnames(pData(cds))
table(pData(cds)$Cluster)
table(pData(cds)$Cluster,pData(cds)$new)

plot_cell_clusters(cds, 1, 2 )

## 我们这里并不能使用 monocle的分群
# 还是依据前面的 seurat分群, 其实取决于自己真实的生物学意图
pData(cds)$Cluster <- pData(cds)$celltype#cell type为Seurat分群
table(pData(cds)$Cluster)

#理想情况下，我们希望尽可能少地使用正在研究的系统生物学的先验知识。我们希望从数据中发现重要的排序基因，而不是依赖于文献和教科书，因为这可能会在排序中引入偏见。
#我们将从一种更简单的方法开始，但是我们通常推荐一种更复杂的方法，称为“dpFeature”。
#方法1：简单方法
#Step 1:基因的选择
Sys.time()
diff_test_res <- differentialGeneTest(cds[expressed_genes,],
                                      fullModelFormulaStr = "~treatment")
Sys.time()
ordering_genes <- row.names (subset(diff_test_res, qval < 0.01))
cds <- setOrderingFilter(cds, ordering_genes)
plot_ordering_genes(cds)

#step 2: reduce data dimensionality 降维 
cds <- reduceDimension(cds, max_components = 2,
                       method = 'DDRTree')
#step 3: order cells along the trajectory排序
cds <- orderCells(cds)
plot_cell_trajectory(cds, color_by = "treatment")

plot_cell_trajectory(cds, color_by = "State")

plot_cell_trajectory(cds, color_by = "Pseudotime")
#"State" is just Monocle's term for the segment of the tree. The function below is handy for identifying the State which contains most of the cells from time zero
GM_state <- function(cds){
  if (length(unique(pData(cds)$State)) > 1){
    T0_counts <- table(pData(cds)$State, pData(cds)$treatment)[,"0"]
    return(as.numeric(names(T0_counts)[which
                                       (T0_counts == max(T0_counts))]))
  } else {
    return (1)
  }
}
cds_1 <- orderCells(cds, root_state = GM_state(cds))

plot_cell_trajectory(cds_1, color_by = "Pseudotime")

plot_cell_trajectory(cds, color_by = "State") +
  facet_wrap(~State, nrow = 1)
#And if you don't have a timeseries, you might need to set the root based on where certain marker genes are expressed, using your biological knowledge of the system. 
blast_genes <- row.names(subset(fData(cds),
                                gene_short_name %in% c("CCNB2", "MYOD1", "MYOG")))
plot_genes_jitter(cds[blast_genes,],
                  grouping = "State",
                  min_expr = 0.1)
#单个基因的时间变化（要选有意义的基因）
cds_expressed_genes <-  row.names(subset(fData(cds),
                                          num_cells_expressed >= 10))
cds_filtered <- cds[cds_expressed_genes,]
my_genes <- row.names(subset(fData(cds_filtered),
                             gene_short_name %in% c("CDK1", "MEF2C", "MYH3")))
cds_subset <- HSMM_filtered[my_genes,]
plot_genes_in_pseudotime(cds_subset, color_by = "treatment")

#2.复杂方法：dpFeature根据伪时间表达pattern聚类基因(推荐)
#To use dpFeature, we first select superset of feature genes as genes expressed in at least 5% of all the cells
#Step 1:基因的选择
cds <- detectGenes(cds, min_expr = 0.1)
fData(cds)$use_for_ordering <-fData(cds)$num_cells_expressed > 0.05 * ncol(cds)
#然后就是进行PCA分析，根据下面的图选择合适的PC
plot_pc_variance_explained(cds, return_all = F)

##step 2: reduce data dimensionality 降维 ，这里的num_dim我用的是seurat中选择的PC数
cds <- reduceDimension(cds,
                       max_components = 2,
                       norm_method = 'log',
                       num_dim = 6,
                       reduction_method = 'tSNE',
                       verbose = T)

#随后进行密度峰值聚类，鉴定cluster
#The densityPeak algorithm clusters cells based on each cell's local density (Ρ) and the nearest distance (Δ) of a cell to another cell with higher distance
#可以人为自己设定Ρ和Δ
#默认参数进行聚类
cds <- clusterCells(cds,verbose = F)
#查看聚类结果，主要是看一下用聚类分出来的结果和实验分群结果（seurat分群结果）是否一致，如果一致，那就可以进一步找这些cluster中的基因，作为分群的基因，不一致就要自己调参数
plot_cell_clusters(cds, color_by = 'as.factor(Cluster)')
plot_cell_clusters(cds, color_by = 'as.factor(treatment)')

#查看每一个细胞的Ρ和Δ，然后自己设定阈值，图片中的黑色点代表cluster的数量
plot_rho_delta(cds, rho_threshold = 2, delta_threshold = 12)

#按照设定的阈值重新聚类
cds <- clusterCells(cds,
                    rho_threshold = 2,
                    delta_threshold = 12,
                    skip_rho_sigma = T,
                    verbose = F)
#查看聚类结果，两种结果图很一致，那么就找这些群的基因
plot_cell_clusters(cds, color_by = 'as.factor(Cluster)')
plot_cell_clusters(cds, color_by = 'as.factor(treatment)')

#这里的res.0.2是fData(cds)的列名，可以根据自己的需要进行选择，我这里选择的是分辨率为0.2的seurat聚类结果。
#head(fData(cds))
cds_expressed_genes <-  row.names(subset(fData(cds),
                                         num_cells_expressed >= 10))
#这一步耗时很长
clustering_DEG_genes <-differentialGeneTest(cds[cds_expressed_genes,],
                                            fullModelFormulaStr = "~treatment",
                                            cores = 1)#fullModelFormulaStr的值是meta里的

#We will then select the top 1000 significant genes as the ordering genes.
##选取前1000个基因进行拟时间轨迹分析
##和之前的步骤一样，第一步确定合适的基因
cds_ordering_genes <-row.names(clustering_DEG_genes)[order(clustering_DEG_genes$qval)][1:1000]
#将基因添加到cds对象中
cds <-setOrderingFilter(cds,
                        ordering_genes = cds_ordering_genes)
plot_ordering_genes(cds)
#第二步用DDRTree进行降维
cds <-reduceDimension(cds, method = 'DDRTree')
#第三步排序构建拟时间曲线
cds <-orderCells(cds)
#cds <-orderCells(cds, root_state = GM_state(cds))
p1 <- plot_cell_trajectory(cds, color_by = "treatment")+
  theme(panel.background = element_rect(fill='white', colour='black'), 
              panel.grid=element_blank())

p3 <- plot_cell_trajectory(cds, color_by = "State")
p1+p2/p3
ggsave("./figure/t2.png",width = 3500,height = 1900,dpi = 300,units = "px")
#"State" is just Monocle's term for the segment of the tree. The function below is handy for identifying the State which contains most of the cells from time zero
#根据经验人为设定拟时序的起点
GM_state <- function(cds){
  if (length(unique(pData(cds)$State)) > 1){
    T0_counts <- table(pData(cds)$State, pData(cds)$treatment)[,"Control"]#拟时序的根点即0点
    return(as.numeric(names(T0_counts)[which
                                       (T0_counts == max(T0_counts))]))
  } else {
    return (1)
  }
}
cds <- orderCells(cds, root_state = GM_state(cds))

p2 <- plot_cell_trajectory(cds, color_by = "Pseudotime")
p1+p2/p3
ggsave("./figure/t2.png",width = 3500,height = 1900,dpi = 300,units = "px")

plot_cell_trajectory(cds, color_by = "State") +
  facet_wrap(~State, nrow = 1)
#And if you don't have a timeseries, you might need to set the root based on where certain marker genes are expressed, using your biological knowledge of the system. 
blast_genes <- row.names(subset(fData(cds),
                                gene_short_name %in% c("Lyz2", "Kdr")))
plot_genes_jitter(cds[blast_genes,],
                  grouping = "treatment",
                  color_by = "State",
                  min_expr = 0.1,
                  nrow= 1,
                  ncol = NULL,
                  plot_trend = TRUE)
#单个基因的时间变化（要选有意义的基因）
cds_expressed_genes <-  row.names(subset(fData(cds),
                                         num_cells_expressed >= 10))
cds_filtered <- cds[cds_expressed_genes,]
my_genes <- row.names(subset(fData(cds_filtered),
                             gene_short_name %in% c("Lyz2", "Kdr")))
cds_subset <- cds_filtered[my_genes,]
plot_genes_in_pseudotime(cds_subset, color_by = "treatment")

####5.Differential Expression Analysis差异分析####
#官方给出的差异分析有三大方法，我们重点关注第三个
#1、Basic Differential Analysis
#2、Finding Genes that Distinguish Cell Type or State
#3、Finding Genes that Change as a Function of Pseudotime

#1.Basic Differential Analysis
marker_genes <- row.names(subset(fData(cds),
                                 gene_short_name %in% c("MEF2C", "MEF2D", "MYF5",
                                                        "ANPEP", "PDGFRA","MYOG",
                                                        "TPM1",  "TPM2",  "MYH2",
                                                        "MYH3",  "NCAM1", "TNNT1",
                                                        "TNNT2", "TNNC1", "CDK1",
                                                        "CDK2",  "CCNB1", "CCNB2",
                                                        "CCND1", "CCNA1", "ID1")))
diff_test_res <- differentialGeneTest(cds[marker_genes,],
                                      fullModelFormulaStr = "~treatment")

#Select genes that are significant at an FDR < 10%
sig_genes <- subset(diff_test_res, qval < 0.1)

sig_genes[,c("gene_short_name", "pval", "qval")]

MYOG_ID1 <- HSMM_myo[row.names(subset(fData(cds),
                                      gene_short_name %in% c("MYOG", "CCNB2"))),]
plot_genes_jitter(MYOG_ID1, grouping = "treatment", ncol= 2)

#2.Finding Genes that Distinguish Cell Type or State选择对区分细胞类型有意义的基因
#Finding Genes that Distinguish Cell Type or State 
to_be_tested <- row.names(subset(fData(cds),
                                 gene_short_name %in% c("CDC20", "NCAM1", "ANPEP")))
cds_subset <- cds[to_be_tested,]

diff_test_res <- differentialGeneTest(cds_subset,
                                      fullModelFormulaStr = "~treatment")
diff_test_res[,c("gene_short_name", "pval", "qval")]

plot_genes_jitter(cds_subset,
                  grouping = "treatment",
                  color_by = "treatment",
                  nrow= 1,
                  ncol = NULL,
                  plot_trend = TRUE)

full_model_fits <-
  fitModel(cds_subset,  modelFormulaStr = "~treatment")
reduced_model_fits <- fitModel(cds_subset, modelFormulaStr = "~1")
diff_test_res <- compareModels(full_model_fits, reduced_model_fits)
diff_test_res

#3.Finding Genes that Change as a Function of Pseudotime
#Monocle的主要工作是将细胞按照生物过程(如细胞分化)的顺序排列，而不事先知道要查看哪些基因。一旦这样做了，你就可以分析细胞，找到随着细胞进步而变化的基因。
to_be_tested <- row.names(subset(fData(cds),
                                 gene_short_name %in% c("Lyz2", "MEF2C", "CCNB2", "TNNT1")))
cds_subset <- cds[to_be_tested,]

diff_test_res <- differentialGeneTest(cds_subset,
                                      fullModelFormulaStr = "~sm.ns(Pseudotime)")

diff_test_res[,c("gene_short_name", "pval", "qval")]

plot_genes_in_pseudotime(cds_subset, color_by = "seurat_clusters")

#Clustering Genes by Pseudotemporal Expression Pattern
#注意下面的num_clusters 指的是基因可以聚成几个类，而不是细胞
marker_genes <- row.names(subset(fData(cds),
                                 gene_short_name %in% c("Lyz2", "H2-Aa", "Ednrb",
                                                        "ANPEP", "PDGFRA","MYOG",
                                                        "TPM1",  "TPM2",  "MYH2",
                                                        "MYH3",  "NCAM1", "TNNT1",
                                                        "TNNT2", "TNNC1", "CDK1",
                                                        "CDK2",  "CCNB1", "CCNB2",
                                                        "CCND1", "CCNA1", "ID1")))
marker_genes <- as.character(rownames(markers_df))[1:21]
diff_test_res <- differentialGeneTest(cds[marker_genes,],
                                      fullModelFormulaStr = "~sm.ns(Pseudotime)")
sig_gene_names <- row.names(subset(diff_test_res, qval < 0.01))
plot_pseudotime_heatmap(cds[sig_gene_names,],
                        num_clusters = 6,
                        cores = 1,
                        show_rownames = T)

plot_pseudotime_heatmap(cds[cds_ordering_genes,],#前1000个基因
                        num_clusters = 6,
                        cores = 1,
                        trend_formula = "~sm.ns(Pseudotime, df=3)",
                        show_rownames = F)
####6.Multi-Factorial Differential Expression Analysis####
#Monocle可以在多个因素存在的情况下进行差异分析，这可以帮助你减去一些因素来看到其他因素的影响。
#在下面的简单例子中，Monocle测试了三个基因在成肌细胞和成纤维细胞之间的差异表达，同时减去percent.mt的影响。
#为此，我们必须同时指定完整模型和简化模型。完整的模型同时捕捉细胞类型和percent.mt的影响。
to_be_tested <-row.names(subset(fData(cds),
                   gene_short_name %in% c("TPM1", "MYH3", "CCNB2", "GAPDH")))

cds_subset <- cds[to_be_tested,]

diff_test_res <- differentialGeneTest(cds_subset,
                                      fullModelFormulaStr = "~treatment + percent.mt",
                                      reducedModelFormulaStr = "~percent.mt")
diff_test_res[,c("gene_short_name", "pval", "qval")]

plot_genes_jitter(cds_subset,
                  grouping = "seurat_clusters", 
                  color_by = "CellType", 
                  plot_trend = TRUE) +
  facet_wrap( ~ feature_label, scales= "free_y")

####7.Analyzing Branches in Single-Cell Trajectories单细胞轨迹的“分支”分析####
#分支表达式分析建模，或BEAM
plot_cell_trajectory(cds, color_by = "State")
#BEAM进行统计分析
BEAM_res <- BEAM(cds, branch_point = 1, cores = 8)
BEAM_res <- BEAM_res[order(BEAM_res$qval),]
BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]
plot_genes_branched_heatmap(cds[row.names(subset(BEAM_res,
                                                  qval < 1e-4)),],
                            branch_point = 1,
                            num_clusters = 4,
                            cores = 1,
                            use_gene_short_name = T,
                            show_rownames = T)

#该热图显示的是同一时间点两个谱系的变化。热图的列是伪时间的点，行是基因。
#从热图中间往右读，是伪时间的一个谱系；往左是另一个谱系。基因是被按照等级聚类的，所以你看到的基因表达模式和谱系的表达模式是非常相似的
genes <- row.names(subset(fData(cds),
                          gene_short_name %in% c( "MEF2C", "CCNB2", "TNNT1")))

plot_genes_branched_pseudotime(cds[genes,],
                               branch_point = 1,
                               color_by = "State",
                               ncol = 1)









