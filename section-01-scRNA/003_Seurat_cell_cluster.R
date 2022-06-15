####4.细胞分类####
#4.1分类前首先要对数据集进行降维
#Scaling the data
sce <- ScaleData(sce, features = rownames(sce))
#Perform linear dimensional reduction
sce <- RunPCA(sce, features = VariableFeatures(object = sce))
#Examine and visualize PCA results a few different ways
print(sce[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(sce, dims = 1:2, reduction = "pca")
DimPlot(sce, reduction = "pca",group.by = "EXP")
DimHeatmap(sce, dims = 1, cells = 704, balanced = TRUE)
DimHeatmap(sce, dims = 1:2, cells = 704, balanced = TRUE)
#判断选择的主成分数目，根据网络教程选择19个
ElbowPlot(sce)

#4.2定义数据集的“维度”
#这里我们需要选择出主成分的数目，用于后续细胞分类。这里定义的“维度”并不代表细胞类型的数目，而是对细胞分类时需要用到的一个参数
#NOTE: This process can take a long time for big datasets, comment out for expediency. 
#More approximate techniques such as those implemented in ElbowPlot() can be used to reduce computation time
sce <- JackStraw(sce, num.replicate = 100)
sce <- ScoreJackStraw(sce, dims = 1:20)
JackStrawPlot(sce, dims = 1:15)
ElbowPlot(sce)
#JackStraw和Elbow都可以决定数据的“维度”。但是Elbow比较直观，
#我们选择Elbow结果进行解读。可以看到，主成分（PC）7到10之间，数据的标准差基本不在下降。所以我们需要在7到10之间进行选择用于细胞的分类
#但文献中The top 15 principal componentswere chosen for cell clustering and t-SNE projection because no significant changes were observed beyond 15 principal components，所以这里也用15
#4.3 细胞分类
#选择不同的resolution值可以获得不同的cluster数目，值越大cluster数目越多，默认值是0.5
#文献标准cells were clustered using FindClusters function with resolution=0.7. 
#Each cluster was screened for marker genes by differential expression analysis (DEA) between cells inside and outside the cluster using Finscearkers function with parameters min.pct=0.25 (genes expressed in at least 25% of cells either inside or outside ofa cluster) and test.use=“wilcox”
sce <- FindNeighbors(sce, dims = 1:15)
sce <- FindClusters(sce,resolution = 0.4)
#看看每个cluster的细胞数目
table(sce@meta.data$RNA_snn_res.0.4)
#这里我们设置了dims = 1:15 即选取前15个主成分来分类细胞。分类的结果如下,可以看到，细胞被分为5个类别。
#Look at cluster IDs of the first 5 cells
head(Idents(sce), 5)

#4.4可视化分类结果
#TSNE和UMAP两种方法经常被用于可视化细胞类别。

#UMAP
set.seed(1234)#设置随机数种子，这样就不会每次图都变化
sce <- RunUMAP(sce, dims = 1:15, label = T)
head(sce@reductions$umap@cell.embeddings) # 提取UMAP坐标值。
#note that you can set `label = TRUE` or use the LabelClusters function to help label individual clusters
DimPlot(sce[,sce@meta.data$treatment=="Control"], 
              reduction = "umap"
              )
DimPlot(sce[,sce@meta.data$treatment=="Control"], reduction = "umap")+
  theme(panel.background = element_rect(fill='white', colour='black'), 
                                  panel.grid=element_blank()) 

#T-SNE
set.seed(123)
sce <- RunTSNE(sce, dims = 1:15)
head(sce@reductions$tsne@cell.embeddings)
DimPlot(sce[,sce@meta.data$treatment=="Diabetes"], 
              reduction = "tsne"
              )
DimPlot(sce,reduction = "tsne",group.by = "orig.ident",split.by = "treatment")+
  theme(panel.background = element_rect(fill='white', colour='black'), 
        panel.grid=element_blank()) 
ggsave("./figure/exp_trt_tsne.png",width = 3960,height = 1700,dpi = 300,units = "px")

DimPlot(sce,reduction = "tsne",group.by = "EXP")+
    theme(panel.background = element_rect(fill='white', colour='black'), 
      panel.grid=element_blank())
ggsave("./figure/exp_tsne.png",width = 1980,height = 1700,dpi = 300,units = "px")

saveRDS(sce, file = "sce_tutorial.rds")  #保存数据，用于后续个性化分析

#给图片添加meta信息
tsne_pos=Embeddings(sce,'tsne') 
DimPlot(sce,reduction = "tsne",group.by  ='orig.ident')

DimPlot(sce,reduction = "tsne",label=T,split.by ='orig.ident')

#ggplot可视化
dat=cbind(tsne_pos,sce@meta.data)
library(ggplot2)
cluster <- sce@meta.data$seurat_clusters
p=ggplot(dat,aes(x=tSNE_1,y=tSNE_2,color=cluster))+geom_point(size=0.95)
p=p+stat_ellipse(data=dat,aes(x=tSNE_1,y=tSNE_2,fill=cluster,color=cluster),
                 geom = "polygon",alpha=0.2,level=0.9,type="t",linetype = 2,show.legend = F)+coord_fixed()
theme= theme(panel.grid =element_blank()) +  ##删去网格 
  theme(panel.border = element_blank(),panel.background = element_blank()) +   ##删去外层边框
  theme(axis.line = element_line(size=1, colour = "black")) 
p=p+theme+guides(colour = guide_legend(override.aes = list(size=5)))
print(p)
ggsave("tsnegg.png",p,width = 1980,height = 1800,units = "px",dpi = 300)
