#BiocManager::valid()#查看过时的包
#6.利用先验知识定义细胞类型(官网教程)
#通过对比我们鉴定的marker gene与已发表的细胞类型特意的基因表达marker，可以定义我们划分出来的细胞类群。
#最后，给我们定义好的细胞类群加上名称
#参考细胞marker数据库：http://xteam.xbio.top/CellMarker/

#寻找需要注释疾病细胞的相关文献作为参考，找到文中的marker基因(参考文献中注释的细胞类型)
EC <- c("Flt1","Tie1","Pecam1", "Kdr","Emcn","Cdh5")
MC <- c("Pdgfrb", "Gata3", "Des", "Itga8")
POD <- c("Nphs1", "Nphs2", "Pdpn", "Wt1", 
         "Mafb", "Synpo", "Cdkn1c", "Ptpro")
IC <- c("Ptprc", "Lyz1", "Csf1r", "Itgam", "Cd3d", 
        "Ms4a1","Lyz2","Cd74","H2-Aa")
PEC <- c("Cldn1", "Pax8")
SMC <- c("Acta2", "Myh11", "Tagln")
TEC <- c("Fxyd2","Slc12a1", "Slc12a3", "Slc14a2","Aqp1",
         "Aqp2","Umod","Atp1b1","Cdh16")
#合并
markers <- c(EC,MC,POD,IC,SMC,TEC,PEC)

#使用CaseMatch去除不存在的基因
markers <- CaseMatch(markers,rownames(sce))
markers <- as.character(markers)

#针对marker gene做点图
DotPlot(sce, features = unique(markers),
        group.by = "seurat_clusters")+RotatedAxis()+
  scale_x_discrete("")+scale_y_discrete("")
#根据点图对细胞群进行注释
new.cluster.ids <- c("EC", "MC", "POD", "IC", "TC")
names(new.cluster.ids) <- levels(sce)
sce <- RenameIdents(sce, new.cluster.ids)
sce$new_clusters <- sce@active.ident
DimPlot(sce, reduction = "tsne",split.by = "treatment", 
        label = TRUE, pt.size = 1) + NoLegend()+
  theme(panel.background = element_rect(fill='white', colour='black'), 
        panel.grid=element_blank()) 
ggsave("./figure/tsne_trt.png",p,width = 3960,height = 1700,dpi = 300,units = "px")

sce.markers <- FindAllMarkers(object = sce, only.pos = TRUE, min.pct = 0.25, 
                               thresh.use = 0.25)#only.pos = TRUE：只返回positive的marker基因
print(x = head(sce.markers))

DT::datatable(sce.markers)
write.csv(sce.markers,file=paste0('_sce.markers_tsne.csv'))

sce@meta.data$new_clusters <- sce@active.ident
saveRDS(sce, file = "sce_tutorial.rds")

#如果想合并某几个细胞群
#方法1
levels(Idents(sce)) #查看细胞亚群，与上述结果一致
## 根据levels(Idents(sce)) 顺序重新赋予对应的 B、DC、Mono、Platelet、T 这5个细胞亚群名称，顺序不能乱
new.cluster.ids <- c("EC", "MC","EC ", "IC", "POD", "TC")
names(new.cluster.ids) <- levels(sce)
sce <- RenameIdents(sce, new.cluster.ids)
levels(sce) #查看是否已改名

DimPlot(sce, reduction = 'umap', 
        label = TRUE, pt.size = 0.5) + NoLegend()

#方法2,推荐方法2
cluster2celltype <- c("0"="EC",
                      "1"="MC", 
                      "2"="EC", 
                      "3"= "IC", 
                      "4"= "POD", 
                      "5"= "TC")
                    
sce[['cell_type']] = unname(cluster2celltype[sce@meta.data$seurat_clusters])
DimPlot(sce[,sce@meta.data$treatment=="Diabetes"], reduction = 'umap', group.by = 'cell_type',
        label = TRUE, pt.size = 0.5) + NoLegend()
table(sce[['cell_type']])

table(sce[,sce@meta.data$treatment=="Diabetes"][['cell_type']])
cell_type
#EC  IC  MC POD  TC 
#200  77  43  18  20 

#cell annotation-----自动注释
# 对肿瘤细胞来说，分群后的细胞亚群注释是不可行的

#一、Celltypist
#参考网址：https://mp.weixin.qq.com/s/CdhAp7lO5nRJW4cnLWGwLQ
library(reticulate)
#设定python环境
use_condaenv("C:/Users/bdgcx/anaconda3")
py_config()#安装的python版本环境查看，显示anaconda和numpy的详细信息。
py_available()#[1] TRUE   #检查您的系统是否安装过Python
py_module_available("pandas")#检查“pandas”是否安装
## 在R中载入python模块
scanpy <-  import("scanpy")
celltypist <-  import("celltypist")
pandas <- import("pandas")
numpy <-  import("numpy")

#### 下载参考数据集
celltypist$models$download_models(force_update = T) 


# 二、SingleR 做 cell annotation的流程(准确度较低)
library(SingleR)
??singleR
refdata <- get(load("C:/Users/bdgcx/OneDrive/single cell/analysis of scRNA-seq/rawdata/HumanPrimaryCellAtlasData.Rdata"))
assay(refdata)[1:4,1:4]
head(refdata@colData)
head(refdata)
ref <- celldex::HumanPrimaryCellAtlasData()
#参考数据库，等待时间较长。建议下载成功后，储存为Rdata，以后方便使用。
testdata <- GetAssayData(sce, slot="data")
clusters <- sce@meta.data$seurat_clusters
cellpred <- SingleR(test = testdata, ref = ref, labels = ref$label.fine, 
                    # label.fine耗时比较长一点
                    method = "cluster", clusters = clusters, 
                    assay.type.test = "logcounts", assay.type.ref = "logcounts")
save(cellpred,file = "../../cellpred.Rdata")
load("../../cellpred.Rdata")
rm(refdata, HumanPrimaryCellAtlasData, testdata) #珍惜内存
table(cellpred$labels)
celltype = data.frame(ClusterID=rownames(cellpred), celltype=cellpred$labels, stringsAsFactors = F)
table(celltype$ClusterID,celltype$celltype) #如下为singleR的细胞cluster鉴定结果。
#结合上述结果，给scRNA增添celltype注释信息
sce@meta.data$celltype = "NA"
#先新增列celltype，值均为NA，然后利用下一行代码循环填充
for(i in 1:nrow(celltype)){
  sce@meta.data[which(sce@meta.data$seurat_clusters == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}
p1 <- DimPlot(sce, group.by="celltype", label=F , reduction='umap')

ggsave("../../out/3.3celltype_anno.pdf", plot = p1, width = 18, height = 12) 
