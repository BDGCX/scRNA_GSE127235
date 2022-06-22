rm(list = ls())
options(stringsAsFactors = F)
#载入seurat处理后的rds文件
sce <- readRDS(file = "sce_tutorial.rds")

#热图可视化marker基因的表达差异
library(Seurat)
library(pheatmap)
library(ggplot2)
library(scales)
library(dplyr) 
new_clusters <- sce@active.ident
sce@meta.data$new_clusters <- new_clusters
sce.markers <- FindAllMarkers(object = sce, only.pos = TRUE, min.pct = 0.25, 
                              thresh.use = 0.25)#only.pos = TRUE：只返回positive的marker基因
top10 <- sce.markers %>% group_by(cluster) %>% top_n(10, avg_log2FC)

save(sce,sce.markers,top10,file = "sce.RData")

DoHeatmap(sce,top10$gene,size=2,label = T,
          group.by = "new_clusters",angle = 0)+ scale_fill_gradientn(colors = c("blue", "white", "red"))

#使用pheatmap对top10进行重新作图

####方法1#####
# 提取原始表达矩阵
cts <- GetAssayData(sce, slot = "counts")

cts[1:4,1:4]
# 然后对这个矩阵取log对数据进行标准化
cts <- log10(cts +0.001)

#得到小的top10表达矩阵，矩阵的行就按基因名取
#热图结果是按照cluster进行排序展示的，因此我们也要将小表达矩阵的列按cluster从小到大排序：

#原来的cluster分组信息存储在 sce@active.ident 中
head(sce@active.ident)
# 可见并不是按顺序排的
# 排序之后
new_cluster <- sort(sce@active.ident)
head(new_cluster)
#有了行和列的规定，我们就能很轻松地提取出整个小的表达矩阵：
cts <- as.matrix(cts[top10$gene, names(new_cluster)])

#要做出来DoHeatmap的顶部0-8 cluster的展示，需要使用pheatmap的一个参数：annotation_col
#这个参数接收一个数据框作为输入。因为是对列进行注释，所以这个数据框的行是矩阵的列名，而它的列在这里对应的就是cluster分群信息
ac <- data.frame(cluster=new_cluster)
rownames(ac) <- colnames(cts)
table(new_cluster)
#new_cluster
#EC  MC POD  IC  TC 
#379 147  80  68  30

# 列表指定注释行和注释列的颜色,最好与TSNE、UMAP的颜色一致
hue_pal()(5)
#[1] "#F8766D" "#A3A500" "#00BF7D" "#00B0F6" "#E76BF3"
ann_colors <-  list(
  Time = c("white", "firebrick"),
  cluster = c(EC = "#F8766D", MC = "#A3A500",POD="#00BF7D",IC="#00B0F6",TC="#E76BF3"),
  GeneClass = c(Path1 = "#7570B3", Path2 = "#E7298A", Path3 = "#66A61E")
)

pheatmap(cts,show_colnames =F,show_rownames = T,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(nrow(cts)),
         gaps_row = c(10,20,30,40),#设置行之间的空白间隙
         gaps_col = c(379,526,606,674),#设置列之间的间隔
         annotation_colors = ann_colors[2],
         cluster_rows = F,
         cluster_cols = F,
         annotation_col=ac
         )

#进一步使用处理分组treatment进行热图绘制
#提取ac里的处理分组信息，Control或Diabetes
trt <- unlist(lapply(rownames(ac), function(x){
  strsplit(as.character(x),'_')[[1]][4]
}))
ac1 <- ac
ac1$treatment <- trt
ac1 <- ac1[order(ac1[,1],ac1[2]),]

# 提取原始表达矩阵
cts <- GetAssayData(sce, slot = "counts")
cts[1:4,1:4]
# 然后对这个矩阵取log对数据进行标准化
#cts <- scale(cts,center = T,scale = T)
cts <- log10(cts +0.01)

#有了行和列的规定，我们就能很轻松地提取出整个小的表达矩阵：
cts <- as.matrix(cts[top10$gene, rownames(ac1)])

#cts[cts>=4]=2
#cts[cts < (-1)]= -6 #小于负数时，加括号！

ac1 <- as.data.frame(ac1[,-1])
rownames(ac1) <- colnames(cts)
colnames(ac1)[1] <- "Treatment"
ann_colors <-  list(
  Time = c("white", "firebrick"),
  Treatment = c(Control = "#8787FFFF", Diabetes = "#6EE2FFFF"),
  GeneClass = c(Path1 = "#7570B3", Path2 = "#E7298A", Path3 = "#66A61E")
)

pheatmap(cts,show_colnames =F,show_rownames = T,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(nrow(cts)),
         gaps_row = c(10,20,30,40),#设置行之间的空白间隙
         gaps_col = c(379,526,606,674),#设置列之间的间隔
         annotation_colors = ann_colors[2],
         cluster_rows = F,
         cluster_cols = F,
         annotation_col=ac1
)




