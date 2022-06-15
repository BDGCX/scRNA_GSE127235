#5.提取各个细胞类型的marker gene
#利用 Finscearkers 命令，可以找到找到各个细胞类型中与其他类别的差异表达基因，作为该细胞类型的生物学标记基因。
#其中ident.1参数设置待分析的细胞类别，min.pct表示该基因表达数目占该类细胞总数的比例
#5.1 网络教程
#5.1.1 此步可省略，对marker基因进行可视化，这一步会把所有cluster的marker基因展示出来。这步骤运行时间较久
pro='first'
table(sce@meta.data$seurat_clusters) 
for( i in unique(sce@meta.data$seurat_clusters) ){
  markers_df <- FindMarkers(object = sce, ident.1 = i, min.pct = 0.25,test.use = "wilcox")
  print(x = head(markers_df))
  markers_genes =  rownames(head(x = markers_df, n = 5))
  VlnPlot(object = sce, features =markers_genes,log =T )
  ggsave(filename=paste0(pro,'_VlnPlot_subcluster_',i,'_sce.markers_heatmap.pdf'))
  FeaturePlot(object = sce, features=markers_genes )
  ggsave(filename=paste0(pro,'_FeaturePlot_subcluster_',i,'_sce.markers_heatmap.pdf'))
}

#提取marker基因，min.ct：设定研究的基因在两组细胞中的最小占比。
#Each cluster was screened for marker genes by differential expression analysis (DEA) between cells inside and outside the cluster using Finscearkers function with parameters min.pct=0.25 (genes expressed in at least 25% of cells either inside or outside ofa cluster) and test.use=“wilcox” 
sce.markers <- FindAllMarkers(object = sce, only.pos = TRUE, min.pct = 0.25, 
                              thresh.use = 0.25)#only.pos = TRUE：只返回positive的marker基因
print(x = head(sce.markers))

DT::datatable(sce.markers)
write.csv(sce.markers,file=paste0('_sce.markers_tsne.csv'))

#5.1.2探索感兴趣的基因
#可视化marker基因,我们选择ALB基因可视化看看
#应激相关基因
stress <- c("Fos",#"Scos3",
            "Jun","Junb","Atf3","Egr1",
            "Cebpd","Hspa1b","Hsp90aa1")
markers_genes =   rownames(sce.markers[c("Ptprb","Cyyr1","Egfl7"),])
#小提琴图
VlnPlot(object = sce, features =stress,log=T,ncol=3)
#山峦图
features <- c("Ptprb","Sfrp2","Clic3")
RidgePlot(sce, features = features, ncol = 1,stack = T)+
  theme(axis.title.x.bottom = element_text(size=10),
        axis.title.y = element_blank(),
        strip.text.x = element_text(angle = 0))
ggsave("ridgeplot.png",width = 1500,height = 1000,dpi = 300,units = "px")
#点图
features <- c("Ptprb","Cyyr1","Sfrp2","Agtr1a","Clic3","Cdkn1c")
DotPlot(sce, features = features) + RotatedAxis()+
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(hjust = 0.5,angle = 0,size=9),
        axis.text.y = element_text(size=9),
        strip.text.x = element_text(angle = 0))+
  guides(color=guide_colorbar("Average Expression",
                              barwidth = 0.5,
                              label.theme = element_text(size = 7),
                              title.theme = element_text(size=9)),
         size=guide_legend('Percent Expression',
                           label.theme = element_text(size = 7),
                           title.theme = element_text(size=9)),
         shape="none")

ggsave("dotplot.png",width = 1700,height = 1000,dpi = 300,units = "px")

#特征图
plot_list <-  list() 

for (i in 1:7){
  p <- FeaturePlot(object = sce, 
                   features = stress[i],
                   cols=c("lightgrey", "red"),
                   reduction = "tsne",
                   label = F,
                   keep.scale = "all",
                   ncol=1)+NoLegend()+
    theme(axis.title = element_text(size = 10),
          plot.title = element_text(size = 10))
  
  plot_list[[i]] = p
}
library(patchwork)
#+ 是按列排，/是按列排,但/优先,|按行
  (plot_list[[1]]|plot_list[[2]]|plot_list[[3]])/
  (plot_list[[4]]|plot_list[[5]]|plot_list[[6]])/
  (plot_list[[7]]|p8 |plot_spacer())#plot_spacer()空白图

p8 <- FeaturePlot(object = sce, 
            features ="Hsp90aa1",
            cols=c("lightgrey", "red"),
            reduction = "tsne",
            keep.scale = "all",
            ncol=1)+
  guides(color=(guide_colorbar("Expression ",
                title.position="top",
                title.theme = element_text(size = 8),
                barheight = 0.5,
                barwidth = 3,
                direction = "horizontal",
                label.position ="bottom",
                label.theme = element_text(size = 6)
            )))+
  theme(legend.position =c(1.2,-0.3),
        axis.title = element_text(size=10),
        plot.title = element_text(size = 10))#更改图片标题
 #scale_fill_discrete(labels=c("Low", "High"))更改图例标签名字
  
ggsave("feature.png",width = 1980,height = 1980,dpi = 300,units = "px")      
#热图可视化marker基因的表达差异
library(dplyr) 
top10 <- sce.markers %>% group_by(cluster) %>% top_n(10, avg_log2FC)
DoHeatmap(sce,top10$gene,size=2,label = T)
ggsave(filename=paste0('_sce.markers_——stne_heatmap.pdf'))
save(sce,sce.markers,file = 'last_sce_tsne.Rdata')

#5.2 Seurat官网教程
#5.2.1find all markers of cluster 1
cluster1.markers <- Finscearkers(sce, ident.1 = 1, min.pct = 0.25)
head(cluster1.markers, n = 5)
#利用 DoHeatmap 命令可以可视化marker基因的表达
sce.markers <- FindAllMarkers(sce, only.pos = TRUE, min.pct = 0.25)
?Finscearkers
library(dplyr)
top3 <- sce.markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)
DoHeatmap(sce, features = top3$gene) + NoLegend()

#5.2.2.探索感兴趣的基因
VlnPlot(sce, features = c("Mrpl12", "Nfe2l2","Bloc1s1","Sqstm1"))
#我们能够看到，MS4A1和CD79A两个基因在细胞群体3中特异性表达。
#you can plot raw counts as well
VlnPlot(sce, features = c("Mrpl12", "Nfe2l2"), slot = "counts", log = TRUE)

FeaturePlot(sce, features = c("Mrpl12", "Nfe2l2"))
#这种展示方法把基因表达量映射到UMAP结果中，同样可以直观的看到基因表达的特异性。
save(sce,sce.markers,file = 'last_sce_UAMP.Rdata')
