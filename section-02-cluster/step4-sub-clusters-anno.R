rm(list = ls())
options(stringsAsFactors = F)
library(Seurat)
library(ggplot2)
library(dplyr)
#加载数据
load("sce.IC.subset.Rdata")

#####
#1.手动注释
sce.markers <- FindAllMarkers(object = sce, only.pos = TRUE, min.pct = 0.25, 
                              thresh.use = 0.25)#only.pos = TRUE：只返回positive的marker基因
top10 <- sce.markers %>% group_by(cluster) %>% top_n(10, avg_log2FC)
DoHeatmap(sce,top10$gene,size=2,label = T)
#根据已建立的cell markers进行对比注释

#如何手动注释细胞类型举例
#寻找需要注释疾病细胞的相关文献作为参考，找到文中的marker基因（参考文献中注释的细胞类型
#MCA Macrophage
Macro <-c('C1qa', 'Cd74', 'Cd68','Adgre1',"Chil3","Ccl9","S100a4","Lyz2","Thbs1","Ms4a4c","F10","Ly6c2","Gda","Lgals3")

#MCA Kupffer_cell
KC <- c("Vsig4","Cd5l","Fcna","Cfp","C1qc","Clec4f","Ctsc","Adgre1","Fabp7","C1qa")

#MCA Endothelial_cell
np <- c("Cxcl10","Clec4g","Igfbp7","Adamts1","Plpp3","Iigp1","Kdr","Nrp1","Cyp4b1","Socs3")

#MCA DC
DC <- c("Cst3","Ccr7","H2-Eb1","Ccl22","H2-Aa","Naaa","Gm2a","H2-Ab1","Ppt1","Cytip")

#MCA T_cell
T1 <- c("Gzma","Gzmb","Ccl5","Xcl1","Cd7","Gzmc","Il2rb","Nkg7","Klrb1b","Cd3g")
T2 <- c("Trbc2","Trac","Icos","Satb1","Isy1","Ms4a4b","Cd3d","Trbc1","Sh2d2a","Cd28")

#MCA B_cell
B1 <- c("Bank1", "Bcl11a", "Cd19", "Cd22", "Cd37", "Cd74", "Cd79a", "Cd79b", "Cxcr4", "Ebf1", "Fcer2a", "Fcmr", "Fcrla")
PB <- c("Creld2", "Crip1", "Derl3", "Dnajc3", "Eaf2", "Edem1", "Edem2", "Fam46c", "Glipr1", "Gm43291", "H13", "Herpud1", "Hsp90b1", "Igha", "Igkc", "Iglc2", "Jchain")

#MCA Hepatocyte
Hepa <- c("Alb","Apoa1","Fgb","Gc","Ahsg","Kng1","Mup3","Car3","Gsta3","Hpd","Ass1","Mat1a","Bhmt","Fabp1","Aldob","Wfdc21")

#Cholangiocyte
chol <- c("Alcam", "Ambp", "Ankrd1", "Anxa5", "Atp1b1", "Bicc1", "Ces1d", "Cldn3", "Cldn7", "Clu", "Cp", "Cyr61", "Cystm1", "Dbi", "Ddit4l", "Dsg2")

#HSC
hscc <- c("Angptl6", "Bgn", "C3", "C4b", "Col14a1", "Col1a1", "Col1a2", "Col3a1", "Colec11", "Cxcl12", "Cygb", "Dcn", "Dpt", "Ecm1", "Efemp1", "Gsn", "Ifitm1", "Igfbp5", "Igfbp6")

#Neutrophil
N <- c('S100a8', 'Wfdc21', 'Ly6g')

#B cell
B <- c('Cd79a', 'Mzb1', 'Jchain','Cd19')
#合并
markers <- c(PB,B,T1,T2,DC,Macro,N)

markers <-  c('PTPRC', 'CD3D', 'CD3E','CD4', 'CD8A','FOXP3','KLRD1', ## Tcells
                   'CD19', 'CD79A', 'MS4A1' , # B cells
                   'IGHG1', 'MZB1', 'SDC1',  # plasma 
                   'CD68', 'CD163', 'CD14',  'CD86', 'CCL22','S100A4', ## DC(belong to monocyte)
                   'CD68',  'CD163','MRC1','MSR1' ,'S100A8', ## Macrophage (belong to monocyte)
                   'FGF7', 'MME','COL1A1','ACTA2','PECAM1', 'VWF' ,'PROX1', ## Fibroblasts,Endothelial
                   'EPCAM', 'KRT19','KRT7','KRT8','KRT18', 'PROM1'  ## epi or tumor
)
#G127235文献中的cell markers
markers <- c('C1qa', 'Cd74', 'Cd68','Adgre1',#Macrophage
             'Cd79a', 'Mzb1', 'Jchain','Cd19',#B
             'S100a8', 'Wfdc21', 'Ly6g')#Neutrophil

mac_markers <- c("CCR7","FCGR3A","Nos2", "Il12b", "Tnf",#Macrophage1
                 "CD206","CD16","CD200","CD36","Arg", "Il10", "Clec10a","Mrc1","Cd163")#Macrophage1
#使用CaseMatch去除不存在的基因
markers <- CaseMatch(markers,rownames(sce))
markers <- as.character(markers)
#针对marker gene做点图
DotPlot(sce, features = unique(markers),group.by = "seurat_clusters")+RotatedAxis()+
  scale_x_discrete("")+scale_y_discrete("")
#根据点图对细胞群进行注释
new.cluster.ids <- c("B", "Macro", "Macro","DC","T","Macro","Macro", "Macro")
names(new.cluster.ids) <- levels(sce)
sce <- RenameIdents(sce, new.cluster.ids)

DimPlot(sce,reduction = "tsne",label = T)+
  theme(panel.background = element_rect(fill='white', colour='black'), 
        panel.grid=element_blank()) 

#特征图
markers <- c('Cd74', 'Cd68','Adgre1',#Macrophage
             'Mzb1', 'Jchain',#B
             'Wfdc21', 'Ly6g')

plot_list1 <-  list() 
for (i in 1:7){
  p <- FeaturePlot(object = sce, 
                   features = markers[i],
                   cols=c("lightgrey", "red"),
                   reduction = "tsne",
                   label = F,
                   keep.scale = "all",
                   ncol=1)+NoLegend()+
    theme(axis.title = element_blank(),
          axis.text = element_text(size = 10),
          plot.title = element_text(size = 10))
  
  plot_list1[[i]] = p
}

markers <-  c('C1qa','Cd79a', 'S100a8')
ylabs <- c("Macrophage","B cell","Neutrophil")
plot_list2 <- list()
for (i in 1:3){
  p <- FeaturePlot(object = sce, 
                   features = markers[i],
                   cols=c("lightgrey", "red"),
                   reduction = "tsne",
                   label = F,
                   keep.scale = "all",
                   ncol=1)+NoLegend()+
    theme(axis.title.x =element_blank(),
          axis.title.y = element_text(size = 10,face = "bold"),
          axis.text = element_text(size = 10),
          plot.title = element_text(size = 10))+
    ylab(label = ylabs[i])
  
  plot_list2[[i]] = p
}
library(patchwork)
#+ 是按列排，/是按列排,但/优先,|按行
 (plot_list2[[1]]|plot_list1[[1]]|plot_list1[[2]]|plot_list1[[3]])/
  (plot_list2[[2]]|plot_list1[[4]]|plot_list1[[5]]|plot_spacer())/
  (plot_list2[[3]]|plot_list1[[6]]|plot_list1[[7]]|plot_spacer())#plot_spacer()空白图

ggsave("IC-featureplot.png",width = 2500,height = 1980,
       dpi = 300,units = "px")

#####
# 2.SingleR 做 cell annotation的流程(准确度较低)
#singleR自带的7个参考数据集，需要联网才能下载，其中5个是人类数据，2个是小鼠的数据：
#BlueprintEncodeData Blueprint (Martens and Stunnenberg 2013) and Encode (The ENCODE Project Consortium 2012) （人）
#DatabaseImmuneCellExpressionData The Database for Immune Cell Expression(/eQTLs/Epigenomics)(Schmiedel et al. 2018)（人）
#HumanPrimaryCellAtlasData the Human Primary Cell Atlas (Mabbott et al. 2013)（人）
#MonacoImmuneData, Monaco Immune Cell Data - GSE107011 (Monaco et al. 2019)（人）
#NovershternHematopoieticData Novershtern Hematopoietic Cell Data - GSE24759（人）
#ImmGenData the murine ImmGen (Heng et al. 2008) （鼠）
#MouseRNAseqData a collection of mouse data sets downloaded from GEO (Benayoun et al. 2019).鼠）

library(SingleR)
#remotes::install_github("LTLA/celldex")
library(celldex)
sce_for_SingleR <- GetAssayData(sce, slot="data")
clusters <- sce@meta.data$seurat_clusters

#mouseImmu <- ImmGenData()
mouseImmu <- load("../SingleR_ref/ref_Mouse_imm.RData")

pred.mouseImmu <- SingleR(test = sce_for_SingleR, ref = mouseImmu, 
                     labels = mouseImmu$label.main,
                     #因为样本主要为免疫细胞（而不是全部细胞），因此设置为label.fine
                     method = "cluster", clusters = clusters,
                     #这里我们为上一步分的8个cluster注释celltype
                     assay.type.test = "logcounts", assay.type.ref = "logcounts")
table(pred.mouseImmu$labels)

#mouseRNA <- MouseRNAseqData()
mouseRNA <- load("../SingleR_ref/ref_Mouse_all.RData")
pred.mouseRNA <- SingleR(test = sce_for_SingleR, ref = mouseRNA, labels = mouseRNA$label.fine ,#其它参数label.main,label.fine耗时比较长一点
                         method = "cluster", clusters = clusters, 
                         assay.type.test = "logcounts", assay.type.ref = "logcounts")

cellType <- data.frame(ClusterID=levels(sce@meta.data$seurat_clusters),
                    mouseImmu=pred.mouseImmu$labels,
                    mouseRNA=pred.mouseRNA$labels)
head(cellType)
sce@meta.data$singleR=cellType[match(clusters,cellType$ClusterID),'mouseImmu']
sce@meta.data$singleR1=cellType[match(clusters,cellType$ClusterID),'mouseImmu']

save(sce,file = "./sce.IC.subset.Rdata")
DimPlot(sce,reduction = "umap",label=T, group.by = 'singleR1')
#绘制singleR scores热图
plotScoreHeatmap(pred.mouseImmu,
                 show.labels = F,
                 cluster_cols = T,
                 color = colorRampPalette(c("blue","yellow","red"))(50))

plotDeltaDistribution(pred.mouseImmu, ncol = 3)



