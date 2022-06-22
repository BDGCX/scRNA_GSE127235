rm(list = ls())
#加载包
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)

#加载数据
sce <- readRDS("./sce_tutorial.rds")

?FindMarkers #ident.1是病例组，ident.2是对照组
deg <- FindMarkers(object = sce, 
                   ident.1 = "EC",
                   test.use='MAST' )

gene_up <- rownames(deg[deg$avg_log2FC>1,])
gene_down <- rownames(deg[deg$avg_log2FC < -1,])
gene_total <- rownames(deg[deg$avg_log2FC>1|deg$avg_log2FC < -1,])

library(org.Mm.eg.db)
#把SYMBOL改为ENTREZID，这里会损失一部分无法匹配到的
gene_up <- as.character(na.omit(AnnotationDbi::select(org.Mm.eg.db,
                                                      keys = gene_up,
                                                      columns = 'ENTREZID',
                                                      keytype = 'SYMBOL')[,2]))
gene_down <- as.character(na.omit(AnnotationDbi::select(org.Mm.eg.db,
                                                        keys = gene_down,
                                                        columns = 'ENTREZID',
                                                        keytype = 'SYMBOL')[,2]))

# 人家的包就是 entrez ID  设计
library(clusterProfiler)
#1.function里是kegg和GO以及GSEA分析
source('./code/section-02-cluster/functions.R')
#直接利用之前保存的脚本中的代码对正、负相关基因进行kegg富集分析
run_kegg(gene_up,gene_down,pro='test')
# 下面的 run_go 会比较慢，根据需要
# run_go(gene_up,gene_down,pro='shRNA')

#2.不使用function
gene_up <- unique(gene_up)
kk.up <- enrichKEGG(gene         = gene_up,
                    organism     = 'mmu',
                    #universe   = gene_all,
                    pvalueCutoff = 0.9,
                    qvalueCutoff =0.9)
head(kk.up)[,1:6]

dotplot(kk.up)

kegg_up_dt <- as.data.frame(kk.up)
kegg_up_dt <- kegg_up_dt[order(kegg_up_dt$p.adjust),]
kegg_up_dt_top10 <- kegg_up_dt[1:10,]

#排序
kegg_up_dt_top10 <- kegg_up_dt_top10[order(kegg_up_dt_top10$Count,decreasing = F),]
kegg_up_dt_top10$Description <- factor(kegg_up_dt_top10$Description,levels = kegg_up_dt_top10$Description)


ggplot(data = kegg_up_dt_top10,aes(x=Description,y=Count))+
  geom_point(stat="identity",aes(size=Count,color=p.adjust))+
  scale_color_gradient(low="red",high="blue")+
  coord_flip()+#xy轴置换
  xlab(" ")+
  ylab("")+
  theme(axis.text.x = element_text(colour = "black"))+
  theme_light()+
  theme(legend.key.size = unit(8, "pt"))#修改图例legend大小
  

ggsave("kegg.png",width = 2500,height = 1200,
       dpi = 300,units = "px")

#对正相关基因进行富集分析
#go <- enrichGO(gene_total, OrgDb = "org.Mm.eg.db", ont="all") #ont可设置只计算BP,CC或MF

library(stringr)
go_bp <- enrichGO(gene_up, OrgDb = "org.Mm.eg.db", ont="BP") #ont可设置只计算BP,CC或MF

go_bp <- go_bp[order(go_bp@result$p.adjust)]
go_bp_top10 <- go_bp[1:10,]

#排序
go_bp_top10 <- go_bp_top10[order(go_bp_top10$Count,decreasing = F),]
go_bp_top10$Description <- factor(go_bp_top10$Description,levels = go_bp_top10$Description)

labels <- as.character(paste(go_bp_top10$Description,"(",go_bp_top10$ID,")"))

ggplot(data = go_bp_top10,aes(x=Description,y=Count))+
  geom_bar(stat = "identity",width = 0.8,fill="grey",colour="black")+
  coord_flip()+#xy轴置换
  xlab(" ")+
  ylab("-log(p-value)")+
  scale_x_discrete(label=labels)+
  scale_y_continuous(expand = c(0, 0))+
  theme(axis.text.x = element_text(colour = "black"))+
  theme_classic()

ggsave("go.png",width = 2500,height = 1200,
       dpi = 300,units = "px")
#BP,CC,MF都画
barplot(go, split="Description")+ facet_grid(ONTOLOGY~., scale="free") 






