rm(list = ls())

library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)

#加载数据
sce <- readRDS("./sce_tutorial.rds")

deg <- FindMarkers(object = sce, 
                   ident.1 = "EC",
                   test.use='MAST' )

head(deg)
#得到差异基因列表后取出 ，p值和logFC
nrDEG=deg[,c(2,1)]
colnames(nrDEG)=c('log2FoldChange','pvalue')
head(nrDEG) 
library(org.Mm.eg.db)
#Y叔的R包，把SYMBOL转换为ENTREZID，后面可以直接做 KEGG 和 GO
library(clusterProfiler)
gene <- bitr(rownames(nrDEG),     #转换的列是nrDEG的列名
             fromType = "SYMBOL",     #需要转换ID类型
             toType =  "ENTREZID",    #转换成的ID类型
             OrgDb = org.Mm.eg.db)    #对应的物种，小鼠的是org.Mm.eg.db
#让基因名、ENTREZID、logFC对应起来
gene$logFC <- nrDEG$log2FoldChange[match(gene$SYMBOL,rownames(nrDEG))]
#match返回第一个向量中的元素在第二个向量中的位置坐标
#match(c(1, "TRUE"), c(T, 0, "1"))
head(gene)
geneList=gene$logFC
names(geneList)=gene$ENTREZID 
#按照logFC的值来排序geneList
geneList=sort(geneList,decreasing = T)
head(geneList)

#clusterProfiler包的妙用
kk_gse <- gseKEGG(geneList     = geneList,
                  organism     = 'mmu',#小鼠为mmu,人为hsa
                  nPerm        = 1000,
                  minGSSize    = 10,
                  pvalueCutoff = 0.9,
                  verbose      = FALSE)
#提取结果
tmp=kk_gse@result
kk=DOSE::setReadable(kk_gse, OrgDb='org.Mm.eg.db',keyType='ENTREZID')
tmp=kk@result
pro='test_gsea'
write.csv(kk@result,paste0(pro,'_kegg.gsea.csv'))

if(F){
  #按照enrichment score从高到低排序，便于查看富集通路
  kk_gse=kk
  sortkk<-kk_gse[order(kk_gse$enrichmentScore, decreasing = T),]
  head(sortkk)
  dim(sortkk)
  #write.table(sortkk,"gsea_output2.txt",sep = "\\t",quote = F,col.names = T,row.names = F)
  #可以根据自己想要的通路画出需要的图
  library(enrichplot)
  gseaplot(kk_gse, "mmu04612")
  gseaplot2(kk_gse, "mmu04612", color = "firebrick",
            rel_heights=c(1, .2, .6))#改变更多参数，为了美观
  
  #同时展示多个pathways的结果。
  #画出排名前四的通路
  gseaplot2(kk_gse, row.names(sortkk)[1:4])
  
  #上图用的是ES排名前4个画图，也可以用你自己感兴趣的通路画图
  # 自己在刚才保存的txt文件里挑选就行。
  paths <- c("mmu04974", "mmu04612", "mmu04658")
  paths <- row.names(sortkk)[5:8]
  paths
  gseaplot2(kk_gse, paths)
  
  #这里的GSEA分析其实由三个图构成，GSEA分析的runningES折线图
  # 还有下面基因的位置图，还有所谓的蝴蝶图。如果不想同时展示，还可以通过subplots改变。
  gseaplot2(kk_gse, paths, subplots=1)#只要第一个图
  gseaplot2(kk_gse, paths, subplots=1:2)#只要第一和第二个图
  gseaplot2(kk_gse, paths, subplots=c(1,3))#只要第一和第三个图
  
  #如果想把p值标上去，也是可以的。
  gseaplot2(kk_gse, paths, color = colorspace::rainbow_hcl(4),
            pvalue_table = TRUE)
  
  #最后的总结代码
  gseaplot2(kk_gse,#数据
            row.names(sortkk)[39],#画那一列的信号通路
            title = "Prion disease",#标题
            base_size = 15,#字体大小
            color = "green",#线条的颜色
            pvalue_table = TRUE,#加不加p值
            ES_geom="line")#是用线，还是用d点
  
}

###当然，这里不知道具体需要什么通路，就全部都画出来
# 这里找不到显著下调的通路，可以选择调整阈值，或者其它。
down_kegg<-kk_gse[kk_gse$pvalue<0.1 & kk_gse$enrichmentScore < -0.4,];down_kegg$group=-1
up_kegg<-kk_gse[kk_gse$pvalue<0.1 & kk_gse$enrichmentScore > 0.4,];up_kegg$group=1
dat=rbind(up_kegg,down_kegg)
colnames(dat)
dat$pvalue = -log10(dat$pvalue)
dat$pvalue=dat$pvalue*dat$group 
dat=dat[order(dat$pvalue,decreasing = F),]
library(ggplot2)
#条形图
g_kegg<- ggplot(dat, aes(x=reorder(Description,order(pvalue, decreasing = F)), y=pvalue, fill=group)) + 
  geom_bar(stat="identity") + 
  scale_fill_gradient(low="blue",high="red",guide = FALSE) + 
  scale_x_discrete(name ="Pathway names") +
  scale_y_continuous(name ="log10P-value") +
  coord_flip() + theme_bw(base_size = 15)+
  theme(plot.title = element_text(hjust = 0.5),  axis.text.y = element_text(size = 15))+
  ggtitle("Pathway Enrichment") 
g_kegg
print(g_kegg)
ggsave(g_kegg,filename = paste0(pro,'_kegg_gsea.pdf'))

library(enrichplot)
gesa_res=kk@result

###画出每张kegg通路
lapply(1:nrow(down_kegg), function(i){ 
  gseaplot2(kk,down_kegg$ID[i],
            title=down_kegg$Description[i],pvalue_table = T)
  ggsave(paste0(pro,'_down_kegg_',
                gsub('/','-',down_kegg$Description[i])
                ,'.pdf'))
})
lapply(1:nrow(up_kegg), function(i){ 
  gseaplot2(kk,up_kegg$ID[i],
            title=up_kegg$Description[i],pvalue_table = T)
  ggsave(paste0(pro,'_up_kegg_',
                gsub('/','-',up_kegg$Description[i]),
                '.pdf'))
})



ego <- gseGO(geneList     = geneList,
             OrgDb        = org.Hs.eg.db,
             ont          = "ALL",
             nPerm        = 1000,   ## 排列数
             minGSSize    = 100,
             maxGSSize    = 500,
             pvalueCutoff = 0.9,
             verbose      = FALSE)  ## 不输出结果
go=ego@result

write.csv(go,file = 'gse_go.csv') 


