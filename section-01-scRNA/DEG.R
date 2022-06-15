#一键清空
rm(list=ls())
options(stringsAsFactors = F)
#导入数据
load("127235.Rdata")
#数据格式转换
total <- t(total)
colnames(total) <- total[1,]
total <- total[-c(1,2),]
str(total)
name <- rownames(total)
total <- apply(total,2,as.numeric)
#2代表修改的为列，1代表修改的为行
rownames(total) <- name
str(total)
#设置分组
group_list <- unlist(lapply(colnames(total),function(x){
  strsplit(as.character(x),'_')[[1]][2]
})
)
table(group_list)
library(limma)
library(edgeR)
library(ggplot2)
###画图函数
draw_h_v <- function(exprSet,need_DEG,n='DEseq2',group_list,logFC_cutoff){
  ## we only need two columns of DEG, which are log2FoldChange and pvalue
  ## heatmap
  library(pheatmap)
  choose_gene=head(rownames(need_DEG),50) ## 50 maybe better
  choose_matrix=exprSet[choose_gene,]
  choose_matrix[1:4,1:4]
  choose_matrix=t(scale(t(log2(choose_matrix+1)))) 
  ## http://www.bio-info-trainee.com/1980.html
  annotation_col = data.frame( group_list=group_list  )
  rownames(annotation_col)=colnames(exprSet)
  pheatmap(choose_matrix,show_colnames = F,annotation_col = annotation_col,
           filename = paste0(n,'_need_DEG_top50_heatmap.png'))
  
  ###PCA图
  library(ggfortify)
  df=as.data.frame(t(choose_matrix))
  df$group=group_list
  png(paste0(n,'_DEG_top50_pca.png'),res=120)
  p=autoplot(prcomp( df[,1:(ncol(df)-1)] ), data=df,colour = 'group')+theme_bw()
  print(p)
  dev.off()
  
  ####火山图
  if(! logFC_cutoff){
    logFC_cutoff <- with(need_DEG,mean(abs( log2FoldChange)) + 2*sd(abs( log2FoldChange)) )
    
  }
  # logFC_cutoff=1
  
  need_DEG$change = as.factor(ifelse(need_DEG$pvalue < 0.05 & abs(need_DEG$log2FoldChange) > logFC_cutoff,
                                     ifelse(need_DEG$log2FoldChange > logFC_cutoff ,'UP','DOWN'),'NOT')
  )
  this_tile <- paste0('Cutoff for logFC is ',round(logFC_cutoff,3),
                      '\nThe number of up gene is ',nrow(need_DEG[need_DEG$change =='UP',]) ,
                      '\nThe number of down gene is ',nrow(need_DEG[need_DEG$change =='DOWN',])
  )
  library(ggplot2)
  g = ggplot(data=need_DEG, 
             aes(x=log2FoldChange, y=-log10(pvalue), 
                 color=change)) +
    geom_point(alpha=0.4, size=1.75) +
    theme_set(theme_set(theme_bw(base_size=20)))+
    xlab("log2 fold change") + ylab("-log10 p-value") +
    ggtitle( this_tile ) + theme(plot.title = element_text(size=15,hjust = 0.5))+
    scale_colour_manual(values = c('blue','black','red')) ## corresponding to the levels(res$change)
  print(g)
  ggsave(g,filename = paste0(n,'_volcano.png'))
  dev.off()
}

###第一种方法limma包
if(T){
  suppressMessages(library(limma))
  design <- model.matrix(~0+factor(group_list))
  colnames(design)=levels(factor(group_list))
  rownames(design)=colnames(total)
  design
  
  dge <- DGEList(counts=total)
  dge <- calcNormFactors(dge)
  logCPM <- cpm(dge, log=TRUE, prior.count=3)
  
  v <- voom(dge,design,plot=TRUE, normalize="quantile")
  fit <- lmFit(v, design)
  
  group_list
  cont.matrix=makeContrasts(contrasts=c('Diabetes-Control'),levels = design)
  fit2=contrasts.fit(fit,cont.matrix)
  fit2=eBayes(fit2)
  
  tempOutput = topTable(fit2, coef='Diabetes-Control', n=Inf)
  DEG_limma_voom = na.omit(tempOutput)
  head(DEG_limma_voom)
  nrDEG=DEG_limma_voom[,c(1,4)]
  colnames(nrDEG)=c('log2FoldChange','pvalue') 
  draw_h_v(total,nrDEG,'limma',group_list,1)
}

tmp_f=file.path(Rdata_dir,'TOTAL-RNA-DEG_results.Rdata')

if(file.exists(tmp_f)){
  save(DEG_limma_voom, file = tmp_f)
  
}else{
  load(file = tmp_f) 
}
D <- nrDEG[which(abs(nrDEG$log2FoldChange)>=1),]
View(D)
RP <- read.csv(file = "./RP.csv",header = T)
View(RP)
a <- rownames(D)
D1 <- as.matrix(cbind(gene.name=a,D))
View(D1)
merge <- merge(D1,RP,by="gene.name")


nrDEG1=DEG_limma_voom[,c(1,4)]
colnames(nrDEG1)=c('log2FoldChange','pvalue') 