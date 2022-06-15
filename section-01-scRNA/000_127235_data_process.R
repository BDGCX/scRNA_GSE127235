#安装需要的包
options()$repos  
#查看当前工作空间默认的下载包路径
options()$BioC_mirror 
#查看使用BioCManager下载包的默认路径
options(BioC_mirror="https://mirrors.ustc.edu.cn/bioc/") 
# 指定使用BioCManager下载的路径
options("repos" = c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/")) 
# 指定使用install.packages下载包的路径
options()$repos 
options()$BioC_mirror


# https://bioconductor.org/packages/release/bioc/html/GEOquery.html
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager") 
#判断是否存在BiocManger包，没有的话下载该包

BiocManager::install("clusterProfiler",ask = F,update = F)
BiocManager::install("org.Mm.eg.db",ask = F,update = F)

### 一.下载、探索、整理数据----
#https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE127235
#点击下方supplementary files的GSE127235_exp1_counts.txt.gz行中间的http下载
#一键清空
rm(list=ls())
options(stringsAsFactors = F)

#读入并处理数据,a1处理同a
a <- read.table("GSE127235_exp1_counts.txt.gz",sep = "\t",header = T)
a[1:4,1:4]
tail(a$Geneid,10)##可以看是否有一些不是基因的行需要剔除

#行名为Geneid，但ENSMUS开头其实是小鼠的Ensembl ID(library(org.Mm.eg.db))，人的Ensembl ID以ENS开头(library(org.Hs.eg.db))，所以需要进行ID转换为Gene Symble
#列名为sample样本名称
summary(a[,2:5]) 
boxplot(a[,2:5])###可以观察到细胞的分布比较集中，只有个别异常值

#将第1列的Geneid转换为基因名称
library(clusterProfiler)
library(org.Mm.eg.db)
keytypes(org.Mm.eg.db)#使用keytypes可以查看注释包的所有 ID 类型
Ensembl <- a$Geneid

#有版本号，直接转不行的, 先去除版本号
Ensembl<- gsub("\\..*", "",  Ensembl)
gene.symbol <- bitr(geneID = Ensembl, 
                    fromType = "ENSEMBL",
                    toType = c("SYMBOL"),#c("ENTREZID", "SYMBOL", "GENENAME")
                    OrgDb = org.Mm.eg.db)

#将a$Geneid去除版本号
a$Geneid <- unlist(lapply(a$Geneid,function(x){
  strsplit(as.character(x),'[.]')[[1]][1]
})
)
colnames(a)[1]<- "ENSEMBL"
tmp <- merge(gene.symbol,a,by='ENSEMBL')
tmp <- tmp[,-1]
#表达量矩阵是ensembl数据库的id格式，然后需要转为基因的名字，这个时候两个id都转为了同样的名字，后续处理就很尴尬。

#所以先对重复基因取均值，然后进行ID转换，也可以直接删除多余的重复基因id或取最大值
a$Geneid <- unlist(lapply(a$Geneid,function(x){
  strsplit(as.character(x),'[.]')[[1]][1]
})
)
tmp <- as.matrix(tmp)
rownames(tmp) <- tmp[,1]
exp <- tmp[,2:ncol(tmp)]
dimnames <- list(rownames(exp),colnames(exp))
data <- matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
library(limma)
data <- avereps(data)#如果基因存在多行中，为了取得结果就会取均值
#以上操作得到id转换的data和data1

dat <- as.data.frame(cbind(data,data1))

dat1 <- dat
#提取样本的属性信息（即metadata）
#在GEO中点击GSE下载数据下方的下方SRA RUN Selector进去点击total那行的Metadata下载
#分析细胞对应的细胞/组织类型和样本来源
b <- read.table("../test-single cell/raw_data/SraRunTable.txt",
                sep = ",", header = T)
b <- format(b,scientific=F)##数字不表示为科学计数
b[1:4,1:4]
c <- as.data.frame(cbind(Sample.Name=b$Sample.Name,treatment=b$treatment,cell_count=b$cell_count))
#GSE127235.csv是从GEO中下载的sample点+more之后Accession list...打开后export得到
d <- read.csv("../test-single cell/raw_data/GSE127235.csv")
colnames(d)[1] <- "Sample.Name"
c <- merge(c,d,by="Sample.Name")
c <- c[,1:4]
colnames(dat1) <- unlist(lapply(colnames(dat1),function(x){
  strsplit(as.character(x),"_S")[[1]][1]
})
)
#如果是"_S"则在_S <- 处拆分，如果是"[_S]",则既在_处拆分，又在S处拆分
dat1 <- as.data.frame(rbind(colnames(dat1),dat1))
dat1 <- as.data.frame(t(dat1))
colnames(dat1)[1] <- "Title"
dat1 <- merge(c,dat1,by="Title")
dat1 <- as.data.frame(cbind(Sample=paste(dat1$Title,dat1$treatment,sep = "_"),dat1))
dat1 <- dat1[,-c(2,3,4,5)]
rownames(dat1) <- dat1$Sample

#筛选单细胞
#Quality control (QC) filters were applied using the following parameters, similar to what has been reported:
#(1) wells with either no cells or more than one cell when checked
#under the microscope were excluded from downstream analysis; 
#文献中(1)去除cell_count等于0或者大于1个的样本
dat1 <- dat1[dat1$cell_count==1,]
dat1 <- dat1[,-c(1,2)]
dat1 <- as.data.frame(t(dat1))
dim(dat1)
#[1] 29369   833

##提取标本组别来源，对照或糖尿病组
colnames(dat1) #取列名
library(stringr)
trt <- str_split(colnames(dat1),'_',simplify = T)[,4] #取列名，以'_'号分割，提取第三列。
#str_split()函数可以分割字符串
table(trt)
#Control Diabetes 
#407      426 

#提取批次信息
colnames(dat1) #取列名
library(stringr)
exp <- str_split(colnames(dat1),'_',simplify = T)[,1] #取列名，以'_'号分割，提取第三列。
#str_split()函数可以分割字符串
table(exp)
#EXP1 EXP2 
#377  456

#统计每个样本有表达的有多少行（基因）
n_g <-  apply(dat1,2,function(x) sum(x>1)) 
hist(n_g,breaks = 30)

##(样本为行名，列分别为：样本处理信息，样本分组批次信息，样本表达的基因数【注意：不是表达量的和，而是种类数或者说个数】)
metadata <- data.frame(treatment=trt,EXP=exp,n_g=n_g) #新建数据框(细胞的属性信息)
metadata$all='all' #添加列，列名为"all"，没啥意思，就是后面有需要

counts <- dat1

#筛选对照组和糖尿病组
index <- grep("Control",colnames(dat1))#查找带有Control字样的元素
ctrl <-dat1[,index] 
DM <- dat1[,-index]

#提取样本属性信息
colnames(DM) #取列名
library(stringr)
exp <- str_split(colnames(DM),'_',simplify = T)[,1] #取列名，以'_'号分割，提取第三列。
#str_split()函数可以分割字符串
table(exp)

#统计每个样本有表达的有多少行（基因）
n_g <-  apply(DM,2,function(x) sum(x>1)) 
hist(n_g,breaks = 30)

metadata_DM <- data.frame(EXP=exp,n_g=n_g)
metadata_DM$all <- "all"

#提取样本属性信息
colnames(ctrl) #取列名
exp <- str_split(colnames(ctrl),'_',simplify = T)[,1] #取列名，以'_'号分割，提取第三列。
#str_split()函数可以分割字符串
table(exp)

#统计每个样本有表达的有多少行（基因）
n_g <-  apply(ctrl,2,function(x) sum(x>1)) 
hist(n_g,breaks = 30)

metadata_ctrl <- data.frame(EXP=exp,n_g=n_g)
metadata_ctrl$all <- "all"

save(dat,counts,metadata,file = 'input.Rdata')
save(ctrl,metadata_ctrl,DM,metadata_DM,file = "treatment.Rdata")
