#https://github.com/wu-yc/scMetabolism

#一键清空
rm(list = ls())
options(stringsAsFactors = F)

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

#安装所需的包
install.packages(c("devtools", "data.table", 
                   "wesanderson", "Seurat", "devtools",
                   "AUCell", "GSEABase", "GSVA",
                   "ggplot2","rsvd"))
BiocManager::install("AUCell",ask = F,update = F)
devtools::install_github("YosefLab/VISION")
devtools::install_github("wu-yc/scMetabolism")
#加载包
library(scMetabolism)
library(ggplot2)
library(rsvd)
library(Seurat)
library(dplyr)

#读入数据
countexp.Seurat <- readRDS("sce_tutorial.rds") 
DimPlot(countexp.Seurat)

#运行scMetabolism计算代谢通路活性
## Quantify single-cell metabolism with Seurat (Recommended)
countexp.Seurat<-sc.metabolism.Seurat(obj = countexp.Seurat,
                                      method = "VISION",
                                      imputation = F,
                                      ncores = 2,
                                      metabolism.type = "KEGG")

#obj is a Seurat object containing the UMI count matrix.

#method supports VISION, AUCell, ssgsea, and gsva, which VISION is the default method.

#imputation allows users to choose whether impute their data before metabolism scoring.

#ncores is the number of threads of parallel computation.

#metabolism.type supports KEGG and REACTOME, where KEGG contains 85 metabolism pathways and REACTOME contains 82 metabolism pathways.

#To extract the metabolism score, just run 
metabolism.matrix <- countexp.Seurat@assays$METABOLISM$score

metabolism.matrix
#可视化,只能可视化metabolism.matrix里有的通路
DimPlot.metabolism(obj = countexp.Seurat, 
                   pathway = "Oxidative phosphorylation", 
                   dimention.reduction.type = "umap", 
                   dimention.reduction.run = F, 
                   size = 1)

input.pathway<-c("Glycolysis / Gluconeogenesis", "Oxidative phosphorylation", "Citrate cycle (TCA cycle)")
input.pathway<-c( "Oxidative phosphorylation")

DotPlot.metabolism(obj = countexp.Seurat, 
                   pathway = input.pathway, 
                   phenotype = "new_clusters", 
                   norm = "y")

BoxPlot.metabolism(obj = countexp.Seurat, 
                   pathway = input.pathway, 
                   phenotype = "new_clusters", 
                   ncol = 1)


#载入改造的sc.metabolism.Seurat.2函数
source("code/section-07-scMetabolism/function_sc_metabolism_Seurat.R")
#运行上述写好的sc.metabolism.Seurat.2函数
countexp.Seurat<-sc.metabolism.Seurat.2(obj = countexp.Seurat, 
                                        method = "VISION", 
                                        imputation = F,
                                        ncores = 1,
                                        geneList = "./h.all.v7.5.1.symbols.gmt")

input.pathway<-c("HALLMARK_TNFA_SIGNALING_VIA_NFKB",
                 "HALLMARK_HYPOXIA",
                 "HALLMARK_CHOLESTEROL_HOMEOSTASIS")

DotPlot.metabolism(obj = countexp.Seurat,
                   pathway = input.pathway,
                   phenotype = "new_clusters",
                   norm = "y")

#改造相关可视化函数
source("code/section-07-scMetabolism/function_DotPlot_metabolism2.R")

#进行x和y轴的因子水平设置
### y轴
input.pathway<-c("HALLMARK_COMPLEMENT",
                 "HALLMARK_COAGULATION",
                 "HALLMARK_CHOLESTEROL_HOMEOSTASIS") %>%
  factor(levels = c("HALLMARK_COMPLEMENT",
                    "HALLMARK_CHOLESTEROL_HOMEOSTASIS",
                    "HALLMARK_COAGULATION"))

### x轴
countexp.Seurat$new_clusters = factor(countexp.Seurat$new_clusters,
                               levels = c("EC", 
                                          "MC", 
                                          "POD", 
                                          "IC", 
                                          "TC"))
DotPlot.metabolism.2(obj = countexp.Seurat,
                    pathway = input.pathway,
                    phenotype = "new_clusters",
                    norm = "y")







