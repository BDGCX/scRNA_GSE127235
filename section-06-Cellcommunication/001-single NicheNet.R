#一键清空
rm(list=ls())
options(stringsAsFactors = F)

#载入包
#devtools::install_github("saeyslab/nichenetr")
library(nichenetr)
library(Seurat) # please update to Seurat V4
library(tidyverse)
suppressPackageStartupMessages(library(circlize))

#Read in the expression data of interacting cells
sce <- readRDS("sce_tutorial.rds")
#Read in NicheNet’s ligand-target prior model(先验模型), ligand-receptor network and weighted integrated networks(加权整合网络)
ligand_target_matrix <-  readRDS(url("https://zenodo.org/record/3260758/files/ligand_target_matrix.rds"))
#如果上述网址加载不了就使用网址下载后读入
ligand_target_matrix <- readRDS("ligand_target_matrix.rds")

ligand_target_matrix[1:5,1:5] # target genes in rows, ligands in columns
##                 CXCL1        CXCL2        CXCL3        CXCL5         PPBP
## A1BG     3.534343e-04 4.041324e-04 3.729920e-04 3.080640e-04 2.628388e-04
## A1BG-AS1 1.650894e-04 1.509213e-04 1.583594e-04 1.317253e-04 1.231819e-04
## A1CF     5.787175e-04 4.596295e-04 3.895907e-04 3.293275e-04 3.211944e-04
## A2M      6.027058e-04 5.996617e-04 5.164365e-04 4.517236e-04 4.590521e-04
## A2M-AS1  8.898724e-05 8.243341e-05 7.484018e-05 4.912514e-05 5.120439e-05

lr_network <-  readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds"))
head(lr_network)
## # A tibble: 6 x 4
##   from  to    source         database
##   <chr> <chr> <chr>          <chr>   
## 1 CXCL1 CXCR2 kegg_cytokines kegg    
## 2 CXCL2 CXCR2 kegg_cytokines kegg    

weighted_networks <-  readRDS(url("https://zenodo.org/record/3260758/files/weighted_networks.rds"))
head(weighted_networks$lr_sig) # interactions and their weights in the ligand-receptor + signaling network
## # A tibble: 6 x 3
##   from  to     weight
##   <chr> <chr>   <dbl>
## 1 A1BG  ABCC6  0.422 
## 2 A1BG  ACE2   0.101 

head(weighted_networks$gr) # interactions and their weights in the gene regulatory network

#Perform the NicheNet analysis

# indicated cell types should be cell class identities
# check via: 
# seuratObj %>% Idents() %>% table()
nichenet_output <-  nichenet_seuratobj_aggregate(
  seurat_obj = sce, 
  receiver = "EC", 
  condition_colname = "treatment", 
  condition_oi = "Diabetes", 
  condition_reference = "Control", 
  sender = c("MC","POD"), #sender="all" "undefined"
  ligand_target_matrix = ligand_target_matrix, 
  lr_network = lr_network, 
  weighted_networks = weighted_networks, 
  organism = "mouse")
## [1] "Read in and process NicheNet's networks"
## [1] "Define expressed ligands and receptors in receiver and sender cells"
## [1] "Perform DE analysis in receiver cell"
## [1] "Perform NicheNet ligand activity analysis"
## [1] "Infer active target genes of the prioritized ligands"
## [1] "Infer receptors of the prioritized ligands"
## [1] "Perform DE analysis in sender cells"

#print a list
nichenet_output %>% names()
view(nichenet_output)

#Interpret the NicheNet analysis output

#1.Ligand activity analysis results

#The column ‘bona_fide_ligand’ indicates whether the ligand is part of ligand-receptor interactions that are documented in public databases
#NicheNet ranks the ligands based on their pearson correlation coefficient.
nichenet_output$ligand_activities
#To get a list of the 20 top-ranked ligands
nichenet_output$top_ligands
#To see which cell population expresses which of these top-ranked ligands
nichenet_output$ligand_expression_dotplot

DotPlot(sce,features = nichenet_output$top_ligands %>% rev(),
        cols = "RdYlBu")+
  RotatedAxis()

VlnPlot(sce , 
        features = nichenet_output$top_ligands, 
        split.by = "treatment",    
        pt.size = 0, 
        combine = T)
#to see whether some of these ligands are differentially expressed from different group(不同组别对照和DM组配体的差异表达)
nichenet_output$ligand_differential_expression_heatmap#这是个ggplot对象，可以使用ggplot2进行优化

DotPlot(sce, features = nichenet_output$top_ligands %>% rev(), 
        split.by = "treatment")+
  RotatedAxis()
#2.Inferred active ligand-target links
nichenet_output$ligand_target_heatmap
#这是个ggplot对象，可以使用ggplot2进行优化
nichenet_output$ligand_target_heatmap + 
  scale_fill_gradient2(low = "whitesmoke",  
                       high = "royalblue", 
                       breaks = c(0,0.0045,0.009)) + 
  xlab("Predict target genes in send cells") + 
  ylab("Prioritized  cell ligands")

#extract the ligand-target links and their regulatory potential scores in matrix or data frame format
nichenet_output$ligand_target_matrix %>% .[1:10,1:6]

nichenet_output$ligand_target_df # weight column = regulatory potential

#To get a list of the top-predicted target genes of the 20 top-ranked ligands
nichenet_output$top_targets
#visualize the expression of these target genes(receiver cells-EC) and observe whether they are differentially expressed between groups
DotPlot(sce%>% subset(idents = "EC"), features = nichenet_output$top_targets %>% rev(), 
        split.by = "treatment")+
  RotatedAxis()

VlnPlot(sce %>% subset(idents = "EC"), 
        features = c("Vegfa","Csf1"), 
        split.by = "treatment",    
        pt.size = 0,
        #split.plot = TRUE,
        combine = FALSE)

#To visualize ligand activities, expression, differential expression and target genes of ligands

#important: above figure can be considered as one of the most important summary figures of the NicheNet analysis. 
#Here you can see which ligand-receptor pairs have both high differential expression and ligand activity (=target gene enrichment). 
#These are very interesting predictions as key regulators of your intercellular communication process of interest !
nichenet_output$ligand_activity_target_heatmap

#3.Inferred ligand-receptor interactions for top-ranked ligands
nichenet_output$ligand_receptor_heatmap

#extract the ligand-receptor links and their interaction confidence scores in matrix or data frame format 
nichenet_output$ligand_receptor_matrix %>% .[1:10,1:6]

nichenet_output$ligand_receptor_df # weight column accords to number of data sources that document this interaction

#To get a list of the receptors of the 20 top-ranked ligands
nichenet_output$top_receptors

DotPlot(sce %>% subset(idents = "EC"), 
        features = nichenet_output$top_receptors %>% rev(), 
        split.by = "treatment") + 
  RotatedAxis()

VlnPlot(sce %>% subset(idents="EC") , 
        features = nichenet_output$top_receptors, 
        split.by = "treatment",    
        pt.size = 0, 
        combine = T,
        ncol = 8)

#show ‘bona fide’ ligand-receptor links that are described in the literature and not predicted based on protein-protein interactions
nichenet_output$ligand_receptor_heatmap_bonafide

nichenet_output$ligand_receptor_matrix_bonafide
nichenet_output$ligand_receptor_df_bonafide

#checking which geneset (and background set of genes) was used during the ligand activity analysis
nichenet_output$geneset_oi
nichenet_output$background_expressed_genes %>% length()




