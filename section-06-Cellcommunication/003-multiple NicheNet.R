#一Run multiple NicheNet analyses on different receiver cell populations

#一键清空
rm(list=ls())
options(stringsAsFactors = F)

#载入包
#devtools::install_github("saeyslab/nichenetr")
library(nichenetr)
library(Seurat) # please update to Seurat V4
library(tidyverse)

#Read in the expression data of interacting cells
sce <- readRDS("sce_tutorial.rds")
#Read in NicheNet’s ligand-target prior model(先验模型), ligand-receptor network and weighted integrated networks(加权整合网络)
ligand_target_matrix <-  readRDS(url("https://zenodo.org/record/3260758/files/ligand_target_matrix.rds"))
lr_network <-  readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds"))
weighted_networks <-  readRDS(url("https://zenodo.org/record/3260758/files/weighted_networks.rds"))

receiver_celltypes_oi <-  c("EC", "MC")
# receiver_celltypes_oi = seuratObj %>% Idents() %>% unique() # for all celltypes in the dataset: use only when this would make sense biologically
nichenet_output <-  receiver_celltypes_oi %>% lapply(nichenet_seuratobj_aggregate, 
                                                     seurat_obj = sce, 
                                                     condition_colname = "treatment", 
                                                     condition_oi = "Diabetes", 
                                                     condition_reference = "Control", 
                                                     sender = c("POD","IC", "TC"), 
                                                     ligand_target_matrix = ligand_target_matrix, 
                                                     lr_network = lr_network, 
                                                     weighted_networks = weighted_networks, 
                                                     organism = "mouse")

names(nichenet_output) <-  receiver_celltypes_oi

#Check which ligands were top-ranked for both EC and MC and which ligands were more cell-type specific
common_ligands = intersect(nichenet_output$`EC`$top_ligands, nichenet_output$`MC`$top_ligands)
print("common ligands are: ")
print(common_ligands)

EC_ligands = nichenet_output$`EC`$top_ligands %>% setdiff(nichenet_output$`MC`$top_ligands)
MC_ligands = nichenet_output$`MC`$top_ligands %>% setdiff(nichenet_output$`EC`$top_ligands)

print("Ligands specifically regulating DE in EC: ")
print(EC_ligands)

print("Ligands specifically regulating DE in MC: ")
print(MC_ligands)


#二explain differential expression between two cell populations
#first change the seuratObject of the data
sce@meta.data$celltype <-  paste(sce@meta.data$new_clusters,sce@meta.data$treatment, sep = "_")

sce@meta.data$celltype %>% table()

sce <-  SetIdent(sce,value = "celltype")
#Now perform the NicheNet analysis to explain differential expression between the ‘affected’ cell population ‘CD8 T cells after LCMV infection’ and the reference cell population ‘CD8 T cells in steady-state’ by ligands expressed by monocytes and DCs after LCMV infection.
nichenet_output <-  nichenet_seuratobj_cluster_de(
  seurat_obj = sce, 
  receiver_reference = "EC_Control", receiver_affected = "EC_Diabetes", 
  sender = c("MC_Diabetes","POD_Diabetes"), 
  ligand_target_matrix = ligand_target_matrix, 
  lr_network = lr_network, 
  weighted_networks = weighted_networks, 
  organism = "mouse")
#Check the top-ranked ligands and their target genes
nichenet_output$ligand_activity_target_heatmap



