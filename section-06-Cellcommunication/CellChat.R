#一键清空
rm(list=ls())
options(stringsAsFactors = F)

#载入包
#devtools::install_github("sqjin/CellChat")
library(CellChat)
library(patchwork)
library(Seurat)
library(dplyr)
library(ggplot2)
library(ggalluvial)
library(svglite)

####第一部分####
#1.数据输入及预处理
sce <- readRDS("sce_tutorial.rds")
sce@commands$FindClusters  # 你也看一看作者的其他命令，Seurat是记录其分析过程的

#2.创建CellChat 对象,从Seurat对象直接创建
cellchat <- createCellChat(object = sce,
                           group.by = "new_clusters")
cellchat
summary(cellchat)
str(cellchat)
#将细胞信息添加到对象的meta slot中
#cellchat <- addMeta(cellchat, meta = meta)##增加其他meta信息
#cellchat <- setIdent(cellchat, ident.use = "new_clusters") # 将 "labels" 设为默认细胞标记类型，这个可以根据自己的数据自定义
levels(cellchat@idents)
groupSize <- as.numeric(table(cellchat@idents)) # 每组细胞的数量
#查看每个cluster有多少细胞，之后画图会用
groupSize

#3.设置配体受体交互数据库
#CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
CellChatDB <- CellChatDB.mouse
# 查看数据库具体信息,包含interaction、cofactor、complex、geneInfo四个数据库
colnames(CellChatDB$interaction)
CellChatDB$interaction[1:4,1:4]
head(CellChatDB$interaction)
head(CellChatDB$cofactor)
head(CellChatDB$complex)
head(CellChatDB$geneInfo)

showDatabaseCategory(CellChatDB)
#人类中的CellChatDB包含1,939个经过验证的分子相互作用，
#包括61.8％的旁分泌/自分泌信号相互作用，
#21.7％的细胞外基质（ECM） - 受体受体相互作用和16.5％的细胞细胞接触相互作用
#在CellChat中，我们还可以先择特定的信息描述细胞间的相互作者，这个可以理解为从特定的侧面来刻画细胞间相互作用, 比用一个大的配体库又精细了许多
# 可以选择的子数据库，用于分析细胞间相互作用
unique(CellChatDB$interaction$annotation)
# [1] "Secreted Signaling" "ECM-Receptor"       "Cell-Cell Contact"

# Show the structure of the database
dplyr::glimpse(CellChatDB$interaction)

# use a subset of CellChatDB for cell-cell communication analysis
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
# use all CellChatDB for cell-cell communication analysis
# CellChatDB.use <- CellChatDB # simply use the default CellChatDB

# set the used database in the object
cellchat@DB <- CellChatDB.use

#4.预处理用于细胞通信分析的表达数据
# subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat) 
future::plan("multisession", workers = 4) # do parallel
#相当于Seurat里的FindMarkers,找每个细胞群中高表达的配体和受体
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
# #储存上一步的结果到cellchat@LR$LRsig
cellchat <- projectData(cellchat, PPI.human)

#####第二部分#####
#1.细胞通信网络的推断
#计算通信概率并推断cellchat网络
cellchat <- computeCommunProb(cellchat, 
                              raw.use = FALSE,
                              population.size =TRUE )#如果不想用上一步PPI矫正的结果，raw.use=TRUE
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)
#提取推断的cellchat网络作为数据框架
df.net <- subsetCommunication(cellchat)
df.net <- subsetCommunication(cellchat, sources.use = c(1,2), targets.use = c(4,5))#给出从小区组 1 和 2 到小区组 4 和 5 的推断小区间通信。
df.net <- subsetCommunication(cellchat, signaling = c("WNT", "TGFb"))#给出了由信号传导 WNT 和 TGFb 介导的细胞间通讯。
write.csv(df.net, "test_cellchat_lr.csv",quote = F)

#2.在信号通路水平推断细胞-细胞通信（结果储存在@netP下面，有一个概率值和对应的pval）
cellchat <- computeCommunProbPathway(cellchat)
#NB：每个配体受体对和每个信号通路的推断细胞间通信网络分别存储在插槽"net"和"netP"中
df.netp <- subsetCommunication(cellchat,slot.name = "netP")
write.csv(df.netp,"netP_pathway.csv")

#3.计算整合的细胞通信网络,统计细胞-细胞之间通信的数量(有多少个配体-受体对)和强度(概率)
cellchat <- aggregateNet(cellchat)
#计算每种细胞各有多少个
groupSize <- as.numeric(table(cellchat@idents))
groupSize
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, 
                 vertex.weight = groupSize, 
                 weight.scale = T, 
                 label.edge= F, 
                 title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, 
                 vertex.weight = groupSize, 
                 weight.scale = T, 
                 label.edge= F, 
                 title.name = "Interaction weights/strength")
par(opar)
#左图：外周各种颜色圆圈的大小代表细胞的数量，圆越大，数量越多，
#发出箭头的细胞表达配体，箭头指向的细胞表达受体，配-受体对越多，线越粗
#右图：互作的概率/强度值(强度就是概率值相加)

#4.检查每种细胞组发送的信号
mat <- cellchat@net$weight
par(mfrow = c(2,3), mar=c(1,1,1,1),xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
par(opar)
#如果报错，设置以下参数重新运行上述代码
#par("mar")
#par(mar=c(1,1,1,1))

#####第三部分：细胞通信网络的可视化######
#（层次图、网络图、和弦图、热图）
#层次结构图：用户应定义vertex.receiver，这是一个数字矢量，将细胞群的索引作为层次图左侧的目标。此分层图由两个部分组成：左部分显示自分泌和旁分泌向某些感兴趣的细胞组（即定义的）发出信号，右部分显示自分泌和旁分泌向数据集中剩余的细胞组发出信号。因此，层级图提供了一种信息性和直观的方式来可视化自分泌和旁分泌信号之间的细胞群之间的感兴趣通信。例如，在研究成纤维细胞和免疫细胞之间的细胞-细胞通信时，用户可以定义为所有成纤维细胞组。
#和弦图：CellChat 提供两种功能netVisual_chord_cell和netVisual_chord_gene，并可视化具有不同目的和不同级别的细胞通信。netVisual_chord_cell用于可视化不同细胞群之间的细胞-细胞通信（和弦图中的每个部分是细胞组），netVisual_chord_gene用于可视化由多个配体受体或信号通路调节的细胞-细胞通信（和弦图中的每个部分都是配体、受体或信号通路）。
#边缘颜色/权重、节点颜色/大小/形状的解释：在所有可视化图中，边缘颜色与发送者源一致，边缘权重与交互强度成正比。较厚的边缘线表示信号更强。在层次结构图和圆图中，圆的大小与每个细胞组中的细胞数量成正比。在层次图中，实心和开放的圆分别代表源和目标。在和弦图中，内条颜色表示从相应的外条接收信号的目标。内条大小与目标接收的信号强度成正比。这种内条有助于解释复杂的和弦图。请注意，有一些内条没有与任何一些细胞组链接，请忽略它，因为这是一个本包尚未解决的问题。
#不同层次的细胞通信可视化：可以使用netVisual_aggregate可视化信号通路的推断通信网络，并使用netVisual_individual可视化与该信号通路相关的单个L-R对的推断通信网络。
#所有显示重要通信的信号通路均可通过cellchat@netP$pathways获取。
cellchat@netP$pathways
#在这里，我们以输入一个信号通路为例。
pathways.show <- c("VEGF") 

# 1.Hierarchy plot层次图
levels(cellchat@idents)
# Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 
vertex.receiver = seq(1,4) # a numeric vector.细胞群索引 
netVisual_aggregate(cellchat, 
                    signaling = pathways.show,  
                    vertex.receiver = vertex.receiver)
#在层次图中，实心圆和空心圆分别表示源和目标，圆的大小和细胞数量成比例，线越粗，互作信号越强。
#左图中target是我们选中的细胞群，右图是选中的细胞群之外的细胞

# 2.Circle plot网络图
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, 
                    signaling = pathways.show, 
                    layout = "circle")

#3. Chord diagram和弦图
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, 
                    signaling = pathways.show, 
                    layout = "chord")

# 4.Heatmap热图
par(mfrow=c(1,1))
netVisual_heatmap(cellchat, 
                  signaling = pathways.show, 
                  color.heatmap = "Reds")
#纵轴是发出信号的细胞，横轴是接受信号的细胞，热图颜色深浅代表信号强弱，上侧和右侧是纵轴和横轴强度的积累

#对于和弦图，CellChat 具有独立函数netVisual_chord_cell，通过调整circlize包中的不同参数来灵活可视化信号网络。例如，我们可以定义一个group命名的字符矢量，以创建多组和弦图，将细胞群集分组到不同的细胞类型。
# Chord diagram
group.cellType <- c(rep("EC", 4), rep("MC", 4), rep("TC", 4)) # grouping cell clusters into fibroblast, DC and TC cells
names(group.cellType) <- levels(cellchat@idents)
netVisual_chord_cell(cellchat, signaling = pathways.show, group = group.cellType, title.name = paste0(pathways.show, " signaling network"))

#5.计算每个配体受体对整体信号通路的贡献，并可视化由单个配体受体对调节的细胞通信
netAnalysis_contribution(cellchat, 
                         signaling = pathways.show)
#可视化由单个配体受体对调节的细胞-细胞通信。我们提供一个函数extractEnrichedLR来提取给定信号通路的所有重要相互作用（L-R对）和相关信号基因
#提取对VEGF有贡献的所有配体受体
pairLR.VEGF <- extractEnrichedLR(cellchat, 
                                 signaling = pathways.show, 
                                 geneLR.return = FALSE)
LR.show <- pairLR.VEGF[1,] # show one ligand-receptor pair
# Hierarchy plot
vertex.receiver = c(1,4) # a numeric vector
netVisual_individual(cellchat, 
                     signaling = pathways.show,  
                     pairLR.use = LR.show, 
                     vertex.receiver = vertex.receiver)
# Circle plot
netVisual_individual(cellchat, 
                     signaling = pathways.show, 
                     pairLR.use = LR.show, 
                     layout = "circle")

# Chord diagram
netVisual_individual(cellchat, 
                     signaling = pathways.show, 
                     pairLR.use = LR.show, 
                     layout = "chord")

# 6.自动（批量）保存每个信号通路的互作结果
# Access all the signaling pathways showing significant communications
pathways.show.all <- cellchat@netP$pathways
# check the order of cell identity to set suitable vertex.receiver
levels(cellchat@idents)
vertex.receiver = seq(1,4)
for (i in 1:length(pathways.show.all)) {
  # Visualize communication network associated with both signaling pathway and individual L-R pairs
  netVisual(cellchat, signaling = pathways.show.all[i], vertex.receiver = vertex.receiver, layout = "hierarchy")
  # Compute and visualize the contribution of each ligand-receptor pair to the overall signaling pathway
  gg <- netAnalysis_contribution(cellchat, signaling = pathways.show.all[i])
  ggsave(filename=paste0(pathways.show.all[i], "_L-R_contribution.pdf"), plot=gg, width = 3, height = 2, units = 'in', dpi = 300)
}

#7.可视化由多个配体受体或信号通路调节的细胞通信
# show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
netVisual_bubble(cellchat, sources.use = c(4,5), 
                 targets.use = c(1,2,3), 
                 remove.isolate = FALSE)
#> Comparing communications on a single object

# show all the significant interactions (L-R pairs) associated with certain signaling pathways
#指定信号通路或配体-受体
netVisual_bubble(cellchat, sources.use = c(4,5), 
                 targets.use = c(1,2,3), 
                 signaling = c("VEGF","CXCL"),
                 remove.isolate = FALSE)
#指定信号通路或配体-受体并指定细胞
pairLR.use <- extractEnrichedLR(cellchat,
                                signaling = c("VEGF","CXCL"))
netVisual_bubble(cellchat, sources.use = c(4,5), 
                 targets.use = c(1),
                 pairLR.use = pairLR.use,
                 remove.isolate = FALSE)
#VEGF,CXCL信号通路参与IC和TC对EC的调控作用情况

# show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
# show all the interactions sending from Inflam.FIB
netVisual_chord_gene(cellchat, 
                     sources.use = 4, 
                     targets.use = c(1,3,5), 
                     lab.cex = 0.5,
                     legend.pos.y = 30)

# show all the interactions received by Inflam.DC
netVisual_chord_gene(cellchat, 
                     sources.use = c(1,2,3,4), 
                     targets.use = 5, 
                     legend.pos.x = 15)

# show all the significant interactions (L-R pairs) associated with certain signaling pathways
netVisual_chord_gene(cellchat, 
                     sources.use = c(2), 
                     targets.use = c(3,4,5),
                     signaling = c("VEGF","CXCL"),
                     legend.pos.x = 8)

# show all the significant signaling pathways from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
netVisual_chord_gene(cellchat, 
                     sources.use = c(1,2,3,4), 
                     targets.use = c(5), 
                     slot.name = "netP", 
                     legend.pos.x = 10)

#8.使用小提琴/点图绘制参与某条信号通路的(如VEGF)的所有基因在细胞群中的表达情况
plotGeneExpression(cellchat, signaling = "VEGF")

#默认情况下，用户可以通过plotGeneExpression只显示与推断的重要通信相关的信号基因的表达。
plotGeneExpression(cellchat, signaling = "VEGF", 
                   enriched.only = FALSE)

plotGeneExpression(cellchat, signaling = "VEGF",
                   type = "dot" ,
                   color.use = c("red"))

####第四部分：细胞通信网络系统分析####
#识别细胞组的信号角色（例如，占主导地位的发送器、接收器）以及主要贡献信号
#1.计算和可视化网络中心分数
# Compute the network centrality scores
cellchat <- netAnalysis_computeCentrality(cellchat, 
                                          slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
netAnalysis_signalingRole_network(cellchat, 
                                  signaling = pathways.show, 
                                  width = 8, 
                                  height = 2.5, 
                                  font.size = 10)

#2.在 2D 空间中可视化占主导地位的发送器（源）和接收器（目标）
# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
gg1 <- netAnalysis_signalingRole_scatter(cellchat)
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
# Signaling role analysis on the cell-cell communication networks of interest
gg2 <- netAnalysis_signalingRole_scatter(cellchat, signaling = c("CXCL", "VEGF"))
#> Signaling role analysis on the cell-cell communication network from user's input
gg1 + gg2

#3.识别对某些细胞组的传出或传入信号贡献最大的信号
# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")
ht1 + ht2

# Signaling role analysis on the cell-cell communication networks of interest
ht <- netAnalysis_signalingRole_heatmap(cellchat, signaling = c("CXCL", "VEGF"))

#4.识别和可视化分泌细胞的传出通信模式
library(NMF)
library(ggalluvial)
#运行selectK推断模式的数量
selectK(cellchat, pattern = "outgoing")
#当传出模式数为 5 时，Cophenetic 和Silhouette值都开始突然下降。
nPatterns = 5
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns)

# river plot
netAnalysis_river(cellchat, pattern = "outgoing")

# dot plot
netAnalysis_dot(cellchat, pattern = "outgoing")

#5.识别和可视化目标细胞的传入通信模式
selectK(cellchat, pattern = "incoming")
#当传入模式的数量为 4 时，Cophenetic 值开始下降
nPatterns = 4
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = nPatterns)
# river plot
netAnalysis_river(cellchat, pattern = "incoming")

# dot plot
netAnalysis_dot(cellchat, pattern = "incoming")

#6.信号网络的多重和分类学习分析
#根据信号组的功能相似性识别信号组
cellchat <- computeNetSimilarity(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat, type = "functional")

cellchat <- netClustering(cellchat, type = "functional")

# Visualization in 2D-space
netVisual_embedding(cellchat, type = "functional", label.size = 3.5)

# netVisual_embeddingZoomIn(cellchat, type = "functional", nCol = 2)
#基于结构相似性识别信号组
cellchat <- computeNetSimilarity(cellchat, type = "structural")
cellchat <- netEmbedding(cellchat, type = "structural")

cellchat <- netClustering(cellchat, type = "structural")
# Visualization in 2D-space
netVisual_embedding(cellchat, type = "structural", label.size = 3.5)

netVisual_embeddingZoomIn(cellchat, type = "structural", nCol = 2)

####第五部分：保存cellchat对象####
saveRDS(cellchat, file = "cellchat_mouse.rds")






