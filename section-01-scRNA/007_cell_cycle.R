##周期判断----
#在挑选hvg gene那一步因，可能会找到一些细胞周期相关基因；
#它们会导致细胞聚类发生一定的偏移，即相同类型的细胞在聚类时会因为细胞周期的不同而分开。
#因此有必要查看是否有细胞周期相关基因的存在；若有，则剔除

#细胞周期有关基因
?cc.genes
head(c(cc.genes$s.genes,cc.genes$g2m.genes))
length(c(cc.genes$s.genes,cc.genes$g2m.genes))
# [1] "MCM5" "PCNA" "TYMS" "FEN1" "MCM2" "MCM4"
#查看我们选择的高变基因中有哪些细胞周期相关基因,及打分
CaseMatch(c(cc.genes$s.genes,cc.genes$g2m.genes),VariableFeatures(sce))
#在scRNA@meta.data中添加S.Score、G2M.Score和Phase三列有关细胞周期的信息。
g2m_genes = cc.genes$g2m.genes
g2m_genes = CaseMatch(search = g2m_genes, match = rownames(sce))
s_genes = cc.genes$s.genes
s_genes = CaseMatch(search = s_genes, match = rownames(sce))
scRNA <- CellCycleScoring(object=sce,  g2m.features=g2m_genes,  s.features=s_genes)
head(sce@meta.data)
#观察细胞周期相关基因是否影响聚类
scRNA <- RunPCA(sce, features = c(s_genes, g2m_genes))
p1 <- DimPlot(sce, reduction = "pca", group.by = "Phase")
ggsave("../cell-cycle.pdf", plot = p1) 
#影响不大，基本重合在一起了