rm(list = ls())
options(stringsAsFactors = F)
#加载包
library(Seurat)
library(ggplot2)
library(ggsci)
library(ggrepel)
library(tidyverse)

#加载数据
sce <- readRDS("sce_tutorial.rds")
#差异分析:这里 logfc.threshold 设置为 0 
sce.markers <- FindAllMarkers(sce, only.pos = FALSE,#不只返回正向变化基因
                               min.pct = 0.25,
                               logfc.threshold = 0)
#筛选前五和后 5 的基因作为展示
# top 5 genes
top5pos <- sce.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top5negtive <- sce.markers %>% group_by(cluster) %>% top_n(n = -5, wt = avg_log2FC)
# merge
top10 <- rbind(top5pos,top5negtive)

ggplot(sce.markers,
       aes(x = pct.2 - pct.1,y = avg_log2FC)) +
  geom_point(color = 'grey80') +
  # 添加水平线
  geom_hline(yintercept = c(-0.25,0.25),lty = 'dashed',size = 1,color = 'grey50') +
  # 添加 top5pos 基因标签
  geom_text_repel(data = top5pos,
                  aes(x = pct.2 - pct.1,y = avg_log2FC,
                      label = gene,color = cluster),
                  show.legend = F,direction = 'y',
                  hjust = 1, # 右对齐
                  nudge_y = 0.25, # 文字竖直方向调整
                  force = 5, # 文字重叠调整
                  # 文字靠右竖直对齐
                  nudge_x = 0.8 - (top5pos$pct.2 - top5pos$pct.1)) +
  # 添加 top5negtive 基因标签
  geom_text_repel(data = top5negtive,
                  aes(x = pct.2 - pct.1,y = avg_log2FC,
                      label = gene,color = cluster),
                  show.legend = F,direction = 'y',
                  hjust = 0, # 左对齐
                  force = 2.5, # 文字重叠调整
                  # 文字靠左竖直对齐
                  nudge_x = -0.8 - (top5negtive$pct.2 - top5negtive$pct.1)) +
  # top10 点颜色
  geom_point(data = top10,show.legend = F,
             aes(x = pct.2 - pct.1,y = avg_log2FC,color = cluster)) +
  scale_color_npg(name = '') +
  # x y breaks label
  scale_y_continuous(limits = c(-6,10),breaks = seq(-6,10,2)) +
  scale_x_continuous(limits = c(-1,1),breaks = seq(-1,1,0.5)) +
  theme_bw(base_size = 14) +
  # 主题调整
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45,hjust = 1),
        strip.background = element_rect(color = NA,fill = 'grey90')) +
  # 轴标签
  xlab(expression(Delta~'Percentage Diffrence')) +
  ylab('Log2-Fold Change') +
  # 分面
  facet_wrap(~cluster,nrow = 1,scales = 'fixed')

ggsave("marker_genes.png",width = 2500,height = 1500,
       dpi = 300,units = "px")
