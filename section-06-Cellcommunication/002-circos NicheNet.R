#接001-single NicheNet

#Circos绘图来可视化配体-靶标和配体-受体的相互作用

#这一可视化分组根据最强表达的细胞类型预测活性配体。
#因此，我们需要确定每种细胞类型，它们表达的配体比其他细胞类型更强。
#计算发送细胞中平均配体表达量。

# avg_expression_ligands = AverageExpression(seuratObj %>% subset(subset = aggregate == "LCMV"),features = nichenet_output$top_ligands) # if want to look specifically in LCMV-only cells
avg_expression_ligands <-  AverageExpression(sce, features = nichenet_output$top_ligands)

#分配配体给发送细胞
#为了给发送端细胞类型分配配体，我们可以查找哪个发送端细胞类型的表达式高于平均值+ SD。
sender_ligand_assignment <-  avg_expression_ligands$RNA %>% apply(1, function(ligand_expression){
  ligand_expression > (ligand_expression %>% mean() + ligand_expression %>% sd())
}) %>% t()

sender_ligand_assignment[1:4,1:4]


sender_ligand_assignment <-  sender_ligand_assignment %>% apply(2, function(x){x[x == TRUE]}) %>% purrr::keep(function(x){length(x) > 0})
names(sender_ligand_assignment)

(sender_ligand_assignment)

#现在确定哪些优先配体是由CAFs或内皮细胞表达的
all_assigned_ligands <-  sender_ligand_assignment %>% lapply(function(x){names(x)}) %>% unlist()
unique_ligands <-  all_assigned_ligands %>% table() %>% .[. == 1] %>% names()
general_ligands <-  nichenet_output$top_ligands %>% setdiff(unique_ligands)

EC_specific_ligands <-  sender_ligand_assignment$EC %>% names() %>% setdiff(general_ligands)
MC_specific_ligands <-  sender_ligand_assignment$MC %>% names() %>% setdiff(general_ligands)
POD_specific_ligands <-  sender_ligand_assignment$POD %>% names() %>% setdiff(general_ligands)
IC_specific_ligands <-  sender_ligand_assignment$IC %>% names() %>% setdiff(general_ligands)
TC_specific_ligands <-  sender_ligand_assignment$TC %>% names() %>% setdiff(general_ligands)


ligand_type_indication_df <-  tibble(
  ligand_type = c(rep("EC-specific", times = EC_specific_ligands %>% length()),
                  rep("MC-specific", times = MC_specific_ligands %>% length()),
                  rep("POD-specific", times = POD_specific_ligands %>% length()),
                  rep("IC-specific", times = IC_specific_ligands %>% length()),
                  rep("TC-specific", times = TC_specific_ligands %>% length()),
                  rep("General", times = general_ligands %>% length())),
  ligand = c(EC_specific_ligands, MC_specific_ligands, POD_specific_ligands, IC_specific_ligands, TC_specific_ligands,general_ligands))

ligand_type_indication_df %>% head

#定义感兴趣的配体-目标链接
#为了避免circos图中有太多配体目标链接，我们将只显示权重高于预定义截止值的链接:属于最低分数的40%的链接被删除。
#这并不是说用于这种可视化的边界和其他边界可以根据用户的需要进行更改。
active_ligand_target_links_df <-  nichenet_output$ligand_target_df %>% mutate(target_type = "Diabetes-DE") %>% inner_join(ligand_type_indication_df) # if you want ot make circos plots for multiple gene sets, combine the different data frames and differentiate which target belongs to which gene set via the target type

cutoff_include_all_ligands <-  active_ligand_target_links_df$weight %>% quantile(0.40)

active_ligand_target_links_df_circos <-  active_ligand_target_links_df %>% filter(weight > cutoff_include_all_ligands)

ligands_to_remove <-  setdiff(active_ligand_target_links_df$ligand %>% unique(), active_ligand_target_links_df_circos$ligand %>% unique())
targets_to_remove <-  setdiff(active_ligand_target_links_df$target %>% unique(), active_ligand_target_links_df_circos$target %>% unique())

circos_links <-  active_ligand_target_links_df %>% filter(!target %in% targets_to_remove &!ligand %in% ligands_to_remove)

circos_links
#准备circos可视化:给每个片段配体和目标特定的颜色和顺序
colpalette <- c("#F8766D" ,"#A3A500", "#00BF7D" ,"#00B0F6", "#E76BF3")

grid_col_ligand <- c("General" = "lawngreen",
                   "MC-specific" = "#F8766D",
                   "EC-specific" = "#A3A500",
                   "POD-specific" = "#00BF7D",
                   "IC-specific" = "#00B0F6",
                   "TC-specific" = "#E76BF3")
grid_col_target <- c(
  "Diabetes-DE" = "tomato")

grid_col_tbl_ligand <- tibble(ligand_type = grid_col_ligand %>% names(), color_ligand_type = grid_col_ligand)
grid_col_tbl_target <- tibble(target_type = grid_col_target %>% names(), color_target_type = grid_col_target)

circos_links <- circos_links %>% mutate(ligand = paste(ligand," ")) # extra space: make a difference between a gene as ligand and a gene as target!
circos_links <- circos_links %>% inner_join(grid_col_tbl_ligand) %>% inner_join(grid_col_tbl_target)
links_circle <- circos_links %>% select(ligand,target, weight)

ligand_color <- circos_links %>% distinct(ligand,color_ligand_type)
grid_ligand_color <- ligand_color$color_ligand_type %>% set_names(ligand_color$ligand)
target_color <- circos_links %>% distinct(target,color_target_type)
grid_target_color <- target_color$color_target_type %>% set_names(target_color$target)

grid_col <- c(grid_ligand_color,grid_target_color)

# give the option that links in the circos plot will be transparant ~ ligand-target potential score
transparency <- circos_links %>% mutate(weight =(weight-min(weight))/(max(weight)-min(weight))) %>% mutate(transparency = 1-weight) %>% .$transparency 

#准备可视化的circos:排序配体和目标
target_order <- circos_links$target %>% unique()
ligand_order <- c(TC_specific_ligands,POD_specific_ligands, IC_specific_ligands, MC_specific_ligands,EC_specific_ligands, general_ligands) %>% c(paste(.," ")) %>% intersect(circos_links$ligand)
order <- c(ligand_order,target_order)

#准备circos可视化:定义不同片段之间的间隙
width_same_cell_same_ligand_type <- 0.5
width_different_cell <- 6
width_ligand_target <- 15
width_ligand_receptor <- 15
width_same_cell_same_target_type <- 0.5
width_same_cell_same_receptor_type <- 0.5

#绘图1.
gaps <- c(
  # width_ligand_target,
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "TC-specific") %>% distinct(ligand) %>% nrow() -1)),
  width_different_cell,
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "POD-specific") %>% distinct(ligand) %>% nrow() -1)),
  width_different_cell,
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "IC-specific") %>% distinct(ligand) %>% nrow() -1)),
  width_different_cell,
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "MC-specific") %>% distinct(ligand) %>% nrow() -1)),
  width_different_cell,
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "EC-specific") %>% distinct(ligand) %>% nrow() -1)),
  width_different_cell,
  #rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "General") %>% distinct(ligand) %>% nrow() )),
  #width_ligand_target,
  rep(width_same_cell_same_target_type, times = (circos_links %>% filter(target_type == "Diabetes-DE") %>% distinct(target) %>% nrow() -1)),
  width_ligand_target
)
#渲染circos的情节(所有链接相同的透明度)。只有表明每个靶基因的阻滞的宽度与配体-靶的调控电位成正比(~支持调控相互作用的先验知识)
circos.par(gap.degree = gaps)
chordDiagram(links_circle, 
             directional = 1,
             order=order,
             link.sort = TRUE, 
             link.decreasing = FALSE, 
             grid.col = grid_col,
             transparency = 0, 
             diffHeight = 0.005, 
             direction.type = c("diffHeight", "arrows"),
             link.arr.type = "big.arrow", 
             link.visible = links_circle$weight >= cutoff_include_all_ligands,
             annotationTrack = "grid", 
             preAllocateTracks = list(track.height = 0.075))
# we go back to the first track and customize sector labels
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.55), cex = 1)
}, bg.border = NA) #

circos.clear()

#绘制circos图(透明度由配体-靶标相互作用的调控潜力决定)
circos.par(gap.degree = gaps)
chordDiagram(links_circle, directional = 1,order=order,link.sort = TRUE, link.decreasing = FALSE, grid.col = grid_col,transparency = transparency, diffHeight = 0.005, direction.type = c("diffHeight", "arrows"),link.arr.type = "big.arrow", link.visible = links_circle$weight >= cutoff_include_all_ligands,annotationTrack = "grid", 
             preAllocateTracks = list(track.height = 0.075))
# we go back to the first track and customize sector labels
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.55), cex = 1)
}, bg.border = NA) #

circos.clear()

svg("ligand_target_circos.svg", width = 10, height = 10)
circos.par(gap.degree = gaps)
chordDiagram(links_circle, directional = 1,order=order,link.sort = TRUE, link.decreasing = FALSE, grid.col = grid_col,transparency = transparency, diffHeight = 0.005, direction.type = c("diffHeight", "arrows"),link.arr.type = "big.arrow", link.visible = links_circle$weight >= cutoff_include_all_ligands,annotationTrack = "grid",
             preAllocateTracks = list(track.height = 0.075))
# we go back to the first track and customize sector labels
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.55), cex = 1)
}, bg.border = NA) #
circos.clear()
dev.off()

#绘图2.在circos图中可视化优先配体与受体的相互作用
lr_network_top_df <-  nichenet_output$ligand_receptor_df %>% mutate(receptor_type = "Diabetes_MC_receptor") %>% inner_join(ligand_type_indication_df)
grid_col_ligand <-  c("General" = "lawngreen",
                   "MC-specific" = "#F8766D",
                   "EC-specific" = "#A3A500",
                   "POD-specific" = "#00BF7D",
                   "IC-specific" = "#00B0F6",
                   "TC-specific" = "#E76BF3")
grid_col_receptor <-  c(
  "Diabetes_MC_receptor" = "darkred")

grid_col_tbl_ligand <-  tibble(ligand_type = grid_col_ligand %>% names(), color_ligand_type = grid_col_ligand)
grid_col_tbl_receptor <-  tibble(receptor_type = grid_col_receptor %>% names(), color_receptor_type = grid_col_receptor)

circos_links <-  lr_network_top_df %>% mutate(ligand = paste(ligand," ")) # extra space: make a difference between a gene as ligand and a gene as receptor!
circos_links <-  circos_links %>% inner_join(grid_col_tbl_ligand) %>% inner_join(grid_col_tbl_receptor)
links_circle <-  circos_links %>% select(ligand,receptor, weight)

ligand_color <-  circos_links %>% distinct(ligand,color_ligand_type)
grid_ligand_color <-  ligand_color$color_ligand_type %>% set_names(ligand_color$ligand)
receptor_color <-  circos_links %>% distinct(receptor,color_receptor_type)
grid_receptor_color <-  receptor_color$color_receptor_type %>% set_names(receptor_color$receptor)

grid_col <-  c(grid_ligand_color,grid_receptor_color)

# give the option that links in the circos plot will be transparant ~ ligand-receptor potential score
transparency <-  circos_links %>% mutate(weight =(weight-min(weight))/(max(weight)-min(weight))) %>% mutate(transparency = 1-weight) %>% .$transparency 

#制备可视化的circos:有序配体和受体
receptor_order <-  circos_links$receptor %>% unique()
ligand_order <-  c(TC_specific_ligands,POD_specific_ligands, IC_specific_ligands, MC_specific_ligands,EC_specific_ligands, general_ligands) %>% c(paste(.," ")) %>% intersect(circos_links$ligand)
order <-  c(ligand_order,receptor_order)

#准备circos可视化:定义不同片段之间的间隙
gaps <- c(
  # width_ligand_target,
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "TC-specific") %>% distinct(ligand) %>% nrow() -1)),
  width_different_cell,
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "POD-specific") %>% distinct(ligand) %>% nrow() -1)),
  width_different_cell,
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "IC-specific") %>% distinct(ligand) %>% nrow() -1)),
  width_different_cell,
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "MC-specific") %>% distinct(ligand) %>% nrow() -1)),
  width_different_cell,
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "EC-specific") %>% distinct(ligand) %>% nrow() -1)),
  width_different_cell,
  #rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "General") %>% distinct(ligand) %>% nrow()-1 )),
  #width_ligand_target,
  rep(width_same_cell_same_receptor_type, times = (circos_links %>% filter(receptor_type == "Diabetes_MC_receptor") %>% distinct(receptor) %>% nrow() -1)),
  width_ligand_receptor
)
#渲染circos的情节(所有链接相同的透明度)。只有表明每个受体的阻滞的宽度与配体-受体相互作用的重量成比例(~支持相互作用的先验知识)
circos.par(gap.degree = gaps)
chordDiagram(links_circle, directional = 1,
             order=order,
             link.sort = TRUE, 
             link.decreasing = FALSE, 
             grid.col = grid_col,
             transparency = 0, 
             diffHeight = 0.005, 
             direction.type = c("diffHeight", "arrows"),
             link.arr.type = "big.arrow", 
             link.visible = links_circle$weight >= cutoff_include_all_ligands,
             annotationTrack = "grid", 
             preAllocateTracks = list(track.height = 0.075))
# we go back to the first track and customize sector labels
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.55), cex = 0.8)
}, bg.border = NA) #

circos.clear()

#渲染circos图(透明程度由配体-受体相互作用的先验相互作用权重决定——正如指示每个受体的块的宽度)
circos.par(gap.degree = gaps)
chordDiagram(links_circle, directional = 1,
             order=order,
             link.sort = TRUE, 
             link.decreasing = FALSE, 
             grid.col = grid_col,
             transparency = transparency, 
             diffHeight = 0.005, 
             direction.type = c("diffHeight", "arrows"),
             link.arr.type = "big.arrow", 
             link.visible = links_circle$weight >= cutoff_include_all_ligands,
             annotationTrack = "grid", 
             preAllocateTracks = list(track.height = 0.075))
# we go back to the first track and customize sector labels
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.55), cex = 0.8)
}, bg.border = NA) 

circos.clear()



