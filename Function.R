########## DEG-GO circular visualization #######
library(clusterProfiler)
library(org.Hs.eg.db)  # 人类基因注释数据库
library(ggplot2)
library(dplyr)
library(circlize)

sig_genes <- results %>%
  filter(adj.P.Val < 0.05 & abs(logFC) > 1) %>%
  rownames_to_column("Gene") %>%
  pull(Gene)
gene_symbols <- sig_genes
entrez_ids <- mapIds(org.Hs.eg.db,
                     keys = gene_symbols,
                     column = "ENTREZID",
                     keytype = "SYMBOL",
                     multiVals = "first")
entrez_ids <- entrez_ids[!is.na(entrez_ids)]
# 然后用Entrez ID进行富集分析
go_enrich <- enrichGO(
  gene = entrez_ids,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "ALL",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  readable = TRUE  # 结果中显示Gene Symbol
)
ego1_df <- data.frame(go_enrich) %>%
  group_by(ONTOLOGY) %>%
  arrange(pvalue) %>%
  slice_head(n = 6) %>%
  rowwise() %>%
  mutate(fc = eval(parse(text = GeneRatio))/eval(parse(text = BgRatio)))
#saveRDS(ego1_df,file = 'go_df.rds')
#ego1_df=readRDS('go_df.rds')


library(org.Hs.eg.db)  # 人类基因注释数据库
library(ggplot2)
library(dplyr)
library(circlize)
setwd('/Users/wuhongbo/Desktop/wjx/data_quality/output/JTM_1_cross')
ego1_df=readRDS('go_df.rds')

cairo_pdf("GO_Circos_Enrichment.pdf", width = 4.5, height = 4.5, onefile = TRUE)

layout(matrix(c(1, 2), nrow = 2, ncol = 1), heights = c(0.85, 0.15))
par(mar = c(1, 1, 1, 1))
my_colors <- c("#B6D0E2","#E7C7DC", "#a7dace") # 粉色、浅黄、青绿色
names(my_colors) <- c("BP", "MF", "CC") # 假设这是你的ONTOLOGY类别
col <- my_colors[as.character(ego1_df$ONTOLOGY)]

circos.clear()
circos.initialize(sectors = ego1_df$ID, xlim = c(0, 1))

par(cex = 0.8)
total_tracks <- 5 # 总共有5个轨道

# first GO id - 减小字体并减小高度
circos.track(sectors = ego1_df$ID, ylim = c(0, 1),
             bg.col = col,
             track.height = 0.12, # 调整高度
             panel.fun = function(x, y) {
               circos.text(x = CELL_META$xcenter, y = CELL_META$ycenter,
                           labels = CELL_META$sector.index,
                           cex = 0.6)
             })

# add count track
circos.track(sectors = ego1_df$ID, ylim = c(0, 1), track.height = 0.12,
             panel.fun = function(x, y) {
               # circos.axis(h = "bottom")
             })

for (i in 1:nrow(ego1_df)) {
  # 去除条形边框
  circos.rect(xleft = 0, xright = ego1_df$Count[i]/max(ego1_df$Count),
              ybottom = 0.25, ytop = 0.75,
              sector.index = ego1_df$ID[i],
              col = "#a7b2c7",
              border = NA)
  
  # add xaxis
  circos.axis(h = "bottom", major.at = c(0, 1),
              labels = c(0, max(ego1_df$Count)),
              sector.index = ego1_df$ID[i],
              labels.cex = 0.6)
}

# foldchange enrichment - 调整高度
circos.track(sectors = ego1_df$ID, ylim = c(0, 1), track.height = 0.08)
for (i in 1:nrow(ego1_df)) {
  circos.text(x = 0.5, y = 0.5,
              labels = round(ego1_df$fc[i], digits = 1), 
              sector.index = ego1_df$ID[i],
              cex = 0.6)
}


# 条形图 - 调整高度
circos.track(sectors = ego1_df$ID, ylim = c(0, ceiling(max(ego1_df$fc))),
             track.height = 0.2)
for (i in 1:nrow(ego1_df)) {
  circos.barplot(value = ego1_df$fc[i], pos = 0.5,
                 sector.index = ego1_df$ID[i],
                 col = "#FED7C3",
                 border = NA)
}

# -log10 pvalue 条形图 - 调整高度
circos.track(sectors = ego1_df$ID, ylim = c(0, ceiling(max(-log10(ego1_df$pvalue)))),
             track.height = 0.2)
for (i in 1:nrow(ego1_df)) {
  # 获取当前pvalue的-log10值
  pval_log10 <- -log10(ego1_df$pvalue[i])
  
  # 绘制条形图，去除边框
  circos.barplot(value = pval_log10, pos = 0.5,
                 sector.index = ego1_df$ID[i],
                 col = col[i],
                 border = NA)
  
  # 在条形图内部用黑色显示数值
  text_y <- pval_log10 * 0.8
  
  # 只在条形足够高时显示文本
  if (pval_log10 > max(-log10(ego1_df$pvalue)) * 0.2) {
    circos.text(x = 0.5, y = text_y,
                labels = format(round(pval_log10, 2), nsmall = 2),
                sector.index = ego1_df$ID[i],
                col = "black",
                cex = 0.5,
                facing = "inside",
                adj = c(0.5, 0.5))
  }
}

# 恢复默认字体大小
par(cex = 1)

# 重置底部的边距
par(mar = c(0, 0, 0, 0))
plot.new()

legend(
  x = 0.0, y = 1.0,
  legend = c("GO Categories:", names(my_colors)),
  fill = c(NA, my_colors),
  border = NA,
  text.font = c(2, 1, 1, 1),
  horiz = TRUE,
  bty = "n",
  cex = 0.75, # 🌟 字体调小
  x.intersp = 0.5 # 🌟 强行压缩色块与文字的间距
)

legend(
  x = 0.0, y = 0.65, # y=0.65 稍微靠近上排一点，显得紧凑
  legend = c("Track Metrics:", "Gene Count", "Fold Change"),
  fill = c(NA, "#a7b2c7", "#FED7C3"),
  border = NA,
  text.font = c(2, 1, 1),
  horiz = TRUE,
  bty = "n",
  cex = 0.75, # 🌟 字体调小
  x.intersp = 0.5 # 🌟 强行压缩间距
)

dev.off()




# 表格
library(grid)
library(gridExtra)
library(gtable)

# 你的数据
g_data <- data.frame(
  Ontology = ego1_df$ONTOLOGY,
  ID = ego1_df$ID,
  Description = ego1_df$Description,
  stringsAsFactors = FALSE
)

# 创建一组美观的颜色（10种）
pathway_colors <- rep(c("#B6D0E2","#E7C7DC", "#a7dace"), each=6)

# 更简单的版本，不使用gtable_add_grob
create_simple_table <- function(data, colors) {
  # 设置主题，在主题中直接设置行颜色
  my_theme <- ttheme_default(
    core = list(
      fg_params = list(
        hjust = 0,
        x = 0.1,
        fontsize = 11
      ),
      bg_params = list(
        # 为每一行设置不同的颜色
        fill = colors
      )
    ),
    colhead = list(
      fg_params = list(
        fontsize = 12,
        fontface = "bold",
        col = "white"
      ),
      bg_params = list(fill = "#2c3e50")
    )
  )
  # 创建表格
  table_grob <- tableGrob(
    d = data,
    theme = my_theme,
    rows = NULL
  )
  # 调整列宽
  table_grob$widths <- unit(c(0.11, 0.15,0.74), "npc")
  return(table_grob)
}

# 创建简单表格
table_g_plot <- create_simple_table(g_data, pathway_colors)

cairo_pdf("GO_Pathway_Table.pdf", width = 7, height = 7, onefile = TRUE)
grid.newpage()
# 先绘制表格
grid.draw(table_g_plot)
# 添加标题 (y坐标稍微往下一点，防止被裁掉)
grid.text("GO Pathway Enrichment Analysis",
          x = 0.5, y = 0.96,
          gp = gpar(fontsize = 16, fontface = "bold", col = "#2c3e50"))
# 关闭绘图设备
dev.off()

# 表格
library(grid)
library(gridExtra)
library(gtable)  # 添加gtable包

# 你的数据
g_data <- data.frame(
  Ontology=ego1_df$ONTOLOGY,
  ID = ego1_df$ID,
  Description = ego1_df$Description,
  stringsAsFactors = FALSE
)

# 创建一组美观的颜色（10种）
pathway_colors <- rep(c("#B6D0E2","#FED7C3", "#a7dace"),each=6)
# 方法1：使用gridExtra创建精确表格
create_table_plot <- function(data, colors) {
  # 设置主题
  my_theme <- ttheme_minimal(
    core = list(
      fg_params = list(
        hjust = 0, 
        x = 0.05,
        fontsize = 11
      ),
      bg_params = list(
        fill = c(rep(c("white", "#f8f9fa"), length.out = nrow(data)))  # 斑马纹
      )
    ),
    colhead = list(
      fg_params = list(
        fontsize = 12,
        fontface = "bold",
        col = "white"
      ),
      bg_params = list(fill = "#2c3e50")
    ),
    rowhead = list(
      fg_params = list(fontsize = 11)
    )
  )
  # 创建表格
  table_grob <- tableGrob(
    d = data,
    theme = my_theme,
    rows = NULL
  )
  # 在表格左侧添加颜色条
  for (i in 1:nrow(data)) {
    # 创建颜色矩形
    color_rect <- rectGrob(
      x = 0, 
      y = 1 - (i - 0.5)/nrow(data),
      width = 0.02,
      height = 1/nrow(data),
      gp = gpar(fill = colors[i], col = NA)
    )
    # 正确使用gtable_add_grob函数
    table_grob <- gtable::gtable_add_grob(
      table_grob,
      color_rect,
      t = i + 1,  # +1是因为表头占一行
      l = 1,
      b = i + 1,
      r = 1
    )
  }
  # 调整列宽
  table_grob$widths <- unit(c(0.05, 0.15, 0.8), "npc")
  return(table_grob)
}
# 方法1A：更简单的版本，不使用gtable_add_grob
create_simple_table <- function(data, colors) {
  # 设置主题，在主题中直接设置行颜色
  my_theme <- ttheme_default(
    core = list(
      fg_params = list(
        hjust = 0, 
        x = 0.1,
        fontsize = 11
      ),
      bg_params = list(
        # 为每一行设置不同的左侧颜色
        fill = colors
      )
    ),
    colhead = list(
      fg_params = list(
        fontsize = 12,
        fontface = "bold",
        col = "white"
      ),
      bg_params = list(fill = "#2c3e50")
    )
  )
  # 创建表格
  table_grob <- tableGrob(
    d = data,
    theme = my_theme,
    rows = NULL
  )
  # 调整列宽
  table_grob$widths <- unit(c(0.11, 0.15,0.74), "npc")
  return(table_grob)
}
# 创建并显示简单表格
table_g_plot <- create_simple_table(g_data, pathway_colors)
# 显示表格
grid.newpage()
grid.draw(table_g_plot)
# 添加标题
grid.text("GO Pathway Enrichment Analysis",
          x = 0.5, y = 0.95,
          gp = gpar(fontsize = 16, fontface = "bold", col = "#2c3e50"))


########## DEG-KEGG circular visualization #######
setwd('/Users/wuhongbo/Desktop/wjx/data_quality/output/JTM_1_cross')
kegg_df=readRDS('kegg_df.rds')

# === 🌟 开启 cairo_pdf 高清输出 ===
cairo_pdf("KEGG_Circos_Enrichment.pdf", width = 4.5, height = 4.5, onefile = TRUE)

layout(matrix(c(1, 2), nrow = 2, ncol = 1), heights = c(0.85, 0.15))
par(mar = c(1, 1, 1, 1))

# 为KEGG通路选择一组颜色
kegg_colors <- colorRampPalette(c("#FF6B6B", "#4ECDC4", "#45B7D1", "#96CEB4", "#FFEAA7", "#DDA0DD"))(nrow(kegg_df))
col <- kegg_colors

circos.clear()
circos.initialize(sectors = kegg_df$ID, xlim = c(0, 1))

# 设置较小的字体大小
par(cex = 0.8)

# 第一轨道：KEGG ID
circos.track(sectors = kegg_df$ID, ylim = c(0, 1),
             bg.col = col,
             track.height = 0.1, # 进一步减小高度
             panel.fun = function(x, y) {
               circos.text(x = CELL_META$xcenter, y = CELL_META$ycenter,
                           labels = CELL_META$sector.index,
                           cex = 0.7,
                           facing = "inside",
                           adj = c(0.5, 0.5))
             })

# 第二轨道：基因计数
circos.track(sectors = kegg_df$ID, ylim = c(0, 1), track.height = 0.12)

for (i in 1:nrow(kegg_df)) {
  # 绘制条形图，去除边框
  circos.rect(xleft = 0, xright = kegg_df$Count[i]/max(kegg_df$Count),
              ybottom = 0.25, ytop = 0.75,
              sector.index = kegg_df$ID[i],
              col = "#a7b2c7",
              border = NA)
  
  # 添加计数标签
  circos.text(x = 0.5, y = 0.5,
              labels = kegg_df$Count[i],
              sector.index = kegg_df$ID[i],
              cex = 0.6,
              col = "black")
  
  # 添加坐标轴
  circos.axis(h = "bottom", major.at = c(0, 1),
              labels = c(0, max(kegg_df$Count)),
              sector.index = kegg_df$ID[i],
              labels.cex = 0.5)
}

# 第三轨道：富集倍数（FC）
circos.track(sectors = kegg_df$ID, ylim = c(0, 1), track.height = 0.1)
for (i in 1:nrow(kegg_df)) {
  circos.text(x = 0.5, y = 0.5,
              labels = round(kegg_df$fc[i], digits = 1),
              sector.index = kegg_df$ID[i],
              cex = 0.6)
}

# 第四轨道：富集倍数条形图
circos.track(sectors = kegg_df$ID, ylim = c(0, ceiling(max(kegg_df$fc))),
             track.height = 0.2)
for (i in 1:nrow(kegg_df)) {
  circos.barplot(value = kegg_df$fc[i], pos = 0.5,
                 sector.index = kegg_df$ID[i],
                 col = "#FED7C3",
                 border = NA)
}

# 第五轨道：-log10 pvalue 条形图
circos.track(sectors = kegg_df$ID, ylim = c(0, ceiling(max(-log10(kegg_df$pvalue)))),
             track.height = 0.2)
for (i in 1:nrow(kegg_df)) {
  # 获取当前pvalue的-log10值
  pval_log10 <- -log10(kegg_df$pvalue[i])
  
  # 绘制条形图，去除边框
  circos.barplot(value = pval_log10, pos = 0.5,
                 sector.index = kegg_df$ID[i],
                 col = col[i],
                 border = NA)
  
  # 在条形图内部显示数值
  text_y <- pval_log10 * 0.8
  
  # 只在条形足够高时显示文本
  if (pval_log10 > max(-log10(kegg_df$pvalue)) * 0.2) {
    circos.text(x = 0.5, y = text_y,
                labels = format(round(pval_log10, 2), nsmall = 2),
                sector.index = kegg_df$ID[i],
                col = "black",
                cex = 0.5,
                facing = "inside",
                adj = c(0.5, 0.5))
  }
}

# === 🌟 底部新增：画从里到外的极简图注 ===
# 恢复默认字体大小
par(cex = 1)

# 重置底部的边距
par(mar = c(0, 0, 0, 0))
plot.new()

legend(
  x = 0.0, y = 0.8, # 🌟 绝对左对齐
  legend = c("", "-log10(pvalue)", "FC", "Count"),
  fill = c(NA, "#FF6B6B", "#FED7C3", "#a7b2c7"), # 内红、中橙、外灰
  border = NA,
  text.font = c(2, 1, 1, 1), # 标题加粗
  horiz = TRUE,
  bty = "n",
  cex = 0.75, # 🌟 字体调小
  x.intersp = 0.5 # 🌟 强行压缩间距
)

# === 🌟 关闭设备输出 PDF ===
dev.off()

# 表格
library(grid)
library(gridExtra)
library(gtable)  # 添加gtable包

# 你的数据
kegg_data <- data.frame(
  ID = kegg_df$ID,
  Description = kegg_df$Description,
  stringsAsFactors = FALSE
)

# 创建一组美观的颜色（10种）
pathway_colors <- kegg_colors

# 方法1：使用gridExtra创建精确表格
create_table_plot <- function(data, colors) {
  # 设置主题
  my_theme <- ttheme_minimal(
    core = list(
      fg_params = list(
        hjust = 0, 
        x = 0.05,
        fontsize = 11
      ),
      bg_params = list(
        fill = c(rep(c("white", "#f8f9fa"), length.out = nrow(data)))  # 斑马纹
      )
    ),
    colhead = list(
      fg_params = list(
        fontsize = 12,
        fontface = "bold",
        col = "white"
      ),
      bg_params = list(fill = "#2c3e50")
    ),
    rowhead = list(
      fg_params = list(fontsize = 11)
    )
  )
  
  # 创建表格
  table_grob <- tableGrob(
    d = data,
    theme = my_theme,
    rows = NULL
  )
  
  # 在表格左侧添加颜色条
  for (i in 1:nrow(data)) {
    # 创建颜色矩形
    color_rect <- rectGrob(
      x = 0, 
      y = 1 - (i - 0.5)/nrow(data),
      width = 0.02,
      height = 1/nrow(data),
      gp = gpar(fill = colors[i], col = NA)
    )
    
    # 正确使用gtable_add_grob函数
    table_grob <- gtable::gtable_add_grob(
      table_grob,
      color_rect,
      t = i + 1,  # +1是因为表头占一行
      l = 1,
      b = i + 1,
      r = 1
    )
  }
  
  # 调整列宽
  table_grob$widths <- unit(c(0.05, 0.15, 0.8), "npc")
  
  return(table_grob)
}

# 方法1A：更简单的版本，不使用gtable_add_grob
create_simple_table <- function(data, colors) {
  # 设置主题，在主题中直接设置行颜色
  my_theme <- ttheme_default(
    core = list(
      fg_params = list(
        hjust = 0, 
        x = 0.1,
        fontsize = 11
      ),
      bg_params = list(
        # 为每一行设置不同的左侧颜色
        fill = colors
      )
    ),
    colhead = list(
      fg_params = list(
        fontsize = 12,
        fontface = "bold",
        col = "white"
      ),
      bg_params = list(fill = "#2c3e50")
    )
  )
  
  # 创建表格
  table_grob <- tableGrob(
    d = data,
    theme = my_theme,
    rows = NULL
  )
  
  # 调整列宽
  table_grob$widths <- unit(c(0.15, 0.85), "npc")
  
  return(table_grob)
}

# 创建并显示简单表格
table_plot <- create_simple_table(kegg_data, pathway_colors)

# 显示表格
grid.newpage()
grid.draw(table_plot)

# 添加标题
grid.text("KEGG Pathway Enrichment Analysis",
          x = 0.5, y = 0.95,
          gp = gpar(fontsize = 16, fontface = "bold", col = "#2c3e50"))
