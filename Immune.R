

library(GSVA)      # 用于ssGSEA分析
library(ggplot2)   # 绘图
fibrosis_geneset <- list(FibrosisSignature = frgs)
param <- ssgseaParam(expr_matrix, fibrosis_geneset)
fibrosis_scores <- gsva(param, verbose=TRUE)
fibrosis_score <- as.numeric(fibrosis_scores)
sample_info <- data.frame(
  SampleID = colnames(expr_matrix),
  FibrosisScore = fibrosis_score,
  Group = GROUP$Group,  # 模拟的真实标签
  ClinicalOutcome = ifelse(GROUP$Group=='IPF',1,0),# 0=健康，1=纤维化（用于ROC分析）
  Project=GROUP$Batch
)
is_ipf <- sample_info$Group == "IPF"
sample_info <- sample_info[is_ipf, ]
expr_matrix <- expr_matrix[, sample_info$SampleID]
sample_info$FibrosisGroup <- ifelse(
  sample_info$FibrosisScore > median(sample_info$FibrosisScore), "High", "Low")

library(GSVA)
library(limma)
library(ggplot2)
library(pheatmap)
library(reshape2)
library(dplyr)
library(tibble)
annotation_col=sample_info
rownames(annotation_col)=annotation_col$SampleID

# 1. 执行ssGSEA分析
param <- ssgseaParam(expr_matrix,
                     cellMarker,
                     minSize = 5, 
                     maxSize = 500) 
ssgsea_scores <- gsva(param, verbose = TRUE)

# write.csv(ssgsea_scores, "ssgsea_scores_results.csv")
# 2. 使用limma进行差异分析
annotation_col$FibrosisGroup <- factor(annotation_col$FibrosisGroup,
                                       levels = c("Low", "High"))
design <- model.matrix(~ 0 + FibrosisGroup, data = annotation_col)
colnames(design) <- gsub("FibrosisGroup", "", colnames(design))
fit <- lmFit(ssgsea_scores, design)
contrast.matrix <- makeContrasts(high_vs_low = High - Low, levels = design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
results <- topTable(fit2, number = Inf, adjust.method = "BH", sort.by = "p")

# 3. 筛选显著结果 (p < 0.01)和top10结果
significant_cells <- results[results$adj.P.Val < 0.01, ]

if (nrow(significant_cells) > 0) {
  cat("\nSignificant immune cell types:\n")
  print(significant_cells[, c("logFC", "P.Value", "adj.P.Val")])
  # 分类上调和下调的细胞类型
  upregulated <- rownames(significant_cells[significant_cells$logFC > 0, ])
  downregulated <- rownames(significant_cells[significant_cells$logFC < 0, ])
  cat("\nUpregulated in high fibrosis group (", length(upregulated), "):\n")
  print(upregulated)
  cat("\nDownregulated in high fibrosis group (", length(downregulated), "):\n")
  print(downregulated)
}
top10_cells <- head(results[order(results$adj.P.Val), ], 10)
top_cells <- rownames(top10_cells)

# 4. 可视化
# 4.1 热图
library(Cairo)
sample_order <- order(-annotation_col$FibrosisScore)
annotation_df <- data.frame(
  "Fibrosis Level" = annotation_col$FibrosisGroup[sample_order],
  row.names = colnames(ssgsea_scores)[sample_order],
  check.names = FALSE
)
plot_matrix <- ssgsea_scores[
  rownames(top10_cells),
  sample_order
]
ann_colors <- list(
  "Fibrosis Level" = c(
    "Low" = "#93cdd9",
    "High" = "#fa6647"
  )
)
# 保存 PDF
cairo_pdf("heatmap-t.pdf", width = 6, height = 3)
hm=pheatmap(
  plot_matrix,
  annotation_col = annotation_df,
  annotation_colors = ann_colors,
  scale = "row",
  color = colorRampPalette(
    c("#2874A6", "white", "red")
  )(100),
  show_colnames = FALSE,
  cluster_cols = FALSE,
  clustering_method = "complete",
  fontsize_row = 8,
  border_color = NA
)
dev.off()

heatmap_cell_order <- rownames(plot_matrix)[hm$tree_row$order]

# 4.2 箱线图展示top10免疫细胞
setwd('/Users/wuhongbo/Desktop/wjx/data_quality/JTM_4_ NSCLC')
library(ggpubr)
library(ggplot2)
library(reshape2)
# 1. 准备数据
plot_data <- as.data.frame(t(ssgsea_scores[heatmap_cell_order, ]))
plot_data$Sample <- rownames(plot_data)
plot_data$FibrosisGroup <- annotation_col[plot_data$Sample, "FibrosisGroup"]
plot_data_long <- melt(
  plot_data,
  id.vars = c("Sample", "FibrosisGroup"),
  variable.name = "CellType",
  value.name = "EnrichmentScore"
)
# 2. 按 logFC 排序并设置因子水平
plot_data_long$CellType <- factor(plot_data_long$CellType, levels = rev(heatmap_cell_order))
logFC_vals <- top10_cells[
  match(heatmap_cell_order, rownames(top10_cells)),
  "logFC"
]
cell_type_colors <- ifelse(logFC_vals > 0, "#C0392B", "#2980B9")
cell_type_colors <- rev(cell_type_colors)
# 4. 绘图
p <- ggplot(
  plot_data_long,
  aes(x = CellType, y = EnrichmentScore, fill = FibrosisGroup)
) +
  geom_boxplot(
    width = 0.7, 
    outlier.shape = NA, 
    linewidth = 0.5,
    position = position_dodge(width = 0.8) 
  ) +
  geom_jitter(
    aes(color = FibrosisGroup), 
    size = 0.6,      
    alpha = 0.3,     
    show.legend = FALSE,
    position = position_jitterdodge(
      jitter.width = 0.2, 
      dodge.width = 0.8
    )
  ) +
  coord_flip() +
  scale_fill_manual(values = c("Low" = "#93cdd9", "High" = "#fa6647")) +
  scale_color_manual(values = c("Low" = "#93cdd9", "High" = "#fa6647")) +
  # FDR 校正
  stat_compare_means(
    aes(group = FibrosisGroup), 
    method = "wilcox.test",
    label = "p.signif", 
    hide.ns = TRUE, 
    size = 4,
    p.adjust.method = "fdr"
  ) +
  labs(x = NULL, y = "ssGSEA enrichment score", fill = "Fibrosis") +
  theme_classic(base_size = 12) +
  theme(
    axis.text.y = element_text(color = cell_type_colors),
    axis.text.x = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    legend.position = "top"
  )
print(p)
ggsave(
  "immune_infiltration_boxplot.pdf", 
  p, 
  width = 5.5, 
  height = 4, 
  device = cairo_pdf
)

# 4.3热图 加辐射图-Fibrosis Group ###
library(Hmisc)
library(Cairo)
library(linkET)
setwd('/Users/wuhongbo/Desktop/wjx/data_quality/JTM_4_ NSCLC')
# 数据准备
ssgsea_df <- as.data.frame(t(ssgsea_scores[heatmap_cell_order, ]))
ssgsea_df <- ssgsea_df[, colSums(is.na(ssgsea_df)) == 0]
ssgsea_mat <- as.matrix(ssgsea_df)
cell_names <- colnames(ssgsea_mat)
# 只保留 High
high_vec <- as.numeric(sample_info$FibrosisGroup == "High")
library(vegan)
# Mantel dataframe
mantel_df <- data.frame(
  var1 = rep("Fibrosis Program", length(cell_names)),
  var2 = cell_names,
  r_value = NA,
  p_value = NA,
  stringsAsFactors = FALSE
)
# Mantel test: Fibrosis group vs immune cells
fibrosis_vec <- as.numeric(sample_info$FibrosisGroup == "High")
fibrosis_dist <- dist(fibrosis_vec, method = "euclidean")
set.seed(123)
for (i in seq_along(cell_names)) {
  y <- ssgsea_mat[, i]
  cell_dist <- dist(y, method = "euclidean")
  mt <- mantel(
    fibrosis_dist,
    cell_dist,
    method = "spearman",
    permutations = 9999
  )
  mantel_df$r_value[i] <- mt$statistic
  mantel_df$p_value[i] <- mt$signif
}
# 分箱
mantel_df$r_group <- cut(
  abs(mantel_df$r_value),
  breaks = c(-Inf, 0.1, 0.2, Inf),
  labels = c("<0.1", "0.1–0.2", ">0.2"),
  right = FALSE
)
mantel_df$p_group <- cut(
  mantel_df$p_value,
  breaks = c(-Inf, 0.001, 0.05, Inf),
  labels = c("<0.001", "0.001–0.05", ">0.05"),
  right = TRUE
)
# correlation matrix
corr_mat <- rcorr(ssgsea_mat, type = "spearman")
# 画图
p_mantel <- qcorrplot(
  corr_mat,
  type = "upper",
  diag = FALSE,
  grid_size = 0.4,
  grid_col = "lightgray"
) +
  geom_square(linetype = 0) +
  geom_couple(
    aes(colour = p_group, size = r_group),
    data = subset(mantel_df, p_value < 0.05),
    curvature = 0.1
  ) +
  set_corrplot_style(
    colours = c("#FF8040", "white", "#5BC2CD")
  ) +
  scale_size_manual(
    name = "Mantel's |r|",
    values = c(
      "<0.1" = 0.2,
      "0.1–0.2" = 0.5,
      ">0.2" = 0.8
    )
  ) +
  scale_colour_manual(
    name = "Mantel's p",
    values = c(
      ">0.05" = "gray",
      "0.001–0.05" = "#9ACD32",
      "<0.001" = "#87CEEB"
    )
  ) +
  guides(
    size = guide_legend(
      title = "Mantel's |r|",
      override.aes = list(colour = "grey35"),
      order = 2
    ),
    colour = guide_legend(
      title = "Mantel's p",
      override.aes = list(size = 2),
      order = 1
    )
  ) +
  theme(
    axis.text.x.top = element_text(
      angle = 45,
      hjust = 0,
      vjust = 0
    ),
    plot.margin = margin(t = 10, r = 10, b = 10, l = 80)
  ) +
  coord_cartesian(clip = "off")

cairo_pdf("immune_corrplot.pdf", width = 8.5, height = 5.5)
print(p_mantel)
dev.off()






