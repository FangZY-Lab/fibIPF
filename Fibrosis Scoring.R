# 1. 加载R包
library(GSVA) # 用于ssGSEA分析
library(ggplot2) # 绘图
library(pROC) # ROC分析
library(ggpubr) # 美化图形
library(factoextra)# PCA可视化
library(pheatmap) # 热图绘制
library(survival) # 生存分析（如果需要）
library(survminer) # 生存分析可视化

# 2. 计算FibrosisScore (ssGSEA)
fibrosis_geneset <- list(FibrosisSignature = frgs)
param <- ssgseaParam(expr_matrix, fibrosis_geneset)
fibrosis_scores <- gsva(param, verbose=TRUE)
fibrosis_score <- as.numeric(fibrosis_scores)

# 将分数与样本信息结合
sample_info <- data.frame(
  SampleID = colnames(expr_matrix),
  FibrosisScore = fibrosis_score,
  Group = GROUP$Group,
  ClinicalOutcome = ifelse(GROUP$Group=='IPF',1,0),
  Project=GROUP$Batch
)

# 3. 风险分层（按中位数或最佳cutoff）
sample_info$FibrosisLevel <- ifelse(sample_info$FibrosisScore > median(sample_info$FibrosisScore),
                                    "High", "Low")

# 方法2：使用survminer包寻找最佳cutoff（如果有生存数据）
if(FALSE){ # 假设有生存数据时的代码
  sample_info$OS <- rnorm(n_samples, mean = 60, sd = 20)
  sample_info$Status <- rbinom(n_samples, 1, 0.3)
  cutoff <- surv_cutpoint(sample_info,
                          time = "OS",
                          event = "Status",
                          variables = "FibrosisScore")
  sample_info$FibrosisLevel <- surv_categorize(cutoff)$FibrosisScore
}

# 4. 主成分分析(PCA)
pca_data <- t(expr_matrix[frgs, ])
pca_result <- prcomp(pca_data, scale = TRUE, center = TRUE)
var_explained <- pca_result$sdev^2 / sum(pca_result$sdev^2)

# 创建PCA结果数据框
pca_df <- data.frame(
  PC1 = pca_result$x[, 1],
  PC2 = pca_result$x[, 2],
  FibrosisLevel = sample_info$FibrosisLevel,
  TrueGroup = sample_info$Group,
  FibrosisScore = sample_info$FibrosisScore
)

# 5. ROC曲线分析
roc_result <- roc(response = sample_info$ClinicalOutcome,
                  predictor = sample_info$FibrosisScore,
                  levels = c(0, 1),
                  direction = "<")
auc_value <- auc(roc_result)
ci_auc <- ci.auc(roc_result)

# 6. 可视化
sample_info$FibrosisScore=sample_info$FibrosisScore-min(sample_info$FibrosisScore)

# 6.1 FibrosisScore小提琴图 (新增散点与显著性检验)
p1 <- ggplot(sample_info, aes(x = Group, y = FibrosisScore, fill = Group)) +
  geom_violin(alpha = 0.7, trim = FALSE, linewidth = 0.5) +
  geom_jitter(width = 0.12, size = 1.0, alpha = 0.4, color = "grey20", shape = 16) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA, alpha = 1, linewidth = 0.6) +
  stat_compare_means(comparisons = list(c("Control", "IPF")),
                     method = "wilcox.test",
                     label = "p.signif",
                     size = 5) + # 星号稍微调大一点点
  scale_fill_manual(values = c("IPF" = "#FFBF00", "Control" = "#B6D0E2")) +
  #  顶部增加留白，防止星号贴顶
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.12))) +
  labs(title = " ",
       x = "Group", y = "Fibrosis Score") +
  theme_classic() +
  theme(legend.position = "none",
        
        axis.text = element_text(size = 13, color = "black"), # 刻度文字纯黑
        axis.title = element_text(size = 14),
        # 坐标轴骨架加粗，更有质感
        axis.line = element_line(linewidth = 0.6, color = "black"),
        axis.ticks = element_line(linewidth = 0.6, color = "black"))

print(p1)

ggsave("FibrosisScore_Violin6.pdf", plot = p1, width = 2.7, height = 2.7, device = cairo_pdf)


# 点点线 分布图
library(dplyr)
library(ggthemes)
sample_info_sorted <- sample_info %>%
  arrange(FibrosisScore) %>%
  mutate(SampleIndex = 1:n(),
         FibrosisLevel = factor(FibrosisLevel, levels = c("Low", "High")))

low_max <- sample_info_sorted %>%
  filter(FibrosisLevel == "Low") %>%
  summarise(max_score = max(FibrosisScore)) %>%
  pull(max_score)

cutoff_index <- sample_info_sorted %>%
  filter(FibrosisLevel == "Low", FibrosisScore == low_max) %>%
  summarise(max_index = max(SampleIndex)) %>%
  pull(max_index)

p <- ggplot(sample_info_sorted, aes(x = SampleIndex, y = FibrosisScore, color = FibrosisLevel)) +
  # 1. 稍微缩小点的大小，增加透明度，体现数据的海量感
  geom_point(size = 1.5, alpha = 0.4) +
  # 2. 虚线调细一点，改为深灰色，不抢散点的风头
  geom_vline(xintercept = cutoff_index, linetype = "dashed", color = "grey30", linewidth = 0.6, alpha = 0.8) +
  geom_hline(yintercept = low_max, linetype = "dashed", color = "grey30", linewidth = 0.6, alpha = 0.8) +
  # 3. Cutoff 定位十字架
  annotate("point", x = cutoff_index, y = low_max, size = 4, shape = 3, color = "purple", stroke = 1.2) +
  scale_color_manual(values = c("Low" = "#93cdd9", "High" = "#fa6647")) +
  scale_x_continuous(name = "Sample (sorted by Fibrosis Score)", expand = expansion(mult = c(0.02, 0.02))) +
  scale_y_continuous(name = "Fibrosis Score", limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
  labs(title = " ", color = "Fibrosis Level") +
  theme_few() +
  theme(
    # 4. 字体大小适中，粗细有致
    plot.title = element_text(hjust = 0.5, size = 12),
    # 5. 图例精调：位置微调，去除背景色防止遮挡边框！
    legend.position = c(0.2, 0.85),
    legend.background = element_blank(), # 核心：透明背景不遮挡线条
    legend.key = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 11, color = "black"),
    axis.title = element_text(size = 13)
  )

print(p)

ggsave("FibrosisScore_Distrib.pdf", plot = p, width = 3.3, height = 3, device = cairo_pdf)



# 6.2 顶刊级 PCA 散点图
p2_premium <- ggplot(pca_df, aes(x = PC1, y = PC2, color = FibrosisLevel)) +
  # 优化1：稍微缩小点(1.2)，增加一点透明度(0.7)，让重叠的地方也能看出层次感
  geom_point(aes(shape = TrueGroup), size = 1.2, alpha = 0.7) +
  
  # 优化2：保留淡雅底色和外圈实线
  stat_ellipse(aes(fill = FibrosisLevel), geom = "polygon", alpha = 0.12, 
               level = 0.95, type = "norm", show.legend = FALSE) +
  stat_ellipse(level = 0.95, type = "norm", linewidth = 0.4, show.legend = FALSE) +
  
  scale_color_manual(values = c("High" = "#fa6647", "Low" = "#93cdd9")) +
  scale_fill_manual(values = c("High" = "#fa6647", "Low" = "#93cdd9")) +
  labs(x = sprintf("PC1 (%.1f%%)", var_explained[1]*100),
       y = sprintf("PC2 (%.1f%%)", var_explained[2]*100),
       color = "Fibrosis Level",
       shape = "True Group") +
  theme_classic() +
  theme(
    
    # 优化4：精调图例的字体排版
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 12),
    legend.position = "right",
    legend.background = element_blank(),
    legend.key = element_blank()
  )

print(p2_premium)

ggsave("PCA_Analysis6.pdf", plot = p2_premium, width = 4, height = 2.5, device = cairo_pdf)

# 6.3 ROC曲线
p3 <- ggroc(roc_result, color = "#4DAF4A", size = 1.2) +
  geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1), color = "grey", linetype = "dashed") +
  annotate("text", x = 0.5, y = 0.3,
           label = sprintf("AUC = %.3f
(95%% CI: %.3f-%.3f)", auc_value, ci_auc[1], ci_auc[3]),size = 5) +
  labs(title = " ",
       x = "1 - Specificity", y = "Sensitivity") +
  theme_classic()+
  theme(axis.title = element_text(size = 14),axis.text = element_text(size = 14))

ggsave("ROC.pdf", plot = p3, width = 3.4, height = 3.4, device = cairo_pdf)

#6.4 热图展示11个DE-FRGs表达
sample_order <- order(sample_info$FibrosisScore, decreasing = FALSE)
heatmap_data <- expr_matrix[frgs, sample_order]

sample_info$Project <- factor(sample_info$Project,
                              levels = c("1", "2", "3", "4", "5", "6"),
                              labels = c("GSE92592", "GSE124685", "GSE150910",
                                         "GSE184316", "GSE213001", "GSE231693"))
# 注释条名称变为 FibrosisLevel
annotation_col <- data.frame(
  FibrosisLevel = sample_info$FibrosisLevel[sample_order],
  FibrosisScore = sample_info$FibrosisScore[sample_order],
  TrueGroup = sample_info$Group[sample_order],
  Project=sample_info$Project[sample_order]
)
rownames(annotation_col) <- colnames(heatmap_data)

# 2. 强行重命名列名（加上完美的空格）
colnames(annotation_col) <- c("Fibrosis Level", "Fibrosis Score", "True Group", "Project")

# 3. 颜色字典的名称必须和上面带空格的名字完全一一对应（加上引号）
annotation_colors <- list(
  "True Group" = c("Control" = "#B6D0E2", "IPF" = "#FFBF00"),
  "Fibrosis Level" = c("High" = "#fa6647", "Low" = "#93cdd9"),
  "Fibrosis Score" = colorRampPalette(c("white", "#4daf4a"))(100),
  "Project" = c("GSE92592" = "#8DD3C7", "GSE124685" = "#FFFFB3", "GSE150910" = "#BEBADA",
                "GSE184316" = "#FB8072", "GSE213001" = "#80B1D3", "GSE231693" = "#FDB462")
)

p4 <- pheatmap(heatmap_data,
               scale = "row",
               color = colorRampPalette(c("blue", "white", "red"))(100),
               cluster_rows = TRUE, cluster_cols = FALSE,
               show_colnames = FALSE,
               annotation_col = annotation_col,
               annotation_colors = annotation_colors,
               main = " ",
               fontsize_row = 10, border_color = NA)

ggsave("heatmap6.pdf", plot = p4, width = 8, height = 4.7, device = cairo_pdf)


# 6.5 箱线图
library(reshape2)
FC_seq=c("COMP", "CXCL14", "THY1", "CTHRC1", "COL1A1","POSTN",
         "COL1A2", "IGF1", "CXCL12", "TGFB3","ADAM12")
# 过滤掉矩阵中没有的基因，避免报错
FC_seq <- intersect(FC_seq, rownames(heatmap_data))

exp_t <- as.data.frame(t(heatmap_data[FC_seq,sample_info$SampleID]))
exp_t$SampleID <- rownames(exp_t)

plot_data <- merge(exp_t, sample_info[, c("SampleID", "FibrosisLevel")], by = "SampleID")
plot_data_long <- melt(plot_data, id.vars = c("SampleID", "FibrosisLevel"),
                       variable.name = "Gene", value.name = "Expression")

plot_data_long$FibrosisLevel <- factor(plot_data_long$FibrosisLevel, levels = c("Low", "High"))

# 赋值给 p6 变量
p6 <- ggplot(plot_data_long, aes(x = Gene, y = Expression, fill = FibrosisLevel)) +
  geom_boxplot(outlier.size = 0.5, alpha = 0.8) +
  stat_compare_means(aes(group = FibrosisLevel), method = "wilcox.test",
                     label = "p.signif", hide.ns = TRUE, vjust = 0.5) +
  scale_fill_manual(values = c("Low" = "#93cdd9", "High" = "#fa6647")) +
  labs(x = " ", y = "Expression Level",
       title = " ",
       fill = "Fibrosis Level") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10, color = "black"), # 纯黑字体
    axis.text.y = element_text(size = 10, color = "black"), # 纯黑字体
    axis.title = element_text(size = 10),
    plot.title = element_text(size = 10, hjust = 0.5),
    legend.position = "top",
    legend.title = element_text(size = 10),
    panel.grid.major.x = element_blank(),
    panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5)
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))

print(p6)

# 使用 cairo_pdf 保存，宽度设为 8，高度设为 4.5（适合长条形的基因排布）
ggsave("FRGs_Boxplot.pdf", plot = p6, width = 6, height = 3.5, device = cairo_pdf)





# 7. 统计摘要输出
cat("=== FibrosisScore Model Summary ===
")
cat(sprintf("Total Samples: %d
", nrow(sample_info)))
cat(sprintf("High Fibrosis: %d, Low Fibrosis: %d
",
            sum(sample_info$FibrosisLevel == "High"),
            sum(sample_info$FibrosisLevel == "Low")))
cat(sprintf("FibrosisScore Range: %.3f - %.3f
",
            min(sample_info$FibrosisScore), max(sample_info$FibrosisScore)))
cat(sprintf("ROC AUC: %.3f (95%% CI: %.3f-%.3f)
",
            auc_value, ci_auc[1], ci_auc[3]))

# 8. 保存结果
supp_table2 <- data.frame(
  SampleID = sample_info$SampleID,
  FibrosisScore = round(sample_info$FibrosisScore, 3),
  FibrosisLevel = sample_info$FibrosisLevel,
  TrueGroup = sample_info$Group,
  ClinicalOutcome = sample_info$ClinicalOutcome
)
supp_table2 <- cbind(supp_table2, as.data.frame(pca_result$x[, 1:5]))
write.csv(supp_table2, "Supplementary_Table2_FibrosisScore_Results.csv", row.names = FALSE)

ggsave("FibrosisScore_Distribution.png", p1, width = 8, height = 6)
ggsave("PCA_Analysis.png", p2, width = 10, height = 6)
ggsave("ROC_Curve.png", p3, width = 8, height = 6)

