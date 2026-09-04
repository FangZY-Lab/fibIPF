## LASSO #######################
# 1. 加载R包
library(glmnet)
library(ggplot2)
library(dplyr)

# 2. 定义初始向量
x <- t(expr_data)  # 维度：534个样本 × 13108个基因
y <- ifelse(GROUP$Group == "IPF", 1, 0)  # 二分类响应变量（1=IPF，0=对照）
y <- as.factor(y)  # 分类变量格式

# 3. 构建Lasso回归模型（含交叉验证选择最优λ）
set.seed(123)  # 设置随机种子保证结果可重复
cv_lasso <- cv.glmnet(
  x = x, 
  y = y, 
  family = "binomial",  # 二分类问题选择binomial，连续变量选gaussian
  alpha = 1,  # alpha=1为Lasso回归，alpha=0为Ridge回归
  nfolds = 10,  # 10折交叉验证
  type.measure = "deviance"  # 损失函数：偏差（分类问题适用）
)
best_lasso <- glmnet(
  x = x, 
  y = y, 
  family = "binomial", 
  alpha = 1, 
  lambda = cv_lasso$lambda.min  # 也可替换为cv_lasso$lambda.1se
)

# 4. 筛选重要生物标志物（系数非零的基因）
lasso_coef <- coef(best_lasso)  # 提取模型系数
biomarkers <- as.data.frame(as.matrix(lasso_coef)) %>%
  tibble::rownames_to_column("Gene") %>%
  rename(Coef = s0) %>%
  filter(Gene != "(Intercept)" & Coef != 0) %>%
  arrange(desc(abs(Coef)))

# 5. Lasso系数路径图可视化

# 提取系数路径数据
library(stringr)
lasso_path_data <- as.data.frame(as.matrix(cv_lasso$glmnet.fit$beta)) %>%
  tibble::rownames_to_column("Gene") %>%
  tidyr::pivot_longer(-Gene, names_to = "Lambda_Index", values_to = "Coef") %>%
  mutate(
    Lambda_Index_Num = as.integer(stringr::str_extract(Lambda_Index, "\\d+")),
    Lambda = cv_lasso$glmnet.fit$lambda[Lambda_Index_Num + 1],  # +1是因为R索引从1开始，glmnet列名从s0（对应索引1）
    logLambda = log(Lambda)
  )

# 绘制ggplot2版本系数路径图
# 定义美观的颜色方案
gene_colors <- c("#377EB8", "#89a67a", "#984EA3", "#FF7F00",
                 "#Fdd538", "#A65628", "#Fa8dbe", "#999999", "#33b39f")

p1 <- ggplot(lasso_path_data, aes(x = logLambda, y = Coef, group = Gene)) +
  geom_line(alpha = 0.2, color = "gray60", linewidth = 0.3) +
  geom_line(data = subset(lasso_path_data, Gene %in% biomarkers$Gene[1:9]),
            aes(color = Gene), linewidth = 0.7) +
  scale_color_manual(values = gene_colors) +
  geom_vline(xintercept = log(cv_lasso$lambda.min), 
             linetype = "dashed", color = "black", linewidth = 0.4, alpha = 0.8) +
  geom_vline(xintercept = log(cv_lasso$lambda.1se), 
             linetype = "dashed", color = "#0072B2", linewidth = 0.4, alpha = 0.8) +
  annotate("text", x = log(cv_lasso$lambda.min) + 0.5, y = max(lasso_path_data$Coef),
           label = "lambda.min", color = "black", hjust = 0, 
           size = 4) +
  annotate("text", x = log(cv_lasso$lambda.1se) + 0.5, y = max(lasso_path_data$Coef) - 1,
           label = "lambda.1se", color = "#0072B2", hjust = 0,
           size = 4) +
  labs(title = "Lasso Coefficient Path",
       x = "log(λ)",
       y = "Coefficient",
       color = " ") +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 10),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 10),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = ),
    legend.position = "right",
    panel.grid.major = element_line(color = "gray90", linewidth = 0.3),
    panel.grid.minor = element_blank()
  )
ggsave("Lasso_Coefficient_Path_ggplot.png", p1, width = 12, height = 8, dpi = 300)

# 6. 交叉验证误差曲线可视化
plot(cv_lasso,
     main = " ",
     xlab = "log(λ)",
     ylab = "Binomial Deviance",
     mgp = c(2.2, 0.8, 0))
points(log(cv_lasso$lambda), cv_lasso$cvm, 
       col = "#e74b35", type = "b", pch = 19, cex = 0.9)

## common ##############
# 加载必需包
library(VennDiagram)
library(grid)

# 1. 定义基因集合（你的原始集合）
set1 <- biomarkers$Gene
set2 <- rownames(gene_importance_df)[1:10]
gene_list <- list(
  "LASSO" = set1,  # 第一个集合对应brown
  "RF" = set2      # 第二个集合对应#984EA3
)

p2 <- venn.diagram(
  x = gene_list, # 向量集合
  scaled = F, # 根据比例显示大小
  category.names = c("RF\n(10)","LASSO\n(9)" ),
  fill = c("#7f7f7f", "#6f6db9"), 
  cat.col = c("#6f6db9", "#7f7f7f"),   # 设置类别标签的颜色
  alpha=0.5,  # 着色透明度
  col="black",
  output = FALSE ,
  disable.logging = F,
  lwd = 1.3, # 外圈的颜色和粗细
  lty = 1,
  cex = 1.4, # 圈内字体大小
  cat.cex = 1.4, # 外圈标签字体大小
  cat.default.pos = "outer",
  cat.dist = 0.15,        # 圈的半径
  cat.pos = c(-140, 140),  # 调整标签位置,顺时钟360度
  output=F,
  filename=NULL,
  imagetype="pdf" ,
  height = 180 , 
  width = 180 , 
  resolution = 200,
  margin  = 0.25 # 圈外标签与图形边距的距离
)
grid.newpage()
grid.draw(p2)



## RF Random Forest全部基因 ############
set.seed(12345)
library(randomForest)
library(ggplot2)
library(dplyr)

datadir='/Users/wuhongbo/Desktop/wjx/data_quality/output'
setwd(datadir)
exp=readRDS('Batch_limma_data.rds')
GROUP=readRDS('combined_group.rds')

frgs=c("COMP","CTHRC1","CXCL14","THY1","COL1A1", 
       "POSTN","IGF1","COL1A2","CXCL12","ADAM12","TGFB3")
expr_data=exp[,GROUP$Sample]#全部基因
expr_data <- t(expr_data)
sample_groups=GROUP$Group

rf_model <- randomForest(as.factor(sample_groups) ~ .,
                         data=expr_data,
                         ntree=500,
                         importance=TRUE)
rf_model
plot(rf_model, main="Random Forest Error Rate", lwd=2,
     col=c("black", "#377EB8", "#e74b35"))
legend("topright",
       legend=c("OOB", "Control", "IPF"),
       col=c("black", "#377EB8", "#e74b35"),
       lty=c(1, 2, 3),  # 分别为实线、虚线、点线
       lwd=2) 

optimal_ntree <- which.min(rf_model$err.rate[, 1])
rf_optimal <- randomForest(as.factor(sample_groups) ~ .,
                           data=expr_data,
                           ntree=optimal_ntree,
                           importance=TRUE,
                           proximity=TRUE)
gene_importance <- importance(rf_optimal)


gene_importance_df <- as.data.frame(gene_importance) %>%
  mutate(Gene = rownames(.)) %>%
  arrange(desc(MeanDecreaseAccuracy))
top_genes_accuracy <- gene_importance_df[1:11, ]
lollipop_plot <- ggplot(top_genes_accuracy, 
                        aes(x = reorder(Gene, MeanDecreaseAccuracy),  #按重要性降序排序x轴 
                            y = MeanDecreaseAccuracy)) +
  geom_segment(aes(xend = Gene, yend = 0),color = "#2c3e50", 
               size = 1,alpha = 0.7) +
  geom_point(aes(color = MeanDecreaseAccuracy), 
             size = 4,alpha = 0.8,shape = 16) +
  scale_color_gradientn(colors = c("#add8e6", "#6f6db9", "#00008b"),  # 红→紫→蓝
                        name = "Mean Decrease\nAccuracy") +
  labs(title = "Top Genes by MeanDecreaseAccuracy",
       subtitle = "Fibrosis-related Hub Genes",
       x = "Gene Name",
       y = "Mean Decrease Accuracy") +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 12, color = "#7f8c8d"),
    axis.title.x = element_text(size = 11),
    axis.title.y = element_text(size = 11),
    axis.text.x = element_text(size = 9),  # 旋转x轴标签，避免重叠
    panel.grid.minor = element_blank()  # 隐藏次要网格线
  )+
  coord_flip()

lollipop_plot





## RF Random Forest frgs ############
set.seed(12345)
library(randomForest)
library(ggplot2)
library(dplyr)

datadir='/Users/wuhongbo/Desktop/wjx/data_quality/output'
setwd(datadir)
exp=readRDS('Batch_limma_data.rds')
GROUP=readRDS('combined_group.rds')

frgs=c("COMP","CTHRC1", "CXCL14", "THY1", "COL1A1", 
       "POSTN","IGF1", "COL1A2", "CXCL12", "ADAM12","TGFB3")
expr_data=exp[frgs,GROUP$Sample]
expr_data <- t(expr_data)
sample_groups=GROUP$Group

rf_model <- randomForest(as.factor(sample_groups) ~ .,
                         data=expr_data,
                         ntree=500,
                         importance=TRUE)
rf_model
plot(rf_model, main="Random Forest Error Rate", lwd=2,
     col=c("black", "#3498db", "#e74c3c"))
legend("topright",
       legend=c("OOB", "Control", "IPF"),
       col=c("black", "#3498db", "#e74c3c"),
       lty=c(1, 2, 3),  # 分别为实线、虚线、点线
       lwd=2) 

optimal_ntree <- which.min(rf_model$err.rate[, 1])
rf_optimal <- randomForest(as.factor(sample_groups) ~ .,
                           data=expr_data,
                           ntree=optimal_ntree,
                           importance=TRUE,
                           proximity=TRUE)
gene_importance <- importance(rf_optimal)


gene_importance_df <- as.data.frame(gene_importance) %>%
  mutate(Gene = rownames(.)) %>%
  arrange(desc(MeanDecreaseAccuracy))
top_genes_accuracy <- gene_importance_df[1:11, ]
lollipop_plot <- ggplot(top_genes_accuracy, 
                        aes(x = reorder(Gene, MeanDecreaseAccuracy),  #按重要性降序排序x轴 
                            y = MeanDecreaseAccuracy)) +
  geom_segment(aes(xend = Gene, yend = 0),color = "#2c3e50", 
               size = 1,alpha = 0.7) +
  geom_point(aes(color = MeanDecreaseAccuracy), 
             size = 4,alpha = 0.8,shape = 16) +
  scale_color_gradientn(colors = c("#6495ED", "#9b59b6", "#e74c3c"),  # 红→紫→蓝
                        name = "Mean Decrease\nAccuracy") +
  labs(title = "Top Genes by MeanDecreaseAccuracy",
       subtitle = "Fibrosis-related Hub Genes",
       x = "Gene Name",
       y = "Mean Decrease Accuracy") +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 12, color = "#7f8c8d"),
    axis.title.x = element_text(size = 11),
    axis.title.y = element_text(size = 11),
    axis.text.x = element_text(size = 9),  # 旋转x轴标签，避免重叠
    panel.grid.minor = element_blank()  # 隐藏次要网格线
  )+
  coord_flip()

lollipop_plot




## SVM ###############
datadir='/Users/wuhongbo/Desktop/wjx/data_quality/output'
setwd(datadir)
exp=readRDS('Batch_limma_data.rds')
GROUP=readRDS('combined_group.rds')

frgs=c("COMP","CTHRC1","CXCL14","THY1","COL1A1", 
       "POSTN","IGF1","COL1A2","CXCL12","ADAM12","TGFB3")
expr_data=exp[,GROUP$Sample]

library(e1071)
library(kernlab)
library(caret)
library(ggplot2)

set.seed(123)
data=t(exp[,GROUP$Sample])
group=GROUP$Group
#feature_sizes <- c(2, 4, 6, 8, seq(10, 40, by=3))
feature_sizes <- c(2:11)

Profile <- rfe(
  x = data,                              # 特征矩阵
  y = as.numeric(as.factor(group)),      # 标签
  sizes = feature_sizes,                 # 测试的特征数
  rfeControl = rfeControl(
    functions = caretFuncs,            # 使用caret函数
    method = "cv",                     # 交叉验证
    number = 10                        # 10折
  ),
  method = "svmRadial"                   # RBF核SVM
)

# 最优特征数
optSize <- Profile$optsize

# 筛选的基因
optVariables <- Profile$optVariables



## common ########
setwd('/Users/wuhongbo/Desktop/wjx/data_quality/output')
setwd('./JTM_3_HubGene')
library(ggVennDiagram)
library(ggplot2)
library(UpSetR)

# set1 <- optVariables
# set2 <- rownames(gene_importance_df)[1:10]
# set3=biomarkers$Gene
# gene_sets <- list("SVM" = set1,"RF" = set2,"LASSO" = set3)
# saveRDS(gene_sets,file = 'hub_gene_list.rds')
gene_sets=readRDS('hub_gene_list.rds')
p=ggVennDiagram(gene_sets, 
                set_color = c("#a6c05c", "#d17995",'#bf9f92'),  
                label = "count") +  
  ggplot2::scale_fill_gradient(low = "white", high = "white") +  
  ggplot2::theme(plot.title = element_text(hjust = 0.5, size = 14),
                 plot.margin = margin(l = 20, r = 20, unit = "pt") ) +
  ggplot2::guides(fill = "none", color = "none")  

#ggsave("venn.png", width = 8, height = 6, dpi = 300)
cairo_pdf("venn.pdf", width = 3, height = 3)
print(p)
dev.off()

HubGene=intersect(set1,intersect(set2,set3))



#### hub gene 在免疫中的作用 ################

datadir='/Users/wuhongbo/Desktop/wjx/data_quality/output'
setwd(datadir)
expr_matrix=readRDS('Batch_limma_data.rds')
GROUP=readRDS('combined_group.rds')



library(GSVA)
library(ggplot2)
library(ggpubr)
library(pheatmap)
library(circlize)
library(ComplexHeatmap)
library(reshape2)
library(dplyr)
library(tidyverse)
library(ggcorrplot)



set.seed(123)
hub_genes <- c("CTHRC1", "COMP", "CXCL14")
target_expr <- expr_matrix[hub_genes, ]


# 这里使用Bindea等人2013年定义的28种免疫细胞特征
# 读取免疫特征基因集
load("~/Desktop/wjx/data_quality/input/ssGSEA28.Rdata")
immunity_genesets=cellMarker

# 过滤表达矩阵中存在的基因
immunity_genesets_filtered <- lapply(immunity_genesets, function(x) {
  intersect(x, rownames(expr_matrix))
})

# 执行ssGSEA分析
ssgsea_scores <- gsva(expr_matrix, 
                      immunity_genesets_filtered,
                      method = "ssgsea",
                      kcdf = "Gaussian",  # 对于log2转换的FPKM/TPM数据
                      min.sz = 3,         # 最小基因集大小
                      max.sz = 500,       # 最大基因集大小
                      verbose = TRUE,
                      parallel.sz = 1)    # 并行核心数，根据您的机器调整

# 将ssGSEA分数矩阵转置为样本×免疫细胞类型
immune_scores <- t(ssgsea_scores)



# 相关性分析
# 将目标基因表达与免疫浸润分数合并
target_expr_t <- t(target_expr)  # 转置为样本×基因
colnames(target_expr_t) <- paste0("Expr_", colnames(target_expr_t))

# 合并数据
combined_data <- cbind(target_expr_t, immune_scores)

# 计算相关系数矩阵
cor_matrix <- cor(combined_data, method = "spearman")  # 使用Spearman相关系数

# 提取目标基因与免疫细胞的相关性子矩阵
target_immune_cor <- cor_matrix[
  grep("^Expr_", rownames(cor_matrix)),
  !grepl("^Expr_", colnames(cor_matrix))
]

# 重命名行名
rownames(target_immune_cor) <- gsub("Expr_", "", rownames(target_immune_cor))


library(ComplexHeatmap)
library(circlize)

# 首先计算相关性矩阵和p值矩阵
library(Hmisc)  # 用于计算带p值的相关系数
# 计算相关性矩阵和p值矩阵
cor_result <- rcorr(as.matrix(combined_data), type = "spearman")
# 提取目标基因与免疫细胞的相关性矩阵
target_immune_cor <- cor_result$r[
  grep("^Expr_", rownames(cor_result$r)),
  !grepl("^Expr_", colnames(cor_result$r))
]
# 提取对应的p值矩阵
target_immune_p <- cor_result$P[
  grep("^Expr_", rownames(cor_result$r)),
  !grepl("^Expr_", colnames(cor_result$r))
]
# 重命名行名
rownames(target_immune_cor) <- gsub("Expr_", "", rownames(target_immune_cor))
rownames(target_immune_p) <- gsub("Expr_", "", rownames(target_immune_p))

heatmap_plot <- Heatmap(
  target_immune_cor,
  name = "Correlation",
  col = colorRamp2(c(-1, 0, 1), c("#2874A6", "white", "red")),
  row_names_side = "left",
  row_names_gp = gpar(fontsize = 10),
  column_names_side = "bottom",
  column_names_gp = gpar(fontsize = 10),
  column_names_rot = 45,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  cell_fun = function(j, i, x, y, width, height, fill) {
    p_val <- target_immune_p[i, j]
    if (!is.na(p_val)) {
      if (p_val < 0.001) {
        grid.text("***", x, y, gp = gpar(fontsize = 10, col = "black"))
      } else if (p_val < 0.01) {
        grid.text("**", x, y, gp = gpar(fontsize = 10, col = "black"))
      } else if (p_val < 0.05) {
        grid.text("*", x, y, gp = gpar(fontsize = 10, col = "black"))
      }
    }
  },
  column_title = "Hub FRGs Immune Correlation",
  column_title_gp = gpar(fontsize = 10),
  heatmap_legend_param = list(
    title = "Spearman r",
    title_position = "leftcenter-rot",
    legend_height = unit(3, "cm")
  )
)

print(heatmap_plot)




# 可视化免疫浸润聚类热图
heatmap_data <- immune_scaled
annotation_col <- data.frame(
  Cluster = as.factor(immune_clusters$cluster),
  row.names = rownames(immune_scores)
)

# 为每个聚类添加目标基因表达信息
for(gene in hub_genes) {
  if(gene %in% rownames(expr_matrix)) {
    gene_expr <- as.numeric(expr_matrix[gene, ])
    gene_expr_scaled <- scale(gene_expr)
    annotation_col[[paste0(gene, "_expr")]] <- gene_expr_scaled
  }
}

# 创建热图
cluster_colors <- list(
  Cluster = setNames(RColorBrewer::brewer.pal(k_selected, "Set1"), 
                     1:k_selected)
)

for(gene in hub_genes) {
  if(paste0(gene, "_expr") %in% colnames(annotation_col)) {
    cluster_colors[[paste0(gene, "_expr")]] <- 
      colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
  }
}

# 绘制复杂热图
ha <- HeatmapAnnotation(df = annotation_col, col = cluster_colors)

ht <- Heatmap(t(heatmap_data),
              name = "Z-score",
              top_annotation = ha,
              col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
              show_row_names = TRUE,
              show_column_names = FALSE,
              row_names_gp = gpar(fontsize = 8),
              column_split = immune_clusters$cluster,
              cluster_columns = TRUE,
              cluster_rows = TRUE,
              row_title = "免疫细胞类型",
              column_title = paste("免疫浸润聚类 (k =", k_selected, ")"),
              heatmap_legend_param = list(title = "Z-score")
)

png("Immune_Infiltration_Clusters_Heatmap.png", width = 14, height = 10, units = "in", res = 300)
draw(ht)
dev.off()







#### hub gene 在肺功能中的作用 GSE32537 ################
datadir='/Users/wuhongbo/Desktop/wjx/data_cleaning/output'
setwd(datadir)
expr_matrix=readRDS('GSE32537-array.rds')
grodir='/Users/wuhongbo/Desktop/wjx/data_cleaning/input/group'
setwd(grodir)
GROUP=readRDS('GSE32537_anno.rds')
expr_matrix=expr_matrix[,rownames(GROUP)]



library(ggplot2)
library(ggpubr)
library(tidyr)
library(dplyr)
library(gridExtra)
hub_genes <- c("COMP","CTHRC1", "CXCL14", "THY1", "COL1A1", 
               "POSTN","IGF1", "COL1A2", "CXCL12", "ADAM12","TGFB3")

# 1. 提取hub基因的表达数据
hub_expr <- expr_matrix[rownames(expr_matrix) %in% hub_genes, ]
hub_expr <- as.data.frame(t(hub_expr))  # 转置为样本×基因
hub_expr$geo_accession <- rownames(hub_expr)

# 2. 准备临床数据
GROUP$DLCO <- as.numeric(GROUP$DLCO)
GROUP$FVC <- as.numeric(GROUP$FVC)

# 3. 合并数据
merged_data <- merge(hub_expr, GROUP, by = "geo_accession")

# 4. 分析每个基因与肺功能指标的相关性
cor_results <- data.frame()

for(gene in hub_genes) {
  # DLCO相关性
  cor_dlco <- cor.test(merged_data[[gene]], merged_data$DLCO, method = "pearson")
  # FVC相关性
  cor_fvc <- cor.test(merged_data[[gene]], merged_data$FVC, method = "pearson")
  
  cor_results <- rbind(cor_results,
                       data.frame(Gene = gene,
                                  Metric = "DLCO",
                                  r = cor_dlco$estimate,
                                  p = cor_dlco$p.value),
                       data.frame(Gene = gene,
                                  Metric = "FVC",
                                  r = cor_fvc$estimate,
                                  p = cor_fvc$p.value))
}

# 格式化p值
cor_results$p_format <- ifelse(cor_results$p < 0.001, "p < 0.001",
                               ifelse(cor_results$p < 0.01, "p < 0.01",
                                      ifelse(cor_results$p < 0.05, "p < 0.05",
                                             sprintf("p = %.3f", cor_results$p))))

print(cor_results)

# 5. 创建可视化图形
plots <- list()

# 为每个基因创建两个散点图（DLCO和FVC）
for(i in seq_along(hub_genes)) {
  gene <- hub_genes[i]
  
  # DLCO散点图
  p_dlco <- ggplot(merged_data, aes_string(x = gene, y = "DLCO")) +
    geom_point(aes(color = Group), size = 3, alpha = 0.7) +
    geom_smooth(method = "lm", se = TRUE, color = "blue", fill = "lightblue") +
    stat_cor(method = "pearson", label.x.npc = "left", label.y.npc = "top") +
    labs(title = paste(gene, "vs DLCO"),
         x = paste(gene, "Expression (log2)"),
         y = "DLCO (%)") +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  # FVC散点图
  p_fvc <- ggplot(merged_data, aes_string(x = gene, y = "FVC")) +
    geom_point(aes(color = Group), size = 3, alpha = 0.7) +
    geom_smooth(method = "lm", se = TRUE, color = "red", fill = "pink") +
    stat_cor(method = "pearson", label.x.npc = "left", label.y.npc = "top") +
    labs(title = paste(gene, "vs FVC"),
         x = paste(gene, "Expression (log2)"),
         y = "FVC (%)") +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  plots[[2*i-1]] <- p_dlco
  plots[[2*i]] <- p_fvc
}

# 6. 合并所有图形
grid.arrange(grobs = plots, ncol = 4, 
             top = "Hub Genes Expression vs Lung Function Parameters")

# 7. 可选：创建相关性热图
cor_matrix <- matrix(NA, nrow = length(hub_genes), ncol = 2,
                     dimnames = list(hub_genes, c("DLCO", "FVC")))

for(gene in hub_genes) {
  cor_matrix[gene, "DLCO"] <- cor(merged_data[[gene]], merged_data$DLCO, use = "complete.obs")
  cor_matrix[gene, "FVC"] <- cor(merged_data[[gene]], merged_data$FVC, use = "complete.obs")
}

# 转换为长格式用于ggplot
cor_long <- as.data.frame(cor_matrix) %>%
  tibble::rownames_to_column("Gene") %>%
  pivot_longer(cols = c("DLCO", "FVC"), 
               names_to = "Parameter", 
               values_to = "Correlation")

# 热图
heatmap_plot <- ggplot(cor_long, aes(x = Parameter, y = Gene, fill = Correlation)) +
  geom_tile(color = "white", size = 1) +
  geom_text(aes(label = sprintf("%.3f", Correlation)), color = "white", size = 5) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                       midpoint = 0, limits = c(-1, 1)) +
  labs(title = "Correlation Heatmap: Hub Genes vs Lung Function",
       x = "Lung Function Parameter",
       y = "Gene",
       fill = "Pearson r") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
        plot.title = element_text(hjust = 0.5))

print(heatmap_plot)

# 8. 可选：分组比较（IPF vs Control）
group_plots <- list()

for(gene in hub_genes) {
  p_group <- ggplot(merged_data, aes_string(x = "Group", y = gene, fill = "Group")) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    geom_jitter(width = 0.2, size = 1, alpha = 0.5) +
    stat_compare_means(method = "t.test", label = "p.format") +
    labs(title = paste(gene, "Expression by Group"),
         x = "Group",
         y = paste(gene, "Expression (log2)")) +
    theme_minimal() +
    theme(legend.position = "none")
  
  group_plots[[gene]] <- p_group
}

# 显示分组比较图
grid.arrange(grobs = group_plots, ncol = 3,
             top = "Hub Genes Expression in IPF vs Control Groups")

# 9. 保存结果
# 保存相关性结果
write.csv(cor_results, "hub_genes_lung_function_correlation.csv", row.names = FALSE)

# 保存合并后的数据
write.csv(merged_data, "merged_expression_clinical_data.csv", row.names = FALSE)

cat("分析完成！\n")
cat("发现：\n")
for(i in 1:nrow(cor_results)) {
  cat(sprintf("%s与%s的相关性: r = %.3f, %s\n",
              cor_results$Gene[i], cor_results$Metric[i],
              cor_results$r[i], cor_results$p_format[i]))
}

######### hub gene 分组差异通路GSEA ######################

datadir='/Users/wuhongbo/Desktop/wjx/data_quality/output'
setwd(datadir)
expr_matrix=readRDS('Batch_limma_data.rds')
GROUP=readRDS('combined_group.rds')
setwd('/Users/wuhongbo/Desktop/wjx/data_quality/JTM_4_ NSCLC')


library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(dplyr)
library(tidyr)

library(patchwork)
library(cowplot)
# Publication Theme
theme_pub <- theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(
      size = 16,
      face = "bold",
      hjust = 0.5,
      colour = "black"
    ),
    axis.title = element_text(
      size = 14,
      face = "bold"
    ),
    axis.text = element_text(
      size = 12,
      colour = "black"
    ),
    legend.title = element_text(
      size = 12,
      face = "bold"
    ),
    legend.text = element_text(
      size = 11
    ),
    
    panel.grid.major = element_line(
      colour = "grey90",
      linewidth = 0.3
    ),
    panel.grid.minor = element_blank(),
    strip.text = element_text(
      size = 13,
      face = "bold"
    ),
    plot.margin = margin(10, 10, 10, 10)
  )

# Hub基因列表
hub_genes <- c("CTHRC1", "COMP", "CXCL14")
# 存储最终过滤结果的列表
all_gsea_results <- list()

# 主循环：对每个 Hub 基因执行分组 GSEA
for (gene in hub_genes) {
  cat("\n========== 正在分析基因:", gene, "==========\n")
  
  # --- 1. 根据该基因表达量中位数分组 ---
  gene_expr <- expr_matrix[gene, ]
  median_val <- median(gene_expr, na.rm = TRUE)
  group_labels <- ifelse(gene_expr >= median_val, "High", "Low")
  names(group_labels) <- colnames(expr_matrix)
  
  high_samples <- names(group_labels)[group_labels == "High"]
  low_samples  <- names(group_labels)[group_labels == "Low"]
  
  cat(paste0("  分组完成: High (n=", length(high_samples), "), Low (n=", length(low_samples), ")\n"))
  
  # --- 2. 计算每个基因的 t 统计量（排序指标）---
  gene_stats <- apply(expr_matrix, 1, function(x) {
    x_high <- x[high_samples]
    x_low  <- x[low_samples]
    # 避免单样本组导致 t.test 报错
    if (length(x_high) < 2 || length(x_low) < 2) return(0)
    t_res <- t.test(x_high, x_low, var.equal = FALSE)
    return(t_res$statistic)
  })
  
  ranked_list <- sort(gene_stats, decreasing = TRUE)
  
  # --- 3. 基因 ID 转换 (SYMBOL → ENTREZID) ---
  symbols <- names(ranked_list)
  entrez_map <- bitr(symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db, drop = TRUE)
  
  if (nrow(entrez_map) == 0) {
    cat("  警告: 无基因成功转换为 Entrez ID，跳过 GSEA。\n")
    all_gsea_results[[gene]] <- data.frame()
    next
  }
  
  entrez_ranked <- ranked_list[entrez_map$SYMBOL]
  names(entrez_ranked) <- entrez_map$ENTREZID
  entrez_ranked <- na.omit(entrez_ranked)
  
  cat(paste0("  成功转换 ", length(entrez_ranked), " 个基因用于 GSEA。\n"))
  
  # --- 4. 执行 GSEA (KEGG) ---
  gsea_res <- gseKEGG(
    geneList = entrez_ranked,
    organism = "hsa",
    keyType = "ncbi-geneid",
    nPerm = 10000,
    minGSSize = 10,
    maxGSSize = 500,
    pvalueCutoff = 1,      # 后续统一过滤
    pAdjustMethod = "BH",
    verbose = FALSE,
    seed = 123
  )
  
  # --- 5. 过滤显著通路（按您的标准）---
  if (is.null(gsea_res) || nrow(gsea_res) == 0) {
    cat("  GSEA 未返回有效结果。\n")
    all_gsea_results[[gene]] <- data.frame()
    next
  }
  
  gsea_df <- as.data.frame(gsea_res)
  filtered_df <- gsea_df %>%
    filter(pvalue < 0.05, abs(NES) > 1.0, p.adjust < 0.25) %>%
    arrange(desc(abs(NES)))
  
  cat(paste0("  发现 ", nrow(filtered_df), " 个显著富集通路。\n"))
  
  # 保存过滤结果
  all_gsea_results[[gene]] <- filtered_df
  
  # --- 6. 可视化（仅当有显著通路时）---
  if (nrow(filtered_df) > 0) {
    top_n <- min(10, nrow(filtered_df))
    top_ids <- head(filtered_df$ID, top_n)
    
    # 6.1 GSEA 富集图（每个通路一张）
    # 6.1 Publication-quality GSEA plot
    for (pid in top_ids) {
      
      pathway_name <- as.character(gsea_res[pid, "Description"])
      
      # 文件名安全化
      pathway_short <- substr(gsub("[^A-Za-z0-9]", "", pathway_name), 1, 3)
      
      # 生成基础图
      p <- gseaplot2(
        gsea_res,
        geneSetID = pid,
        
        title = pathway_name,
        
        base_size = 16,
        
        rel_heights = c(1.5, 0.5, 1),
        
        subplots = 1:3,
        
        pvalue_table = TRUE,
        
        ES_geom = "line",
        
        color = "#00AF54"
      )
      
      # IMPORTANT:
      # gseaplot2返回的是patchwork/list对象
      # 不能直接 + theme()
      
      p[[1]] <- p[[1]] +
        
        labs(
          subtitle = paste0(gene, " High vs Low")
        ) +
        
        theme_bw(base_family = "Arial") +
        
        theme(
          
          plot.title = element_text(
            size = 18,
            face = "bold",
            hjust = 0.5,
            colour = "black"
          ),
          
          plot.subtitle = element_text(
            size = 13,
            hjust = 0.5,
            face = "italic",
            colour = "grey30"
          ),
          
          axis.title = element_text(
            size = 14,
            face = "bold"
          ),
          
          axis.text = element_text(
            size = 12,
            colour = "black"
          ),
          
          panel.grid.major = element_line(
            colour = "grey90",
            linewidth = 0.3
          ),
          
          panel.grid.minor = element_blank(),
          
          panel.border = element_rect(
            colour = "black",
            linewidth = 0.8
          ),
          
          legend.position = "none",
          
          plot.margin = margin(10, 10, 10, 10)
        )
      
      # Mac下推荐 cairo_pdf
      ggsave(
        filename = paste0(
          "GSEA_Group_",
          gene,
          "_",
          pathway_short,
          ".pdf"
        ),
        
        plot = p,
        
        device = cairo_pdf,
        
        width = 8,
        height = 6,
        
        dpi = 600,
        
        bg = "white"
      )
      
      # PNG
      ggsave(
        filename = paste0(
          "GSEA_Group_",
          gene,
          "_",
          pathway_short,
          ".png"
        ),
        
        plot = p,
        
        width = 8,
        height = 6,
        
        dpi = 600,
        
        bg = "white"
      )
    }
    
    # 6.2 点图（Dot plot）
    # 6.2 Publication-quality Dotplot
    if (nrow(filtered_df) > 1) {
      
      dotp <- dotplot(
        gsea_res,
        showCategory = top_n,
        x = "NES",
        color = "p.adjust",
        size = "setSize",
        font.size = 14,
        title = paste0(gene, ": KEGG pathway enrichment")
      ) +
        
        scale_color_gradient(
          low = "#D73027",
          high = "#4575B4",
          name = "Adjusted p"
        ) +
        
        theme_pub +
        
        theme(
          axis.text.y = element_text(
            size = 12,
            face = "bold"
          ),
          
          axis.text.x = element_text(
            size = 11,
            angle = 45,
            hjust = 1
          ),
          
          legend.position = "right"
        )
      
      ggsave(
        paste0("GSEA_DotPlot_", gene, ".pdf"),
        dotp,
        width = 9,
        height = 6,
        dpi = 600,
        bg = "white"
      )
      
      ggsave(
        paste0("GSEA_DotPlot_", gene, ".png"),
        dotp,
        width = 9,
        height = 6,
        dpi = 600,
        bg = "white"
      )
    }
  }
  
  cat("--------------------------------------------------\n")
}
