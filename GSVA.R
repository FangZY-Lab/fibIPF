
library(GSVA)      
library(ggplot2)   
fibrosis_geneset <- list(FibrosisSignature = frgs)
param <- ssgseaParam(expr_matrix, fibrosis_geneset)
fibrosis_scores <- gsva(param, verbose=TRUE)

fibrosis_score <- as.numeric(fibrosis_scores)
sample_info <- data.frame(
  SampleID = colnames(expr_matrix),
  FibrosisScore = fibrosis_score,
  Group = GROUP$Group,  
  ClinicalOutcome = ifelse(GROUP$Group=='IPF',1,0),# 0=健康，1=纤维化
  Project=GROUP$Batch
)
is_ipf <- sample_info$Group == "IPF"
sample_info <- sample_info[is_ipf, ]
expr_matrix <- expr_matrix[, sample_info$SampleID]
sample_info$FibrosisGroup <- ifelse(
  sample_info$FibrosisScore > median(sample_info$FibrosisScore), "High", "Low")


##### GSVA HALLMARK #######
library(GSVA)
library(GSEABase)
library(limma)
library(pheatmap)
library(msigdbr)
library(dplyr)
# 获取 HALLMARK 基因集
hallmark_df <- msigdbr(species = "Homo sapiens", category = "H") %>%
  dplyr::select(gs_name, gene_symbol)
hallmark_list <- split(hallmark_df$gene_symbol, hallmark_df$gs_name)
hallmark_list <- lapply(hallmark_list, function(genes) {
  genes[genes %in% rownames(expr_matrix)]
})
hallmark_list <- hallmark_list[sapply(hallmark_list, length) > 0]
# GSVA 分析
param_hallmark <- ssgseaParam(expr_matrix, hallmark_list)
gsva_score <- gsva(param_hallmark, verbose = FALSE)
# 差异分析
group <- sample_info$FibrosisGroup # "High" / "Low"
group <- factor(group, levels = c("Low", "High")) # High - Low
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)
fit <- lmFit(gsva_score, design)
contrast <- makeContrasts(High - Low, levels = design)
fit2 <- contrasts.fit(fit, contrast)
fit2 <- eBayes(fit2)
diff_result <- topTable(fit2, number = Inf, sort.by = "P")
sig_pathways <- diff_result[diff_result$adj.P.Val < 0.01, ]
sig_pathways <- sig_pathways[order(-abs(sig_pathways$logFC)), ]
top20_pathways <- head(rownames(sig_pathways), 20)
# 数据与图例清理
top_pathways <- head(rownames(sig_pathways), 10)
mat <- gsva_score[top_pathways, ]
# 样本排序：High 在左，Low 在右
ord <- order(sample_info$FibrosisGroup == "Low")
mat <- mat[, ord]
# 自定义短通路名
pathway_short <- c(
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION" = "EMT",
  "HALLMARK_ANGIOGENESIS" = "Angiogenesis",
  "HALLMARK_KRAS_SIGNALING_UP" = "KRAS Signaling Up",
  "HALLMARK_COAGULATION" = "Coagulation",
  "HALLMARK_APICAL_JUNCTION" = "Apical Junction",
  "HALLMARK_MYOGENESIS" = "Myogenesis",
  "HALLMARK_HEDGEHOG_SIGNALING" = "Hedgehog Signaling",
  "HALLMARK_APICAL_SURFACE" = "Apical Surface",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB" = "TNFα/NF-κB",
  "HALLMARK_WNT_BETA_CATENIN_SIGNALING" = "Wnt/β-catenin",
  "HALLMARK_IL6_JAK_STAT3_SIGNALING" = "IL6–JAK–STAT3",
  "HALLMARK_INFLAMMATORY_RESPONSE" = "Inflammatory Response",
  "HALLMARK_SPERMATOGENESIS" = "Spermatogenesis"
)
pathway_names=rownames(mat)
display_names <- pathway_short[pathway_names]
rownames(mat) <- display_names
annotation_col <- data.frame(
  "Fibrosis Level" = sample_info$FibrosisGroup[ord],
  row.names = colnames(mat),
  check.names = FALSE
)
annotation_colors <- list(
  "Fibrosis Level" = c(
    "Low" = "#93cdd9",
    "High" = "#fa6647"
  )
)
# 绘制并保存热图
ph=pheatmap(
  mat,
  annotation_col = annotation_col,
  annotation_colors = annotation_colors,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  scale = "row",
  show_rownames = TRUE,
  show_colnames = FALSE,
  color = colorRampPalette(c("#2874A6", "white", "red"))(100),
  fontsize_row = 9,
  border_color = NA,
  treeheight_row = 30
)
heatmap_order <- rownames(mat)[ph$tree_row$order]
cairo_pdf("hallmark_heatmap.pdf", width = 6.2, height = 3.5, onefile = TRUE)
print(ph)
dev.off()

##### 4.4.2 HAllmark 箱+辐射 #########
library(reshape2)
library(ggpubr)
library(ggplot2)
# 提取 top hallmark score
plot_data <- as.data.frame(t(gsva_score[pathway_names, ]))
plot_data$Sample <- rownames(plot_data)
plot_data$FibrosisGroup <- sample_info[
  match(plot_data$Sample, sample_info$SampleID),
  "FibrosisGroup"
]
# 转长表
plot_data_long <- melt(
  plot_data,
  id.vars = c("Sample", "FibrosisGroup"),
  variable.name = "Pathway",
  value.name = "GSVA_Score"
)

plot_data_long$FibrosisGroup <- factor(
  plot_data_long$FibrosisGroup,
  levels = c("Low","High")
)
plot_data_long$Pathway <- factor(
  pathway_short[as.character(plot_data_long$Pathway)],
  levels = rev(heatmap_order)   # coord_flip 后需要反转
)

# y轴文字颜色也按箱线图显示顺序匹配
logFC_vals <- sig_pathways[names(pathway_short)[match(heatmap_order, pathway_short)], "logFC"]
pathway_colors <- ifelse(rev(logFC_vals) > 0,"#C0392B","#2980B9")
# 作图
p <- ggplot(
  plot_data_long,
  aes(
    x = Pathway,
    y = GSVA_Score,
    fill = FibrosisGroup
  )
) +
  geom_boxplot(
    width = 0.7,
    outlier.shape = NA,
    linewidth = 0.5
  ) +
  geom_jitter(
    aes(color = FibrosisGroup),
    position = position_jitterdodge(
      jitter.width = 0.15,
      dodge.width = 0.75
    ),
    size = 0.6,
    alpha = 0.35,
    show.legend = FALSE
  ) +
  coord_flip() +
  scale_fill_manual(breaks = c("Low","High"),
    values = c(
      "Low" = "#93cdd9",
      "High" = "#fa6647"
    )
  ) +
  scale_color_manual(breaks = c("Low","High"),
    values = c(
      "Low" = "#93cdd9",
      "High" = "#fa6647"
    )
  ) +
  stat_compare_means(
    aes(group = FibrosisGroup),
    method = "wilcox.test",
    label = "p.signif",
    hide.ns = TRUE,
    p.adjust.method = "fdr",
    size = 4
  ) +
  labs(
    x = NULL,
    y = "ssGSEA enrichment score",
    fill = "Fibrosis Level"
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.y = element_text(
      color = pathway_colors    ),
    legend.position = "top"
  )

p

ggsave(
  "hallmark_boxplot.pdf",
  p,
  width = 5,
  height = 4,
  device = cairo_pdf
)


###辐射图
library(Hmisc)
library(linkET)
library(Cairo)
library(vegan)

### Hallmark Mantel test corrplot

hallmark_df <- as.data.frame(t(gsva_score[top_pathways, , drop = FALSE]))
hallmark_df <- hallmark_df[, colSums(is.na(hallmark_df)) == 0, drop = FALSE]
hallmark_mat <- as.matrix(hallmark_df)

# 通路名称转换，避免 pathway_short 不全时出现 NA
new_names <- pathway_short[colnames(hallmark_mat)]
new_names[is.na(new_names)] <- gsub(
  "_", " ",
  gsub("^HALLMARK_", "", colnames(hallmark_mat)[is.na(new_names)])
)
colnames(hallmark_mat) <- new_names

## 1. 指定与热图一致的通路顺序
hallmark_order <- c("Angiogenesis",
  "EMT",
  "Coagulation",
  "KRAS Signaling Up",
  "TNFα/NF-κB","Apical Surface","Myogenesis","Wnt/β-catenin",
  "Apical Junction",
  
  "Hedgehog Signaling"
  
  
  
)

hallmark_order <- hallmark_order[hallmark_order %in% colnames(hallmark_mat)]
hallmark_mat <- hallmark_mat[, hallmark_order, drop = FALSE]
pathway_names <- colnames(hallmark_mat)

# Fibrosis group distance matrix: Low = 0, High = 1
fibrosis_vec <- as.numeric(sample_info$FibrosisGroup == "High")
fibrosis_dist <- dist(fibrosis_vec, method = "euclidean")

mantel_df <- data.frame(
  var1 = rep("Fibrosis Program", length(pathway_names)),
  var2 = pathway_names,
  r_value = NA_real_,
  p_value = NA_real_,
  stringsAsFactors = FALSE
)

set.seed(123)

for (i in seq_along(pathway_names)) {
  y <- hallmark_mat[, i]
  pathway_dist <- dist(y, method = "euclidean")
  
  mt <- mantel(
    fibrosis_dist,
    pathway_dist,
    method = "spearman",
    permutations = 9999
  )
  
  mantel_df$r_value[i] <- mt$statistic
  mantel_df$p_value[i] <- mt$signif
}

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

mantel_sig <- subset(mantel_df, p_value < 0.05)

cairo_pdf(
  "hallmark_corrplot.pdf",
  width = 7,
  height = 6
)

p_mantel <- qcorrplot(
  rcorr(hallmark_mat, type = "spearman"),
  type = "upper",
  diag = FALSE,
  grid_size = 0.4,
  grid_col = "lightgray"
) +
  geom_square(linetype = 0) +
  geom_couple(
    aes(
      colour = p_group,
      size = r_group
    ),
    data = mantel_sig,
    curvature = 0.1
  ) +
  set_corrplot_style(
    colours = c("#FF8040", "white", "#5BC2CD")
  ) +
  scale_colour_manual(
    name = "Mantel's p",
    values = c(
      "<0.001" = "#87CEEB",
      "0.001–0.05" = "#9ACD32"
    ),
    breaks = c("<0.001", "0.001–0.05"),
    drop = FALSE,
    guide = guide_legend(order = 1)
  ) +
  scale_size_manual(
    name = "Mantel's |r|",
    values = c(
      "<0.1" = 0.2,
      "0.1–0.2" = 0.5,
      ">0.2" = 0.9
    ),
    breaks = c("<0.1", "0.1–0.2", ">0.2"),
    drop = FALSE,
    guide = guide_legend(order = 2)
  ) +
  guides(
    fill = guide_colorbar(
      title = "r",
      order = 3
    )
  ) +
  
  theme(
    axis.text.x.top = element_text(
      angle = 45,
      hjust = 0,
      vjust = 0
    ),
    legend.box = "vertical"
  )

print(p_mantel)
dev.off()

