


suppressMessages({
  library(dplyr)
  library(Seurat)
  library(patchwork)
  library(DoubletFinder)
  library(optparse)
  library(ggplot2)
  library(S4Vectors)
})


#-------------------------------------------------------------------------------
# 参数传递
#-------------------------------------------------------------------------------

option_list <- list(
  make_option("--RDSfile", type="character", default=NULL, help="scdata.rds"),
  make_option("--resolution", type="double", default=NULL, help="根据上一步选出来的最佳分辨率"),
  make_option("--OutPath", type = "character", default="./", help = "输出文件的路径")
)

args <- parse_args(OptionParser(option_list=option_list))

SampleFile <- args$SampleFile
resolution <- args$resolution
OutPath <- args$OutPath



## 选择合适的分辨率后需要显式指定分辨率

#-------------------------------------------------------------------------------
# 整合后降维聚类
#-------------------------------------------------------------------------------

resolution <- paste0("RNA_snn_res.", as.character(resolution))
Idents(scdata) <- resolution
scdata <- RunUMAP(scdata, reduction = "harmony", dims = 1:30)


#-------------------------------------------------------------------------------
# 可视化
#-------------------------------------------------------------------------------

# Visualization
p <- DimPlot(scdata, reduction = "umap", label = TRUE)
ggplot2::ggsave(file = "umap.cluster.pdf", plot = p , width = 6, height = 5)

p.splitbysample <- DimPlot(scdata, reduction = "umap", label = TRUE, split.by = "orig.ident")  # split.by 指定的meta.data中的列
ggplot2::ggsave(file = add_path(OutPath,"umap.cluster.splitbysample.pdf"), plot = p.splitbysample,  width = 10, height = 5)


# 每个样本细胞数量统计
# cluster_count <- scdata@meta.data %>%
#   group_by(seurat_clusters, orig.ident) %>%
#   summarise(n = n(), .groups = "drop")

# write.table(x = cluster_count, file = "cluster_cells_count.tsv", row.names = F, quote = F, sep = "\t")


#-------------------------------------------------------------------------------
# 提取各 cluster 的 marker 基因
#-------------------------------------------------------------------------------
message("\n\nFinding markers for all clusters\n\n")

all_markers <- FindAllMarkers(scdata, 
                              only.pos = TRUE,
                              min.pct = 0.25,
                              logfc.threshold = 0.25,
                              verbose = FALSE) |>
  dplyr::select(gene, cluster, everything())

write.table(all_markers, add_path(OutPath,"all_markers.tsv"), row.names = FALSE, quote = F, sep = "\t")

# 每个 cluster 取 top10 marker
# top10 <- all_markers %>%
#   group_by(cluster) %>%
#   slice_max(order_by = avg_log2FC, n = 10) %>%
#   ungroup() |>
#   dplyr::select(gene, cluster, everything())


# 核糖体蛋白基因：RPL13A、RPS27A、RPLP0、MRPL12 等；写成精确匹配，避免误伤 RPS6KA1 这类激酶
RIBO_PATTERN <- "^(RPL[0-9]+[A-Z]{0,2}|RPS[0-9]+[A-Z]{0,2}|RPLP[0-9]|MRPL[0-9]+|MRPS[0-9]+[A-Z]?)$"
 
 
select_top_markers <- function(markers, n = 10,
                               padj_max        = 0.05,
                               lfc_min         = 0.5,    # 最低 log2FC
                               pct1_min        = 0.25,   # 本 cluster 中至少多少比例细胞表达
                               diff_min        = 0.10,   # pct.1 - pct.2 的最低差值
                               exclude_pattern = RIBO_PATTERN,
                               exclude_genes   = NULL,   # 例如线粒体基因、血红蛋白基因
                               unique_genes    = FALSE,
                               warn_excluded   = 0.3) {  # 显著基因中被排除基因占比超过此值时提示
 
  is_excluded <- function(g) grepl(exclude_pattern, g) | g %in% exclude_genes
 
  # 所有显著上调基因（用于诊断）
  sig <- markers |>
    dplyr::filter(p_val_adj < padj_max, avg_log2FC >= lfc_min) |>
    dplyr::mutate(pct.diff = pct.1 - pct.2,
                  score    = avg_log2FC * pct.diff,
                  excluded = is_excluded(gene))
 
  # 通过特异性过滤、且不在排除名单中的基因
  kept <- sig |>
    dplyr::filter(!excluded, pct.1 >= pct1_min, pct.diff >= diff_min)
 
  # 一个基因只保留在得分最高的 cluster
  if (unique_genes) {
    kept <- kept |>
      dplyr::group_by(gene) |>
      dplyr::slice_max(score, n = 1, with_ties = FALSE) |>
      dplyr::ungroup()
  }
 
  top <- kept |>
    dplyr::group_by(cluster) |>
    dplyr::slice_max(score, n = n, with_ties = FALSE) |>
    dplyr::ungroup() |>
    dplyr::arrange(cluster, dplyr::desc(score)) |>
    dplyr::select(gene, cluster, score, avg_log2FC, pct.1, pct.2, pct.diff, p_val_adj)
 
  # 诊断表：每个 cluster 的显著基因数、被排除基因占比、最终入选数
  cl_levels <- levels(factor(markers$cluster))
  diag <- data.frame(cluster = cl_levels) |>
    dplyr::left_join(
      sig |>
        dplyr::group_by(cluster = as.character(cluster)) |>
        dplyr::summarise(n_sig = dplyr::n(),
                         frac_excluded = round(mean(excluded), 3),
                         .groups = "drop"),
      by = "cluster") |>
    dplyr::left_join(
      top |> dplyr::count(cluster = as.character(cluster), name = "n_selected"),
      by = "cluster") |>
    dplyr::mutate(n_sig      = dplyr::coalesce(n_sig, 0L),
                  n_selected = dplyr::coalesce(n_selected, 0L))
 
  few <- diag$cluster[diag$n_selected < n]
  if (length(few) > 0) {
    message("通过筛选的 marker 少于 ", n, " 个的 cluster：", paste(few, collapse = ", "),
            "\n  可能是过度聚类拆出来的，也可能是低质量或双胞 cluster")
  }
  hi_ex <- diag$cluster[!is.na(diag$frac_excluded) & diag$frac_excluded > warn_excluded]
  if (length(hi_ex) > 0) {
    message("显著基因中核糖体/线粒体等占比 > ", warn_excluded, " 的 cluster：", paste(hi_ex, collapse = ", "),
            "\n  这类 cluster 常由细胞状态或质量驱动，注释前建议回看它的 QC 指标")
  }
 
  attr(top, "diagnosis") <- diag
  top
}
 
top10_markers <- select_top_markers(all_markers)

write.table(top10_markers, add_path(OutPath,"top10_markers.tsv"), row.names = FALSE, quote = F, sep = "\t")


# 每个 cluster 取 top5 marker（用于热图展示，太多会挤）
topn_marker <- all_markers %>%
  group_by(cluster) %>%
  # slice_max(order_by = avg_log2FC, n = 3) %>%
  top_n(n = 5, wt = avg_log2FC) %>%
  ungroup()

p_dot <- DotPlot(scdata, features = unique(topn_marker$gene)) +
  RotatedAxis() +
  theme(axis.text.x = element_text(size = 7))

ggsave(filename = add_path(OutPath,"marker.cluster.pdf"), plot = p_dot, width = 16, height = 5)




