


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
top10 <- all_markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 10) %>%
  ungroup() |>
  dplyr::select(gene, cluster, everything())

write.table(top10, add_path(OutPath,"top10_markers.tsv"), row.names = FALSE, quote = F, sep = "\t")


# 每个 cluster 取 top5 marker（用于热图展示，太多会挤）
topn_marker <- all_markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 3) %>%
  ungroup()

p_dot <- DotPlot(scdata, features = unique(topn_marker$gene)) +
  RotatedAxis() +
  theme(axis.text.x = element_text(size = 7))

ggsave(filename = add_path(OutPath,"marker.cluster.pdf"), plot = p_dot, width = 16, height = 5)




