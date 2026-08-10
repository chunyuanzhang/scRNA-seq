################################################################################
# 单细胞分析流程
# 1、先按照线粒体比例把最差的细胞过滤掉
# 2、双胞检测
# 3、整合样本，降维
################################################################################


# setwd("~/Desktop/04.湘湖实验室/姜雨鸡单细胞")
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
  make_option("--SampleFile", type="character", default=NULL, help="config/samples.csv"),
  make_option("--MTpattern", type="character", default=NULL, help="线粒体基因前缀"),
  make_option("--percentMT", type="double", default=NULL, help="基础指控要求线粒体比例")
)

args <- parse_args(OptionParser(option_list=option_list))

SampleFile <- args$SampleFile
MTpattern <- args$MTpattern
percentMT <- args$percentMT

# SampleFile <- "~/Desktop/04.湘湖实验室/姜雨鸡单细胞/sampleandpopulation.csv"
# MTpattern <- "J6367"
# percentMT <- 15


plus_one <- function(r){
  r <<- r + 1
}

r <- 0
#-------------------------------------------------------------------------------
# 读取数据
#-------------------------------------------------------------------------------

sampletable <- read.delim(file = SampleFile, sep = ",", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, comment.char = "#")
samples <- sampletable$SampleID
#dirlist <- paste0("~/Desktop/04.湘湖实验室/姜雨鸡单细胞/", samples) 
dirlist <- paste0("result/02.Count/", samples, "/outs/filter_matrix/") 

message("\n\n", plus_one(r), ": Reading samples: ", paste0(dirlist, collapse = ", "), "\n\n")
names(dirlist) <- samples
scdata.data <- Read10X(data.dir = dirlist)
scdata <- CreateSeuratObject(counts = scdata.data, project = "jiangyu", min.cells = 3, min.features = 200)



#-------------------------------------------------------------------------------
# 基础质量控制
#-------------------------------------------------------------------------------

message("\n\n", plus_one(r) ,": Quality control\n\n")
# 计算线粒体比例
scdata[["percent.mt"]] <- PercentageFeatureSet(scdata, pattern = MTpattern)

# Visualize QC metrics as a violin plot
# VlnPlot(scdata, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# 质量控制
upperlimit <- scdata@meta.data |> 
  group_by(orig.ident) |> 
  summarise(upperlimit_nFeature_RNA = nFeature_RNA |> quantile(0.95)) |> 
  ungroup() |> 
  as.data.frame()
scdata@meta.data$upperlimit_nFeature_RNA <- scdata@meta.data |> 
  left_join(upperlimit, by = "orig.ident") |> 
  pull(upperlimit_nFeature_RNA)
scdata <- subset(scdata, cells = scdata@meta.data |> filter(nFeature_RNA > 500 & nFeature_RNA < upperlimit_nFeature_RNA & percent.mt < percentMT) |> rownames())


#-------------------------------------------------------------------------------
# 双胞检测
# DoubletFinder 需要在每个样本上独立跑，所以这里先拆再合
#   nFeature_RNA 太高，或者 percent.mt 太高都有可能是双胞
#-------------------------------------------------------------------------------

message("\n\n",plus_one(r),": Double cell check\n\n")

seu_list <- SplitObject(scdata, split.by = "orig.ident")

## 对细胞数量进行计数
sample_cellnumber <- lapply(seu_list, function(x) {
  ncol(x)
}) |> as.data.frame() |> t() |> as.data.frame() |> setNames("NumberOfCells")

write.table(x = sample_cellnumber, file = "NumberOfCells.tsv", quote = F, sep = "\t")
if(any(sample_cellnumber$NumberOfCells < 100 )){
  message("\n\n\n存在细胞数量过低样本，请核查\n\n\n")
  quit()
}


seu_list <- lapply(seu_list, function(x) {
  SampleID <- x@meta.data$orig.ident |> unique() |> unfactor()
  message("\n\n当前样本 ",SampleID, "\n\n")
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
  x <- ScaleData(x)
  
  max_pcs <- min(ncol(x), nrow(x), 51) - 1
  use_dims <- min(max_pcs, 10)
  
  x <- RunPCA(x, npcs = max_pcs)
  x <- RunUMAP(x, dims = 1:use_dims)
  x <- doubletFinder(x, PCs = 1:10, pN = 0.25, pK = 0.09,
                     nExp = ncol(x) * 0.075, reuse.pANN = NULL, sct = FALSE)
  df_col <- grep("^DF", colnames(x@meta.data), value = TRUE)
  x <- subset(x, cells = rownames(x@meta.data[x@meta.data[[df_col]] == "Singlet", ]))
  x
})

message("\n\nMerge samples\n\n")
if (length(seu_list) == 1) {
  scdata <- seu_list[[1]]
} else {
  scdata <- merge(seu_list[[1]], y = seu_list[-1], project = "jiangyu")
}
remove(seu_list)


#-------------------------------------------------------------------------------
# 整合样本
# Seurat v5 标准流程：归一化 → HVG → Scale → PCA
# v5 中 merge 后 RNA assay 已按 orig.ident 自动分 layer
#-------------------------------------------------------------------------------

message("\n\n",plus_one(r),": Normalization & PCA\n\n")
scdata <- NormalizeData(scdata)
scdata <- FindVariableFeatures(scdata, selection.method = "vst", nfeatures = 2000)
scdata <- ScaleData(scdata)
scdata <- RunPCA(scdata)

# saveRDS(object = scdata, file = "~/Desktop/04.湘湖实验室/姜雨鸡单细胞/scdata.filterbydoublecell.rds")


#-------------------------------------------------------------------------------
# CCA 整合（Seurat v5 IntegrateLayers）
#-------------------------------------------------------------------------------

message("\n\n",plus_one(r),": CCA Integration\n\n")
# scdata <- IntegrateLayers(
#   object         = scdata,
#   method         = CCAIntegration,
#   orig.reduction = "pca",
#   new.reduction  = "integrated.cca",
#   verbose        = FALSE
# )

scdata <- IntegrateLayers(
  object         = scdata,
  method         = HarmonyIntegration,
  orig.reduction = "pca",
  new.reduction  = "harmony",
  verbose        = FALSE
)

scdata[["RNA"]] <- JoinLayers(scdata[["RNA"]])

#-------------------------------------------------------------------------------
# 整合后降维聚类
#-------------------------------------------------------------------------------

message("\n\n", plus_one(r),": Clustering & UMAP\n\n")
scdata <- FindNeighbors(scdata, reduction = "harmony", dims = 1:30)
scdata <- FindClusters(scdata, resolution = 1)
scdata <- RunUMAP(scdata, reduction = "harmony", dims = 1:30)

saveRDS(object = scdata, file = "scdata.rds")


#-------------------------------------------------------------------------------
# 可视化
#-------------------------------------------------------------------------------

# Visualization
p <- DimPlot(scdata, reduction = "umap", label = TRUE)
ggplot2::ggsave(file = "umap.cluster.pdf", plot = p , width = 6, height = 5)

p.splitbysample <- DimPlot(scdata, reduction = "umap", label = TRUE, split.by = "orig.ident")  # split.by 指定的meta.data中的列
ggplot2::ggsave(file = "umap.cluster.splitbysample.pdf", plot = p.splitbysample,  width = 10, height = 5)


# 每个样本细胞数量统计
cluster_count <- scdata@meta.data %>%
  group_by(seurat_clusters, orig.ident) %>%
  summarise(n = n(), .groups = "drop")

write.table(x = cluster_count, file = "cluster_cells_count.tsv", row.names = F, quote = F, sep = "\t")


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

write.table(all_markers, "all_markers.tsv", row.names = FALSE, quote = F, sep = "\t")

# 每个 cluster 取 top10 marker
top10 <- all_markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 10) %>%
  ungroup() |>
  dplyr::select(gene, cluster, everything())

write.table(top10, "top10_markers.tsv", row.names = FALSE, quote = F, sep = "\t")


# 每个 cluster 取 top5 marker（用于热图展示，太多会挤）
topn_marker <- all_markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 3) %>%
  ungroup()

p_dot <- DotPlot(scdata, features = unique(topn_marker$gene)) +
  RotatedAxis() +
  theme(axis.text.x = element_text(size = 7))

ggsave(filename = "marker.cluster.pdf", plot = p_dot, width = 16, height = 5)




