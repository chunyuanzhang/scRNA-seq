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
})


#-------------------------------------------------------------------------------
# 参数传递
#-------------------------------------------------------------------------------

option_list <- list(
  make_option("--SampleFile", type="character", default=NULL, help="config/samples.csv")
)

args <- parse_args(OptionParser(option_list=option_list))

SampleFile <- args$SampleFile


#-------------------------------------------------------------------------------
# 读取数据
#-------------------------------------------------------------------------------

sampletable <- read.delim(file = SampleFile, sep = ",", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, comment.char = "#")
samples <- sampletable$SampleID
#dirlist <- paste0("~/Desktop/04.湘湖实验室/姜雨鸡单细胞/", samples) 
dirlist <- paste0("result/02.Count/", samples, "/outs/filter_matrix/") 

message("\n\ntReading samples: ", paste0(dirlist, collapse = ", "), "\n\n")
names(dirlist) <- samples
scdata.data <- Read10X(data.dir = dirlist)
scdata <- CreateSeuratObject(counts = scdata.data, project = "jiangyu", min.cells = 3, min.features = 200)



#-------------------------------------------------------------------------------
# 基础质量控制
#-------------------------------------------------------------------------------

message("\n\nQuality control\n\n")
# 计算线粒体比例
scdata[["percent.mt"]] <- PercentageFeatureSet(scdata, pattern = "J6367")

# Visualize QC metrics as a violin plot
# VlnPlot(scdata, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# 质量控制
upperlimit <- scdata@meta.data |> group_by(orig.ident) |> summarise(upperlimit_nFeature_RNA = nFeature_RNA |> quantile(0.95)) |> ungroup() |> as.data.frame()
scdata@meta.data$upperlimit_nFeature_RNA <- scdata@meta.data |> left_join(upperlimit, by = "orig.ident") |> pull(upperlimit_nFeature_RNA)
scdata <- subset(scdata, cells = scdata@meta.data |> filter(nFeature_RNA > 500 & nFeature_RNA < upperlimit_nFeature_RNA & percent.mt < 15) |> rownames())


#-------------------------------------------------------------------------------
# 双胞检测
#-------------------------------------------------------------------------------

message("\n\nDouble cell check\n\n")
# nFeature_RNA 太高，或者 percent.mt 太高都有可能是双胞
# 下面进行双胞检验
## Pre-process Seurat object (standard) --------------------------------------------------------------------------------------
scdata <- NormalizeData(scdata)
scdata <- FindVariableFeatures(scdata, selection.method = "vst", nfeatures = 2000)
scdata <- ScaleData(scdata)
scdata <- RunPCA(scdata)
scdata <- RunUMAP(scdata, dims = 1:10)

## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
scdata <- doubletFinder(scdata, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = ncol(scdata) * 0.075, reuse.pANN = NULL, sct = FALSE)


# 提取过滤双胞后的数据
df_col <- colnames(scdata@meta.data) |> grep(pattern = "DF", value = T)
scdata <- subset(scdata, cells = scdata@meta.data[scdata@meta.data[[df_col]] == "Singlet",] |> rownames() )

saveRDS(object = scdata, file = "scdata.filterbydoublecell.rds")


#-------------------------------------------------------------------------------
# 整合样本
#-------------------------------------------------------------------------------

message("\n\nIntegrateLayers\n\n")

scdata <- IntegrateLayers(object = scdata, 
    method = CCAIntegration, 
    orig.reduction = "pca", 
    new.reduction = "integrated.cca",
    verbose = FALSE)

# re-join layers after integration
scdata[["RNA"]] <- JoinLayers(scdata[["RNA"]])

#-------------------------------------------------------------------------------
# 整合后的样本降维聚类【整合后不能再运行normalized】
#-------------------------------------------------------------------------------

message("\n\n Rerun UMAP\n\n") 
scdata <- FindNeighbors(scdata, reduction = "integrated.cca", dims = 1:30)
scdata <- FindClusters(scdata, resolution = 1)
scdata <- RunUMAP(scdata, dims = 1:30)

saveRDS(object = scdata, file = "scdata.rds")

#-------------------------------------------------------------------------------
# 可视化
#-------------------------------------------------------------------------------

# Visualization
# p1 <- DimPlot(scdata, reduction = "umap", group.by = c("stim", "seurat_annotations"))
# p2 <- DimPlot(scdata, reduction = "umap", split.by = "stim")





