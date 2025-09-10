* 本仓库的内容修改自赵老师的代码仓库，部署在家禽团队自己的服务器上
* 由于gitlab账号被封，当前仓库位于github上
* 【注意】由于 celllranger 版本升级，输出bam文件的参数发生了变化，当前使用 --create-bam

## 下载流程

```
git clone git@github.com:chunyuanzhang/scRNA-seq.git

```



## 提交运行
* 湘湖实验室家禽团队自己的服务器只是一个单独的服务器，没有集群，也没有调度系统，snakemake -j32 提交即可


# 单细胞分析

## 配置 config/sampleandpopulation.tsv 文件
任何分析都需要现配置 config/sampleandpopulation.tsv 文件，该文件共两类，第一列为样本名，第二列为群体名，中间用tab键分隔，每行一个样本，样本不允许重复

```
JH1,JH
JH10,JH
JH2,JH
JH3,JH
JH4,JH
HBei1,HBM
HBei10,HBM
HBei2,HBM
HBei3,HBM
```

## 配置 config/config.yaml 文件



* **单细胞分析使用ggplot2包版本不得超过 3.4.4, 否则无法正常绘制图片**

## 单细胞分析步骤

单细胞分析步骤须从 count,create,cluster,celltype,marker,diffexp,alternativesplincing 中进行选择

* count:  cellranger比对测序数据到参考基因组上，得到基因表达矩阵
* create: 创建Seurat对象，将矫批次效应后将数据整合到一起
* cluster: 降维聚类
* celltype: 确定细胞类型 【注意：自动鉴定的细胞类型不靠谱，必须自行手动鉴定】
* marker: 找标记基因 - 找标记基因的策略是在cluster间进行 1 VS all othrers 的比较
* diffexp: 差异表达可以在clutter间进行，也可以在样本间进行，其他任何在meta.data表格中存在的分类间的差异均可进行差异分析

cluster、celltype、marker分析的内容均会保存在Seurat对象中，完成对应部分分析后也在会输出表格 meta.data.csv 中添加响应的信息。

## 差异分析

差异分析需要配置指定差异分析的文件 config/scRNA_differential_expression.csv     
文件内容示例如下：

```
case,control,difftype
A1312-1-10XSC3,A1412-1-10XSC3,sampleid    # 样本间比较
A1314-3-10XSC3,A1414-3-10XSC3,sampleid    
A1312-1-10XSC3,A1314-1-10XSC3,sampleid
A1412-1-10XSC3,A1414-3-10XSC3,sampleid
13,14,groupid                             # 组间比较
1,2,seurat_clusters                       # cluster间比较
0,2,seurat_clusters
```

* 差异分析可以在任何meta.data.csv表格中出现过的分组之间进行比较

