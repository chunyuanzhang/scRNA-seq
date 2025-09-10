import pandas as pd
import numpy as np
import os
import sys
import re
import yaml
import json
# from snakemake.utils import validate
from collections import defaultdict
import shutil
from itertools import combinations
from itertools import permutations
import gzip

### load config file======================================================================================================
configfile: "config/config.yaml"
configfile: "config/cluster.yaml"

cwd = os.getcwd()

### step =================================================================================================================
step = config["step"]
if step is not None:
    step = step.lower()
    step = step.replace(" ","").replace("\t","").split(",")
else:
    step = ""


### species ==============================================================================================================
species = config["species"]

### samples and populations ==============================================================================================
samplefile = config["samplefile"]
if os.path.exists(samplefile):
    # 样本和群体的对应表格
    samples = pd.read_csv(samplefile, sep=",", dtype=str, index_col = False)
    samples.dropna(how='all', inplace=True) # 删除可能存在的空行
    samples = samples[ ~np.array([s.startswith("#") for s in samples.sampleid.to_list()])]
    samples.index = samples.sampleid
    # 将群体以字典的方式存储，方便任何时候提取
    group_dict = samples.groupby("groupid")['sampleid'].apply(list).to_dict()
    groups = list(group_dict.keys())

    wildcard_constraints:
        samples = "|".join(samples.sampleid)
    wildcard_constraints:
        groups = "|".join(group_dict.keys())

else:
    print(f"您尚未准备好样本文件{samplefile}")
    os._exit(0)

####################################################
# single cell analysis
####################################################

scrna_reference = config["scrna_reference"]
scrna_gtf = config["scrna_gtf"]
scrna_genomename = config["scrna_genomename"]

# 检查参考基因组是否存在
if not scrna_reference.endswith("fasta/genome.fa") or not scrna_gtf.endswith("genes/genes.gtf.gz"):
    mkref = True

# 参数提取和传递
include_introns = config["analysis"]["count"]["include_introns"]
include_introns = str(include_introns).lower()
output_bam = config["analysis"]["count"]["output_bam"]
output_bam = str(output_bam).lower()
percentMT = config["analysis"]["create"]["percentMT"]
cluster_resolution = config["analysis"]["cluster"]["resolution"]
topn = config["analysis"]["marker"]["topn"]
by = config["analysis"]["marker"]["by"]
p_val = config["analysis"]["diffexp"]["p_val"]
avg_log2FC = config["analysis"]["diffexp"]["avg_log2FC"]
min_pct1 = config["analysis"]["diffexp"]["min_pct1"]
Orthologous = config["analysis"]["celltype"]["Orthologous"]




version_10X_chemistry = config["analysis"]["AlternativeSplincing"]["version_10X_chemistry"]
if version_10X_chemistry == "V3":
    soloUMIlen = 12
if version_10X_chemistry == "V2":
    soloUMIlen = 10


##############################
# 文件夹路径
##############################

scrna_reference_path = "result/reference/"
genome_path = scrna_reference_path + scrna_genomename + "/"
scrna_path = "result/"
scrna_fastq_path = "result/01.fastq/"
scrna_count_path = "result/02.Count/"
scrna_creat_path = "result/03.Create/"
scrna_cluster_path = "result/04.Cluster/"
scrna_celltype_path = "result/05.CellType/"
scrna_marker_path = "result/06.Marker/"
scrna_DE_path = "result/07.DiffExp/"
scrna_AS_path = "result/10.AlternativeSplincing/"
scrna_STAR_index = genome_path + "STAR_index/"

##############################
# 声明工具
##############################

cellranger = config["tools"]["cellranger"]
