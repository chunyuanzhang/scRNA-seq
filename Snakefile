import os
include: "rules/0.Common.smk"

def outfiles():
    files=[]

    if scRNA_platform == "10X":
        include: "rules/scRNA.smk"
        files.append(genome_path + "fasta/genome.fa.fai")
        if "count" in step:
            files.append(expand(scrna_count_path + "{sample}/outs/filtered_feature_bc_matrix.h5", sample=samples.index))
        if "create" in step:
            files.append(scrna_creat_path + "creat.done")
        if "cluster" in step:
            files.append(scrna_cluster_path + "cluster.done")
            files.append(scrna_cluster_path + "01.umap_unintegrated.png")
        if "celltype" in step:
            files.append(scrna_celltype_path + "celltype.done")
            files.append(scrna_celltype_path + "01.umap_integrated_celltype.png")
        if "marker" in step:
            files.append(scrna_marker_path + "marker.done")
            files.append(scrna_marker_path +  "01.topn_DoHeatmap.png")
        if "diffexp" in step:
            #files.append(expand(scrna_DE_path + "{diff_group}/{diff_group}.rds", diff_group = diff_groups))
            files.append(expand(scrna_DE_path + "{diff_group}/01.DEgenes.csv", diff_group = diff_groups))
            files.append(expand(scrna_DE_path + "{diff_group}/02.topn_DEgenes.csv", diff_group = diff_groups))
            files.append(expand(scrna_DE_path + "{diff_group}/01.top_Marker_VlnPlot_high.png", diff_group = diff_groups))
        # if "subcelltype" in step:
        #     files.append(subcelltype_path + "subset.rds")

    if scRNA_platform == "BGI":
        include: "rules/scRNA.BGI.smk"
        if "count" in step:
            files.append(expand(scrna_count_path + "{sample}/outs/filter_feature.h5ad", sample=samples.index))
    
    # print(files)
    return files

rule all:
    input: outfiles()
