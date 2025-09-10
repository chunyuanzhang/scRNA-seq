

if mkref == True:

    rule filter_gtf:
        input:
            gtf = scrna_gtf
        output:
            gtf = scrna_reference_path + scrna_genomename + ".filtered.gtf"
        threads:
            1
        log:
            "logs/scRNA/filter_gtf.log"
        shell:
            """
            {cellranger}
            cellranger mkgtf \
                {input.gtf} \
                {output.gtf} \
                --attribute=gene_biotype:protein_coding \
                --attribute=gene_biotype:lncRNA \
                --attribute=gene_biotype:antisense \
                --attribute=gene_biotype:IG_LV_gene \
                --attribute=gene_biotype:IG_V_gene \
                --attribute=gene_biotype:IG_V_pseudogene \
                --attribute=gene_biotype:IG_D_gene \
                --attribute=gene_biotype:IG_J_gene \
                --attribute=gene_biotype:IG_J_pseudogene \
                --attribute=gene_biotype:IG_C_gene \
                --attribute=gene_biotype:IG_C_pseudogene \
                --attribute=gene_biotype:TR_V_gene \
                --attribute=gene_biotype:TR_V_pseudogene \
                --attribute=gene_biotype:TR_D_gene \
                --attribute=gene_biotype:TR_J_gene \
                --attribute=gene_biotype:TR_J_pseudogene \
                --attribute=gene_biotype:TR_C_gene >>{log} 2>&1
            """


    rule cellranger_mkref:
        input:  
            reference = scrna_reference,
            genes = scrna_reference_path + scrna_genomename + ".filtered.gtf"
        output:
            fa = genome_path + "fasta/genome.fa",     # 提前创建文件夹会导致报错
            fai = genome_path + "fasta/genome.fa.fai",
            genes = genome_path + "genes/genes.gtf.gz",
            star = genome_path + "star/geneInfo.tab"
        log:
            "logs/scRNA/cellranger_mkref.log"
        params:
            genomename = scrna_genomename,
            outdir = scrna_reference_path,
            genes = scrna_genomename + ".filtered.gtf",
        threads:
            32
        shell:
            """
            {cellranger}
            cd {params.outdir} &&
            if [[ -d {params.genomename} ]]; then rm -rf {params.genomename} ; fi && 
            cellranger mkref \
                --genome={params.genomename} \
                --fasta={input.reference} \
                --genes={params.genes} \
                --nthreads {threads} && 
            cd - 
            """

else:
    rule mkln_index:
        input:
            reference = scrna_reference
        output:
            fa = genome_path + "fasta/genome.fa",     # 提前创建文件夹会导致报错
            fai = genome_path + "fasta/genome.fa.fai",
            genes = genome_path + "genes/genes.gtf.gz",
            star = genome_path + "star/geneInfo.tab"
        params:
            genome_path = genome_path,
            raw_path = os.path.abspath(scrna_reference)
        threads:
            1
        shell:
            """
            cd {params.genome_path}
            ln -s {params.raw_path}
            cd -
            """


rule count:
    input:
        folder = os.path.abspath(scrna_fastq_path + "{sample}"),
        fa = genome_path + "fasta/genome.fa"
    output:
        web_summary = scrna_count_path + "{sample}/outs/web_summary.html",   # 提前创建目录会导致报错
        metrics_summary = scrna_count_path + "{sample}/outs/metrics_summary.csv",
        cloupe = scrna_count_path + "{sample}/outs/cloupe.cloupe",
        raw = scrna_count_path + "{sample}/outs/raw_feature_bc_matrix.h5",
        filtered = scrna_count_path + "{sample}/outs/filtered_feature_bc_matrix.h5",
        molecule_info = scrna_count_path + "{sample}/outs/molecule_info.h5"
    params:
        genome_path = os.path.abspath(genome_path),
        outdir = "{sample}",
        scrna_count_path = scrna_count_path,
        include_introns = include_introns,
        output_bam = f"{output_bam}"
    log:
        "logs/scRNA/count.{sample}.log"
    threads:
        32
    shell:
        """
        {cellranger}
        cd {params.scrna_count_path} 
        if [[ -d {params.outdir} ]]; then rm -rf {params.outdir} ; fi 
        if [[ -f __{wildcards.sample}.mro ]] ; then rm __{wildcards.sample}.mro ; fi 
        cellranger count \
            --id {wildcards.sample} \
            --create-bam {params.output_bam} \
            --transcriptome {params.genome_path} \
            --fastqs {input.folder} \
            --include-introns {params.include_introns} \
            --localcores {threads} 
        cd - 
        """


localrules: cellranger_csv
rule cellranger_csv:
    input:
        filtered_h5 = expand(scrna_count_path + "{sample}/outs/filtered_feature_bc_matrix.h5", sample = samples.index),
        #molecule_info = expand(scrna_count_path + "{sample}/outs/molecule_info.h5", sample = samples.index),
        samplefile = samplefile
    output:
        csv = cwd + "/config/aggregation.csv"
    params:
        scrna_count_path = scrna_count_path,
        cwd = cwd
    shell:
        """
        echo "sample_id,molecule_h5" > {output.csv}
        for file in {input.filtered_h5}; do 
            echo {params.cwd}/$file | awk -F '/' '{{print $(NF-2)","$0}}' >> {output.csv}
        done
        """


# rule Seurat_sample:
#     input:
#         filtered_h5 = scrna_count_path + "{sample}/outs/filtered_feature_bc_matrix.h5"
#     output:
#         seurat_h5 = scrna_cluster_path + "{sample}/{sample}.rds"
#     params:
#         resolution = resolution,
#         outdir = scrna_cluster_path + "{sample}"
#     log:
#         "logs/scRNA/Seurat_sample.{sample}.log"
#     threads:
#         1
#     shell:
#         """
#         {Seurat5}
#         Rscript {Seurat} \
#             --sampleid {wildcards.sample} \
#             --resolution {params.resolution} \
#             --h5file {input.filtered_h5} \
#             --percentMT 5 \
#             --outdir {params.outdir} >>{log} 2>&1
#         """


rule Creat:
    input:
       h5 = expand(scrna_count_path + "{sample}/outs/filtered_feature_bc_matrix.h5", sample = samples.index),
       csv = cwd + "/config/aggregation.csv",
       samplefile = samplefile
    output:
        donefile = scrna_creat_path + "creat.done"
    log:
        "logs/scRNA/Creat.log"
    params:
        outdir = scrna_creat_path,
        percentMT = percentMT,
        cluster_resolution = cluster_resolution
    threads:
        1
    shell:
        """
        {Seurat5}
        Rscript {Seurat_Create} \
            --samplefile {input.samplefile} \
            --aggrtable {input.csv} \
            --outdir {params.outdir} \
            --percentMT {params.percentMT} \
            --cluster_resolution {params.cluster_resolution}  >{log} 2>&1
        """ 


rule Cluster:
    input:
        donefile = scrna_creat_path + "creat.done"
    output:
        donefile = scrna_cluster_path + "cluster.done"
    params:
        rds = scrna_path + "allsample.rds",
        outdir = scrna_cluster_path
    log:
        "logs/scRNA/Cluster.log"
    threads:
        1
    shell:
        """
        {Seurat5}
        Rscript {Seurat_Cluster} \
            --rds {params.rds} \
            --outdir {params.outdir}  >{log} 2>&1
        """

rule Seurat_Cluster_allsample_vis:
    input:
        done = scrna_cluster_path + "cluster.done",
        sampleandpopulation = samplefile
    output:
        csv = scrna_cluster_path + "01.umap_unintegrated.png"
    log:
        "logs/scRNA/Cluster.vis.log"
    params:
        rds = scrna_path +  "allsample.rds",
        outdir = scrna_cluster_path,
        by = by
    threads:
        1
    shell:
        """
        {Seurat5} 
        Rscript {Seurat_vis} \
            --sampleandpopulation {input.sampleandpopulation} \
            --rds {params.rds} \
            --by {params.by} \
            --step Cluster \
            --outdir {params.outdir}  >{log} 2>&1
        """

if species == "human" or species == "mouse":
    rule SingleR_celltyping:
        input:
            done = scrna_cluster_path + "cluster.done"
        output:
            done = scrna_celltype_path + "celltype.done"
        log:
            "logs/scRNA/celltype.log"
        params:
            rds = scrna_path + "allsample.rds",
            species = species,
            outdir = scrna_celltype_path
        threads:
            1
        shell:
            """
            {Seurat5}
            Rscript {CellTyping} \
                --rds {params.rds} \
                --species {params.species} \
                --outdir {params.outdir}  >{log} 2>&1
            """
else:
    rule SingleR_celltyping:
        input:
            done = scrna_cluster_path + "cluster.done",
            Orthologous = {Orthologous}
        output:
            done = scrna_celltype_path + "celltype.done"
        log:
            "logs/scRNA/celltype.log"
        params:
            rds = scrna_path + "allsample.rds",
            outdir = scrna_celltype_path,
            species = species
        threads:
            1
        shell:
            """
            {Seurat5}
            Rscript {CellTyping} \
                --rds {params.rds} \
                --species {params.species} \
                --Orthologous {input.Orthologous} \
                --outdir {params.outdir}  >{log} 2>&1
            """


rule SingleR_celltyping_vis:
    input:
        done = scrna_celltype_path + "celltype.done"
    output:
        png = scrna_celltype_path + "01.umap_integrated_celltype.png"
    log:
        "logs/scRNA/celltype.vis.log"
    params:
        rds = scrna_path + "allsample.rds",
        outdir = scrna_celltype_path
    threads:
        1
    shell:
        """
        {Seurat5} 
        Rscript {Seurat_vis} \
            --rds {params.rds} \
            --outdir {params.outdir} \
            --step CellType  >{log} 2>&1
        """


# 含有marker的表格不保存在Seurat对象中，该步骤不输出rds文件
rule Marker:
    input:
        done = scrna_cluster_path + "cluster.done"
    output:
        done = scrna_marker_path + "marker.done",
        Marker_table = scrna_marker_path + "01.global_DEGs.csv",
        topn_Marker_table = scrna_marker_path +  "02.topn_markers.csv"
    log:
        "logs/scRNA/marker.log"
    params:
        rds = scrna_path + "allsample.rds",
        topn = topn,
        by = by,
        p_val = p_val,
        avg_log2FC = avg_log2FC,
        min_pct1 = min_pct1,
        outdir = scrna_marker_path 
    threads:
        1
    shell:
        """
        {Seurat5}
        Rscript {Seurat_Marker} \
            --rds {params.rds} \
            --by {params.by} \
            --topn {params.topn} \
            --p_val {params.p_val} \
            --avg_log2FC {params.avg_log2FC} \
            --min_pct1 {params.min_pct1} \
            --outdir {params.outdir}  >{log} 2>&1
        """

rule Seurat_Marker_allsample_vis:
    input:
        done = scrna_marker_path + "marker.done",
        Marker_table = scrna_marker_path + "01.global_DEGs.csv",
        topn_Marker_table = scrna_marker_path +  "02.topn_markers.csv"
    output:
        heatmap = scrna_marker_path +  "01.topn_DoHeatmap.png"
    log:
        "logs/scRNA/marker.vis.log"
    params:
        rds = scrna_path + "allsample.rds",
        outdir = scrna_marker_path,
        by = by
    threads:
        1
    shell:
        """
        {Seurat5}
        Rscript {Seurat_vis} \
            --rds {params.rds} \
            --by {params.by} \
            --step Marker \
            --Marker_table {input.Marker_table} \
            --topn_Marker_table {input.topn_Marker_table} \
            --outdir {params.outdir}  >{log} 2>&1
        """


rule Seurat_Differential_Expression:
    input:
        done = scrna_celltype_path + "celltype.done"
    output:
        # rds = scrna_DE_path + "{diff_group}/{diff_group}.rds",  #差异分析的表格没有添加到Seurat对象中，无需重复输出rds文件
        DE = scrna_DE_path + "{diff_group}/01.DEgenes.csv",
        topnDE = scrna_DE_path + "{diff_group}/02.topn_DEgenes.csv"
    params:
        rds = scrna_path + "allsample.rds",
        case = lambda w: diff_params.loc["{}".format(w.diff_group)]["case"],
        control = lambda w: diff_params.loc["{}".format(w.diff_group)]["control"],
        difftype = lambda w: diff_params.loc["{}".format(w.diff_group)]["difftype"],
        outdir = scrna_DE_path + "{diff_group}",
        topn = topn,
        p_val = p_val,
        min_pct1 = min_pct1,
        avg_log2FC = avg_log2FC
    log:
        "logs/scRNA/Seurat_Differential_Expression.{diff_group}.log"
    threads:
        1
    shell:
        """
        {Seurat5}
        Rscript {Seurat_DE} \
            --rds {params.rds} \
            --case {params.case} \
            --control {params.control} \
            --difftype {params.difftype} \
            --topn {params.topn} \
            --p_val {params.p_val} \
            --avg_log2FC {params.avg_log2FC} \
            --min_pct1 {params.min_pct1} \
            --outdir {params.outdir}  >{log} 2>&1
        """

rule Seurat_Differential_Expression_vis:
    input:
        # rds = scrna_DE_path + "{diff_group}/{diff_group}.rds",
        DE = scrna_DE_path + "{diff_group}/01.DEgenes.csv",
        topnDE = scrna_DE_path + "{diff_group}/02.topn_DEgenes.csv"
    output:
        png = scrna_DE_path + "{diff_group}/01.top_Marker_VlnPlot_high.png"
    log:
        "logs/scRNA/DE.{diff_group}.log"
    threads:
        1
    params:
        rds = scrna_path + "allsample.rds",
        case = lambda w: diff_params.loc["{}".format(w.diff_group)]["case"],
        control = lambda w: diff_params.loc["{}".format(w.diff_group)]["control"],
        difftype = lambda w: diff_params.loc["{}".format(w.diff_group)]["difftype"],
        outdir = scrna_DE_path + "{diff_group}",
        cluster_DEgenes_all = scrna_DE_path + "{diff_group}/03.DEgenes_Clusters.csv"
    shell:
        """
        {Seurat5}
        Rscript {Seurat_vis} \
            --rds {params.rds} \
            --Marker_table {input.DE} \
            --topn_Marker_table {input.topnDE} \
            --step DiffExp \
            --case {params.case} \
            --control {params.control} \
            --difftype {params.difftype} \
            --cluster_DEgenes_all {params.cluster_DEgenes_all} \
            --outdir {params.outdir}  >{log} 2>&1
        """



# rule subcelltype:
#     input:
#         rds = scrna_path + "allsample.rds"
#     output:
#         rds = subcelltype_path + "subset.rds",
#         Marker_table = subcelltype_path + "01.global_DEGs.csv",
#         topn_Marker_table = subcelltype_path + "02.topn_markers.csv"
#     params:
#         cluster_resolution = cluster_resolution,
#         target_cluster = target_cluster,
#         outdir = subcelltype_path,
#         topn = topn,
#         p_val = p_val,
#         avg_log2FC = avg_log2FC,
#         min_pct1 = min_pct1
#     shell:
#         """
#         {Seurat5} 
#         Rscript {subcelltype} --rds {input.rds} \
#             --target_cluster {params.target_cluster} \
#             --cluster_resolution {params.cluster_resolution} \
#             --outdir {params.outdir}

#         Rscript {Seurat_Marker} \
#             --rds {output.rds} \
#             --topn {params.topn} \
#             --p_val {params.p_val} \
#             --avg_log2FC {params.avg_log2FC} \
#             --min_pct1 {params.min_pct1} \
#             --outdir {params.outdir}  

#         Rscript {Seurat_vis} \
#             --rds {output.rds} \
#             --step Marker \
#             --Marker_table {output.Marker_table} \
#             --topn_Marker_table {output.topn_Marker_table} \
#             --outdir {params.outdir} 
#         """









