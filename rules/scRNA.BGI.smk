
rule dnbc4tools:
    input:
        fastq = scrna_fastq_path + "{sample}/http-{sample}"
    output: 
        barcodes = scrna_count_path + "{sample}/outs/filter_matrix/barcodes.tsv.gz",
        features = scrna_count_path + "{sample}/outs/filter_matrix/features.tsv.gz",
        matrix = scrna_count_path + "{sample}/outs/filter_matrix/matrix.mtx.gz",
        # h5ad = temp(scrna_count_path + "{sample}/outs/filter_feature.h5ad"),
        bam = temp(scrna_count_path + "{sample}/outs/anno_decon_sorted.bam"),
        bai = temp(scrna_count_path + "{sample}/outs/anno_decon_sorted.bam.bai")
        # rawbarcodes = temp(scrna_count_path + "{sample}/outs/raw_matrix/barcodes.tsv.gz"),
        # rawfeatures = temp(scrna_count_path + "{sample}/outs/raw_matrix/features.tsv.gz"),
        # rawmatrix = temp(scrna_count_path + "{sample}/outs/raw_matrix/matrix.mtx.gz")
    threads: 
        64
    params:
        genomeDir = "~/zhangchunyuan/reference/bGalGal1_mat_broiler_GRCg7b/dnbc4/Chicken",
        name = "{sample}"
    shell:
        """
        module load dnbc4tools
        cd result/02.Count/

        dnbc4tools rna run \
            --fastqs ../01.fastq/{params.name} \
            --genomeDir {params.genomeDir} \
            --name {params.name} \
            --threads {threads}

        cd -
        """



