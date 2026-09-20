
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
        40
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


rule Seurat5:
    input:
        barcodes = expand(scrna_count_path + "{sample}/outs/filter_matrix/barcodes.tsv.gz", sample=samples.index),
        features = expand(scrna_count_path + "{sample}/outs/filter_matrix/features.tsv.gz", sample=samples.index),
        matrix = expand(scrna_count_path + "{sample}/outs/filter_matrix/matrix.mtx.gz", sample=samples.index)
    output:
        rds = "result/scdata.rds"
    threads:
        16
    params:
        percentMT = percentMT,
        MTpattern = MTpattern,
        samplefile = samplefile
    shell:
        """
        ~/tools/Seurat/bin/Rscript scripts/Seurat.R \
            --SampleFile {params.samplefile} \
            --MTpattern {params.MTpattern} \
            --percentMT {params.percentMT}
        """

