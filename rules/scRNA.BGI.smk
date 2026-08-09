
rule dnbc4tools:
    input:
        fastq = scrna_fastq_path + "{sample}/http-{sample}"
    output: 
        h5ad = scrna_count_path + "{sample}/outs/filter_feature.h5ad"
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



