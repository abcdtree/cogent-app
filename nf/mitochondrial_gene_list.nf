nextflow.enable.dsl=2

process MITOCHONDRIAL_GENE_LIST {
    storeDir "${params.genome_dir}"

    input:
        path(gene_list)

    output:
        path("${params.genome_name}/mito_reads/*.txt")

    script:

        """
        mkdir -p ${params.genome_name}/mito_reads

        cp ${gene_list} ${params.genome_name}/mito_reads/. 

        """
}
