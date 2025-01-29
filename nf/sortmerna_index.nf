nextflow.enable.dsl=2

conda_env_path = "${params.cogent_root}" + '/config/env/sortmerna.yaml'


process SORTMERNA_INDEX {
    conda conda_env_path
    storeDir "${params.genome_dir}"

    input:
        path(fastas)

    output:
        path("${params.genome_name}/sortmernadb/idx")
        path("${params.genome_name}/sortmernadb/*_index.log")
        path("${params.genome_name}/sortmernadb/ref")

    script:

        """
        mkdir -p ${params.genome_name}/sortmernadb/idx
        mkdir -p ${params.genome_name}/sortmernadb/ref

        (sortmerna \\
            ${'--ref '+fastas.join(' --ref ')} \\
            --threads ${task.cpus} \\
            --idx-dir ${params.genome_name}/sortmernadb/idx \\
            --index 1
        ) 2>&1 | tee ${params.genome_name}/sortmernadb/sortmerna_index.log

        cp ${fastas} ${params.genome_name}/sortmernadb/ref/.

        """
}
