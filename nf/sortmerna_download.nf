nextflow.enable.dsl=2

conda_env_path = "${params.cogent_root}" + '/config/env/wget.yaml'

process SORTMERNA_DOWNLOAD {
    conda conda_env_path
    storeDir "${params.genome_dir}"

    input:
        each(fasta_file)

    output:
        path("sortmernadb/ref/*.fasta"), emit: sortmerna_ref

    script:
        """
        mkdir -p sortmernadb/ref
        wget ${fasta_file} -P sortmernadb/ref
        """
}
