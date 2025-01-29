nextflow.enable.dsl=2

conda_env_path = "${params.cogent_root}" + '/config/env/gffread.yaml'

process GFFREAD_INTRON {
    conda conda_env_path
    storeDir "${params.genome_dir}"

    input:
        path (genome)
        path (gtf)

    output:
        path("${params.genome_name}/${params.genome_name}.intron.transcripts.fa"), emit: intron_transcriptome_file

    script:
        """
        mkdir -p ${params.genome_name}

        gffread -w ${params.genome_name}/${params.genome_name}.intron.transcripts.fa -g ${genome} ${gtf}
        """
}
