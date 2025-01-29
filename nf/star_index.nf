nextflow.enable.dsl=2

conda_env_path = "${params.cogent_root}" + '/config/env/star.yaml'
star_index_script = "${params.cogent_root}" + '/bin/star_index.sh'

process STAR_INDEX {
    conda conda_env_path
    storeDir "${params.genome_dir}"

    input:
        path (genome)
        path (gtf)

    output:
        path("${params.genome_name}/star_index"), emit: star_index

    script:
        """
        bash ${star_index_script} ${params.genome_name} ${task.cpus} \
        ${genome} ${gtf} "star_index" \
        ${task.memory.toBytes() - 1024*1024*1024}
        """
}
