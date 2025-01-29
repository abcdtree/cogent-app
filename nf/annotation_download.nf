nextflow.enable.dsl=2

conda_env_path = "${params.cogent_root}" + '/config/env/gzip.yaml'
genome_download_script = "${params.cogent_root}" + '/bin/genome_download.sh'

process ANNOTATION_DOWNLOAD {
    conda conda_env_path
    storeDir "${params.genome_dir}"

    input:
        val(source_url)

    output:
        path("${params.genome_name}/${params.genome_name}.annotation.gtf"), emit: gtf

    script:
        """
        bash ${genome_download_script} \
        ${params.genome_name} \
        ${source_url} \
        "annotation.gtf"
        """
}
