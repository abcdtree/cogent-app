nextflow.enable.dsl=2

conda_env_path = "${params.cogent_root}" + '/config/env/fastqc.yaml'
fastqc_script = "${params.cogent_root}" + '/bin/fastqc_demux.sh'

process FASTQC_DEMUX {
    conda conda_env_path
    publishDir "${params.output_dir}", mode: 'move'

    input:
        path(read1)
        path(read2)

    output:
        path ("00a-fastqc/fastqc_logs/*.zip"), emit: fastqc_zip
        path ("00a-fastqc/fastqc_logs/*.html"), emit: fastqc_html

    script:

        """
        bash ${fastqc_script} "00a-fastqc" "${read1}" "${read2}"
        """
}
