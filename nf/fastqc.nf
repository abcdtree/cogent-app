nextflow.enable.dsl=2

conda_env_path = "${params.cogent_root}" + '/config/env/fastqc.yaml'
fastqc_script = "${params.cogent_root}" + '/bin/fastqc.sh'

process FASTQC {
    conda conda_env_path
    publishDir "${params.output_dir}", mode: 'move'

    input:
        tuple val(sample_id), path(reads)

    output:
        tuple val(sample_id), path ("01-fastqc_pre_trim/fastqc_${sample_id}_logs/*.zip"), emit: fastqc_zip
        tuple val(sample_id), path ("01-fastqc_pre_trim/fastqc_${sample_id}_logs/*.html"), emit: fastqc_html

    script:

        """
        bash ${fastqc_script} "01-fastqc_pre_trim" "$sample_id" "${reads[0]}" "${reads[1]}"
        """
}
