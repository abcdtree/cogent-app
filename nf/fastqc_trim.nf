nextflow.enable.dsl=2

conda_env_path = "${params.cogent_root}" + '/config/env/fastqc.yaml'
fastqc_script = "${params.cogent_root}" + '/bin/fastqc.sh'

process FASTQC_TRIM {
    conda conda_env_path
    publishDir "${params.output_dir}", mode: 'move'

    input:
        tuple val(sample_id), path(reads)

    output:
        tuple val(sample_id), path ("03-fastqc_post_trim/fastqc_${sample_id}_logs/*.zip"), emit: fastqc_zip
        tuple val(sample_id), path ("03-fastqc_post_trim/fastqc_${sample_id}_logs/*.html"), emit: fastqc_html

    script:

        def read1 = "${reads[0]}"
        def read2 = "${reads[1]}"
        """
        bash ${fastqc_script} "03-fastqc_post_trim" "${sample_id}" "$read1" "$read2"
        """
}
