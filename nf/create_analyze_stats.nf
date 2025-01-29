nextflow.enable.dsl=2

stats_script = "${params.cogent_root}" + '/lib/report/create_analyze_stats.py'
conda_env_path = "${params.cogent_root}" + '/config/env/analyze_report.yaml'

process CREATE_ANALYZE_STATS {
    publishDir "${params.output_dir}", mode: 'copy'
    conda conda_env_path 

    input:
        path (cutadapt_log)
        path (sortmerna_log)
        path (salmon_log)
        path (star_log)
        path (star_reads_log)
        path (mito_reads)
        path (salmon_intron_matrix)
        path (salmon_gene_matrix)
        path (salmon_transcript_matrix)
        path (umitools_log)
        path (umitools_umionly_log)
        path (salmon_gene_matrix_umionly)
        path (salmon_transcript_matrix_umionly)
        path (salmon_intron_matrix_umionly)


    output:
        path "analyze_stats.csv", emit: stats_file
        path "analyze_stats_5pUMI_stats.csv", emit: umi_stats_file, optional: true

    script:
        def exp_type = "${params.type_of_experiment}"
        """
        python3 ${stats_script} --input_dir . -n ${params.gene_count_threshold} -t ${exp_type}
        """
}
