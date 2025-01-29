nextflow.enable.dsl=2

fusion_report_script = "${params.cogent_root}" + '/lib/fusion/fusion_report.py'
conda_env_path = "${params.cogent_root}" + '/config/env/analyze_report.yaml'

process STAR_FUSION_REPORT {
    publishDir "${params.output_dir}", mode: 'copy'
    conda conda_env_path 

    input:
        path (star_predictions_file)

    output:
        path("count_matrices/star_junction_matrix.csv"), emit: junction_matrix
        path("count_matrices/star_spanning_matrix.csv"), emit: spanning_matrix

    script:
        """
        python3 ${fusion_report_script} --input_dir \$PWD --output_dir count_matrices
        """
}
