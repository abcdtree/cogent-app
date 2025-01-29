
conda_env_path = "${params.cogent_root}" + '/config/env/salmon_report.yaml'
salmon_matrix_script = "${params.cogent_root}" + '/lib/salmon/create_genematrix.py'

process SALMON_TRANSCRIPT_MATRIX {
    conda conda_env_path

    input:
        path (input_dir)
        val (umi_only)

    output:
        path("prelim_analyze_transcript_matrix*.csv"), emit: transcript_matrix

    script:

        def output_file = umi_only ? "prelim_analyze_transcript_matrix_umis.csv" : "prelim_analyze_transcript_matrix.csv"

        """
        python3 ${salmon_matrix_script} --input_dir \$PWD --type 'transcripts' \
            -o ${output_file}
        """
}
