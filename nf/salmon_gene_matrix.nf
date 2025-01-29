conda_env_path = "${params.cogent_root}" + '/config/env/salmon_report.yaml'
salmon_matrix_script = "${params.cogent_root}" + '/lib/salmon/create_genematrix.py'


process SALMON_GENE_MATRIX {
    conda conda_env_path
    publishDir "${params.output_dir}", mode: 'copy'
    

    input:
        path (salmon_results)
        path (gene_info)
        val (umi_only)

    output:
        path("count_matrices/analyze_gene_matrix*.csv"), emit: gene_matrix
        path("count_matrices/gene_info.csv"), emit: salmon_gene_info, optional:true

    script:

        def output_file = umi_only ? "analyze_gene_matrix_umis.csv" : "analyze_gene_matrix.csv"
        def quants_folder = umi_only ? "*_trimmed_umionly" : "*_trimmed"
        def cp_gene_cmd = umi_only ? " " : "cp ${gene_info} count_matrices/."

        """
        python3 ${salmon_matrix_script} --input_dir \$PWD --type 'genes' -r count_matrices -o ${output_file}

        ${cp_gene_cmd}
        """
}
