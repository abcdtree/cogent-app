conda_env_path = "${params.cogent_root}" + '/config/env/salmon_report.yaml'
salmon_matrix_script = "${params.cogent_root}" + '/lib/salmon/create_genematrix.py'


process SALMON_INTRON_MATRIX {
    conda conda_env_path
    publishDir "${params.output_dir}", mode: 'copy'

    input:
        path (input_dir)
        val(umi_only)

    output:
        path("count_matrices/analyze_incl_introns_genematrix*.csv"), emit: intron_gene_matrix

    script:

        def output_file = umi_only ? "analyze_incl_introns_genematrix_umis.csv" : "analyze_incl_introns_genematrix.csv"

        """
       
        python3 ${salmon_matrix_script} --input_dir \$PWD --type 'introns' -r count_matrices -o ${output_file}
        """
}