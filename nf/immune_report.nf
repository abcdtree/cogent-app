nextflow.enable.dsl=2

immune_report_script = "${params.cogent_root}" + '/lib/immune/immune_report.py'

process IMMUNE_REPORT {
    publishDir "${params.output_dir}", mode: 'copy'

    input:
        path(report)

    output:
        path("10-immune_profiling/immune_clonotype_matrix.csv"), emit: clonotype_matrix
        path("10-immune_profiling/immune_top3_clonotype_matrix.csv"), emit: top3_clonotype_matrix
        path("10-immune_profiling/immune_summary.csv"), emit: immune_summary
        path("10-immune_profiling/immune_top3_summary.csv"), emit: immune_top3_summary
        path("10-immune_profiling/immune_metadata.csv"), emit: immune_metadata
        path("10-immune_profiling/immune_top3_metadata.csv"), emit: immune_top3_metadata

    script:
        """
        python3 ${immune_report_script} --input_dir \$PWD ${report}
        """
}
