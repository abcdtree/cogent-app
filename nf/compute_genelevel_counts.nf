nextflow.enable.dsl=2

compute_counts = "${params.cogent_root}" + '/lib/salmon/compute_genelevel_counts.R'
conda_env_path = "${params.cogent_root}" + '/config/env/gene_level_counts.yaml'

process COMPUTE_GENELEVEL_COUNTS {
    publishDir "${params.output_dir}", mode: 'copy'
    conda conda_env_path 

    input:
        path (salmon_transcript_matrix)
        path (gene_transcript_map)
        path (gene_info)
        path (transcript_info)
        val (umi_only)

    output:
        path ("count_matrices/analyze_incl_introns_genematrix*.csv"), emit: intron_gene_matrix
        path ("count_matrices/analyze_gene_matrix*.csv"), emit: gene_matrix
        path ("count_matrices/analyze_transcript_matrix*.csv"), emit: transcript_matrix
        path ("count_matrices/gene_info.csv"), emit: salmon_gene_info, optional:true
        path ("count_matrices/transcript_info.csv"), emit: salmon_transcript_info, optional:true

    script:

        def umi_flag = umi_only ? "true" : "false"
        def cp_transcript_cmd = umi_only ? " " : "cp ${transcript_info} count_matrices/."
        def cp_gene_cmd = umi_only ? " " : "cp ${gene_info} count_matrices/."
        
        """
        mkdir -p count_matrices 

        Rscript ${compute_counts} \
            --gene_transcript_Map ${gene_transcript_map} \
            --salmon_transcriptCounts ${salmon_transcript_matrix} \
            --outPath count_matrices \
            --umiOnly ${umi_flag}
        
        ${cp_gene_cmd}
        ${cp_transcript_cmd}
        
        """
}
