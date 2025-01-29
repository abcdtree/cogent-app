nextflow.enable.dsl=2

conda_env_path = "${params.cogent_root}" + '/config/env/salmon.yaml'

process SALMON_QUANT {
    conda conda_env_path
    publishDir "${params.output_dir}", mode: 'copy'

    input:
        tuple val(sample_id), path(bam), path(transcript_fa)
        val(umi_only)

    output:
        path("05-salmon_counts/${sample_id}*_quant.sf"), emit: salmon_transcript_results
        path("logs/salmon/${sample_id}*_salmon_quant.log"), emit: salmon_logs
        path("logs/salmon/${sample_id}*_meta_info.json"), emit:salmon_meta_info

    script:

        def strandedness =  'A'

        def output_dir = umi_only ? "${sample_id}_umionly" : "${sample_id}"
        def log_dir = umi_only ? "${sample_id}_salmon_umionly" : "${sample_id}_salmon"
        def min_frags = 1


        """
        mkdir -p 05-salmon_counts
        mkdir -p logs/salmon

        salmon quant \\
            --threads ${task.cpus} \\
            --libType=${strandedness} \\
            --useVBOpt \\
            --writeUnmappedNames \\
            --minAssignedFrags ${min_frags} \\
            -a ${bam} \\
            -t ${transcript_fa} \\
            -o 05-salmon_counts/${output_dir} \\
        
        mv 05-salmon_counts/${output_dir}/quant.sf 05-salmon_counts/${output_dir}_quant.sf

        mv 05-salmon_counts/${output_dir}/aux_info/meta_info.json \\
            logs/salmon/${output_dir}_meta_info.json
        
        mv 05-salmon_counts/${output_dir}/logs/salmon_quant.log \\
            logs/salmon/${output_dir}_salmon_quant.log 

        """
}
