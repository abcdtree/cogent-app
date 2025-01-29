nextflow.enable.dsl=2

conda_env_path = "${params.cogent_root}" + '/config/env/star_fusion.yaml'

process STAR_FUSION {
    conda conda_env_path
    publishDir "${params.output_dir}", mode: 'copy'

    input:
        tuple val(sample_id), path(junction_file), path(fusion_index)

    output:
        path("06-star_fusion/${sample_id}_fusion_predictions.abridged.tsv"), emit: star_fusion_abridged_tsv
        path("06-star_fusion/${sample_id}_fusion_predictions.tsv"), emit: star_fusion_tsv
        path("logs/star_fusion/*.log"), emit: fusion_log

    script:
        """
        (STAR-Fusion \
            --genome_lib_dir "${fusion_index}" \
            -J "${junction_file}" \
            --output_dir "06-star_fusion/fusion_out_${sample_id}" \
            --CPU ${task.cpus} \
            --max_sensitivity \
        ) 2>&1 | tee ${sample_id}_fusion.log

        mv 06-star_fusion/fusion_out_${sample_id}/star-fusion.fusion_predictions.abridged.tsv \
            06-star_fusion/${sample_id}_fusion_predictions.abridged.tsv

        mv 06-star_fusion/fusion_out_${sample_id}/star-fusion.fusion_predictions.tsv \
            06-star_fusion/${sample_id}_fusion_predictions.tsv
        
        mkdir -p logs/star_fusion
        mv ${sample_id}_fusion.log logs/star_fusion/. 

        """
}
