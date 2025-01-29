nextflow.enable.dsl=2

conda_env_path = "${params.cogent_root}" + '/config/env/analyze_report.yaml'
modify_fusion_script = "${params.cogent_root}" + '/lib/fusion/modify_fusion_results.py'

process MODIFY_FUSION_OUT {
    conda conda_env_path

    input:
        tuple val(sample_id), path(junction_file)

    output:
        tuple val(sample_id), path ("*_updated_Chimeric.out.junction"), emit: modified_junction_file

    script:
        """
        python3 ${modify_fusion_script} --input_file ${junction_file} --output_file ${sample_id}_updated_Chimeric.out.junction
        """
}
