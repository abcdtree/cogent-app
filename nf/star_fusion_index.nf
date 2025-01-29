nextflow.enable.dsl=2

conda_env_path = "${params.cogent_root}" + '/config/env/tar.yaml'
star_fusion_index_script = "${params.cogent_root}" + '/bin/star_fusion_index.sh'

process STAR_FUSION_INDEX {
    conda conda_env_path
    storeDir "${params.genome_dir}"

    output:
        path("${params.genome_name}/fusion_index")

    script:
        if(!params.fusion_lib_url) {
            error "Cannot install fusion library because --fusion_lib_url is not set."
        }
        """
        bash ${star_fusion_index_script} ${params.genome_name} \
        ${params.fusion_lib_url} \
        "GRCh38_gencode_v37_CTAT_lib_Mar012021.plug-n-play" \
        """
}
