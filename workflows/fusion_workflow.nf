#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// import fusion related modules
include { MODIFY_FUSION_OUT } from '../modules/modify_fusion_out'
include { STAR_FUSION } from '../modules/star_fusion'
include { STAR_FUSION_REPORT } from '../modules/star_fusion_report'


workflow FUSION_WORKFLOW {

    take:
        junction_file_ch

    main:
        fusion_index = file("${params.genome_dir}/${params.genome}/fusion_index")

        if ( !params.genome_dir.startsWith('s3:') ) {
                if ( !fusion_index.exists() ) {
                    exit 1, "STAR Fusion Index does not exist, please run add_genome with --fusion first."
                }
        }

        Channel.fromPath(fusion_index)
            .ifEmpty { error "STAR Fusion Index does not exist, please run add_genome with --fusion first." }
            .set{ fusion_index_ch }

        MODIFY_FUSION_OUT(junction_file_ch)

        MODIFY_FUSION_OUT.out.modified_junction_file
            .combine(fusion_index_ch)
            .set{ star_fusion_input_ch }

        STAR_FUSION(star_fusion_input_ch)

        STAR_FUSION.out.star_fusion_abridged_tsv
            .collect()
            .set { fusion_results_ch }

        STAR_FUSION_REPORT(fusion_results_ch)

    emit:
        junction_matrix = STAR_FUSION_REPORT.out.junction_matrix
        spanning_matrix = STAR_FUSION_REPORT.out.spanning_matrix

}
