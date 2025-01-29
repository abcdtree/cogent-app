#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// import immune related modules
include { TRUST4 } from '../modules/trust4'
include { IMMUNE_REPORT } from '../modules/immune_report'

workflow IMMUNE_WORKFLOW {

    take:
        fastq_file_ch
    
    main:
        Channel.value(params.genome)
            .set{ genome_ch }

        TRUST4(fastq_file_ch, genome_ch)

        TRUST4.out.trust4_report
            .collect()
            .set { immune_results_ch }

        IMMUNE_REPORT(immune_results_ch)
    
    emit:
        clonotype_matrix = IMMUNE_REPORT.out.clonotype_matrix
        immune_summary = IMMUNE_REPORT.out.immune_summary
        immune_metadata = IMMUNE_REPORT.out.immune_metadata    
}
