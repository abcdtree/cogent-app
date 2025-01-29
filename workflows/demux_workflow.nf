#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// import demux module
include { FASTQC_DEMUX } from '../modules/fastqc_demux'
include { DEMUX } from '../modules/demux'


// set optional demux params

workflow DEMUX_WORKFLOW {

    main:

        Channel.fromPath(params.barcodes_file)
            .set{ bcs_file_ch }

        Channel.fromPath(params.fastq1)
            .set{ fastq1_ch }

        Channel.fromPath(params.fastq2)
            .set{ fastq2_ch }

        if ( params.fastqc ) {
            FASTQC_DEMUX(fastq1_ch, fastq2_ch)
        }

        DEMUX(fastq1_ch, fastq2_ch, bcs_file_ch)

    emit:
        demux_result_ch = DEMUX.out.fastqs

}
