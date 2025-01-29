#!/usr/bin/env nextflow

/*
Cogent Analysis Pipeline
TBUSA BI team | Maintainer: Karthik Ram Padmanabhan
v3.0.18
*/

nextflow.enable.dsl = 2
nextflow.enable.strict = true

/*
Configurations for pipeline parameters need to be set
in the nextflow.config file
*/

log.info """\
        Name: ${params.manifest.name}
        Version: v${params.manifest.version}
        Analysis type: ${params.analysis}
        --
        run as: ${workflow.commandLine}
        started at: ${workflow.start}
        config files: ${workflow.configFiles}
        """
        .stripIndent()

// import workflows
include { DEMUX_WORKFLOW } from './workflows/demux_workflow'
include { ANALYZE_WORKFLOW } from './workflows/analyze_workflow'
include { ADD_GENOME_WORKFLOW } from './workflows/add_genome_workflow'
include { FUSION_WORKFLOW } from './workflows/fusion_workflow'
include { IMMUNE_WORKFLOW } from './workflows/immune_workflow'

/*
 * main script flow
 */

workflow {
    if ( params.analysis == "demux" ){

        DEMUX_WORKFLOW ()

    } else if( params.analysis == "analyze" ){

        fastqs = "${params.input_dir}/*_R{1,2}*.f*q*"

        Channel.fromFilePairs( fastqs )
            .ifEmpty { error "Could not find any reads matching pattern: ${fastqs}" }
            .set { input_fastq_ch }

        ANALYZE_WORKFLOW (input_fastq_ch)

    } else if( params.analysis == "demux_and_analyze") {

        if ( params.demux_mode == 'move' ) {
            { error "demux_mode is set to 'move', change it to 'copy' using --demux_mode before running this analysis type" }
        }

        DEMUX_WORKFLOW ()

        def fastq_reads = DEMUX_WORKFLOW.out.demux_result_ch

        fastq_reads
            .collect()
            .flatten()
            .map{ file ->
                    def sample = file.name.replaceAll("_R[12].fastq.gz", "")
                    tuple(sample, file)
            }
            .groupTuple()
            .ifEmpty { error "Could not find any reads matching pattern: *_R{1,2}.f*q.gz" }
            .set{ input_fastq_ch }

        ANALYZE_WORKFLOW (input_fastq_ch)

    } else if( params.analysis == "add_genome" ){

        ADD_GENOME_WORKFLOW ()

    } else if( params.analysis == "run_fusion_analysis" ) {

        def star_junction_files = "${params.input_dir}/04-star_align_bams/*Chimeric.out.junction"

        Channel.fromPath(star_junction_files)
            .collect()
            .flatten()
            .map{ file ->
                def sample = file.name.replaceAll("_Chimeric.out.junction", "")
                tuple(sample, file)
            }
            .groupTuple()
            .ifEmpty { error "Cannot find junction files in input directory, please run cogent analyze first!" }
            .set { junction_file_ch }

        FUSION_WORKFLOW(junction_file_ch)


    } else if( params.analysis == "run_immune_analysis" ) {

        def trimmed_out = "${params.input_dir}/02-cutadapt_trimmed_fastqs/*_R{1,2}.fastq.gz"
        def ribodepletion_out = "${params.input_dir}/03a-ribodepletion/*_non_rRNA_R{1,2}.fastq.gz"

        Channel.fromFilePairs( trimmed_out )
            .ifEmpty { error "Cannot find fastqs in input directory, please run cogent analyze first!" }
            .set { trimmed_fastq_ch }

        Channel.fromFilePairs( ribodepletion_out )
            .set { ribodepleted_fastq_ch }

        ribodepleted_fastq_path = file("${params.input_dir}/03a-ribodepletion")

        if ( ribodepleted_fastq_path.exists() ) {
            IMMUNE_WORKFLOW (ribodepleted_fastq_ch)
        } else {
            IMMUNE_WORKFLOW (trimmed_fastq_ch)
        }

    } else {
        exit 1, "Invalid analysis type, please check documentation for compatible analyses."
    }
}
