#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// import relevant module
include { FASTQC } from '../modules/fastqc'
include { CUTADAPT } from '../modules/cutadapt'
include { FASTQC_TRIM } from '../modules/fastqc_trim'
include { SORTMERNA } from '../modules/sortmerna'
include { STAR_ALIGN } from '../modules/star_align'
include { STAR_LOAD_GENOME } from '../modules/star_load_genome'
include { DEDUP_UMITOOLS;
    DEDUP_UMITOOLS as DEDUP_UMITOOLS_UMIONLY } from '../modules/dedup_umitools'
include { SALMON_QUANT;
    SALMON_QUANT as SALMON_QUANT_UMIONLY } from '../modules/salmon_quant'
include { SALMON_TRANSCRIPT_MATRIX;
    SALMON_TRANSCRIPT_MATRIX as SALMON_TRANSCRIPT_MATRIX_UMIONLY } from '../modules/salmon_transcript_matrix'
include {COMPUTE_GENELEVEL_COUNTS;
    COMPUTE_GENELEVEL_COUNTS as COMPUTE_GENELEVEL_COUNTS_UMIONLY } from '../modules/compute_genelevel_counts'
include { FUSION_WORKFLOW } from './fusion_workflow'
include { IMMUNE_WORKFLOW } from './immune_workflow'
include { CREATE_ANALYZE_STATS } from '../modules/create_analyze_stats'
include { COGENTDS_REPORT } from '../modules/cogentds_report'

workflow ANALYZE_WORKFLOW {

    take:
        input_fastq_ch

    main:

        transcript_fa = file("${params.genome_dir}/${params.genome}/${params.genome}.intron.transcripts.fa")
        gene_info_file = file("${params.genome_dir}/${params.genome}/gene_info.csv")
        transcript_info_file = file("${params.genome_dir}/${params.genome}/transcript_info.csv")
        gene_transcript_map_file = file("${params.genome_dir}/${params.genome}/${params.genome}.gene_transcript_map.txt")
        star_index = file("${params.genome_dir}/${params.genome}/star_index")

        // only validate genome paths if it is not in a s3 bucket
        if ( !params.genome_dir.startsWith('s3:') ) {
            if ( !transcript_fa.exists() ) {
                exit 1, "Transcriptome file does not exist, please run add_genome first"
            }
            if ( !transcript_info_file.exists() ) {
                exit 1, "Transcript info file does not exist, please run add_genome first"
            }
            if ( !gene_info_file.exists() ) {
                exit 1, "Gene info file does not exist, please run add_genome first"
            }
            if ( !gene_transcript_map_file.exists() ) {
                exit 1, "Gene Transcript Map file does not exist, please run add_genome first"
            }
            if ( !star_index.exists() ) {
                exit 1, "STAR Index path does not exist, please run add_genome first"
            }
        }

        Channel.fromPath(transcript_fa)
            .ifEmpty { error "Transcriptome fasta file does not exist, please run add_genome first" }
            .set{ transcript_fa_ch }

        Channel.fromPath(star_index)
            .ifEmpty { error "STAR Index path does not exist, please run add_genome first" }
            .set{ star_index_ch }

        Channel.fromPath(gene_info_file)
            .ifEmpty { error "Gene info file does not exist, please run add_genome first" }
            .set{ gene_info_file_ch }

        Channel.fromPath(transcript_info_file)
            .ifEmpty { error "Transcript info file does not exist, please run add_genome first" }
            .set{ transcript_info_file_ch }

        Channel.fromPath(gene_transcript_map_file)
            .ifEmpty { error "Gene Transcript Map file does not exist, please run add_genome first" }
            .set{ gene_transcript_map_file_ch }


        CUTADAPT(input_fastq_ch)

        cutadapt_log_ch = Channel.empty()
        CUTADAPT.out.json
            .collect()
            .set{ cutadapt_log_ch }

        STAR_LOAD_GENOME(star_index_ch, cutadapt_log_ch)

        def fastq_trimmed = CUTADAPT.out.trimmed_reads

        fastq_trimmed
            .collect()
            .flatten()
            .map{ file ->
                    def sample = file.name.replaceAll("_R[12].fastq.gz", "")
                    tuple(sample, file)
            }
            .groupTuple()
            .set{ fastq_trimmed_ch }

        if ( params.fastqc ) {
            FASTQC(input_fastq_ch)
            FASTQC_TRIM(fastq_trimmed_ch)
        }

        // check if rRNA removal is needed

        sortmerna_log_ch = Channel.empty()
        clonotype_matrix_ch = Channel.empty()
        clonotype_summary_ch = Channel.empty()
        clonotype_metadata_ch = Channel.empty()

        if ( !params.skip_ribodepletion ) {

            ribo_db = file("${params.genome_dir}/${params.genome}/sortmernadb/idx", checkIfExists: true)

            sortmerna_index_ch = Channel.fromPath(ribo_db)

            Channel.fromPath("${params.genome_dir}/${params.genome}/sortmernadb/ref/*")
                .ifEmpty { error "Could not find ribosomal database fasta files!" }
                .collect()
                .set { sortmerna_fastas_ch }

            fastq_trimmed_ch
                .combine(sortmerna_index_ch)
                .set{ sortmerna_input_ch }

            SORTMERNA(sortmerna_input_ch, sortmerna_fastas_ch)

            SORTMERNA.out.ribodepleted_fastqs
                .combine(star_index_ch)
                .set{ star_align_input_ch }

            SORTMERNA.out.log
                .collect()
                .set{ sortmerna_log_ch }


            if ( params.immune ){

                IMMUNE_WORKFLOW(SORTMERNA.out.ribodepleted_fastqs)

            }

        } else {

            fastq_trimmed_ch
                .combine(star_index_ch)
                .set{ star_align_input_ch }

            if ( params.immune ){

                IMMUNE_WORKFLOW(fastq_trimmed_ch)
            }

        }

        if ( params.immune ){
            clonotype_matrix_ch = IMMUNE_WORKFLOW.out.clonotype_matrix
            clonotype_summary_ch = IMMUNE_WORKFLOW.out.immune_summary
            clonotype_metadata_ch = IMMUNE_WORKFLOW.out.immune_metadata
        }

        star_log_ch = Channel.empty()
        star_reads_ch = Channel.empty()
        umitools_log_ch = Channel.empty()
        salmon_log_ch = Channel.empty()
        salmon_matrix_ch = Channel.empty()
        salmon_transcript_matrix_ch = Channel.empty()
        salmon_intron_matrix_ch = Channel.empty()
        junction_matrix_ch = Channel.empty()
        spanning_matrix_ch = Channel.empty()


        def umi_only = false

        STAR_ALIGN(star_align_input_ch)

        if (params.fusion){
                FUSION_WORKFLOW(STAR_ALIGN.out.junction_file)

                junction_matrix_ch = FUSION_WORKFLOW.out.junction_matrix
                spanning_matrix_ch = FUSION_WORKFLOW.out.spanning_matrix
        }

        STAR_ALIGN.out.log_final_file
            .collect()
            .set{ star_log_ch }

        STAR_ALIGN.out.read_counts_file
            .collect()
            .set{ star_reads_ch }

        if (params.type_of_experiment == 'stranded_umi'){

            DEDUP_UMITOOLS(STAR_ALIGN.out.transcript_bam)

            DEDUP_UMITOOLS.out.dedup_bam
                .combine(transcript_fa_ch)
                .set { salmon_input_ch }

            DEDUP_UMITOOLS.out.umitools_log
                .collect()
                .set{ umitools_log_ch }

        } else {
            STAR_ALIGN.out.transcript_bam
                .combine(transcript_fa_ch)
                .set { salmon_input_ch }

        }

        SALMON_QUANT(salmon_input_ch, umi_only)

        SALMON_QUANT.out.salmon_transcript_results
            .collect()
            .set{ salmon_transcript_results_ch }

        SALMON_QUANT.out.salmon_logs
            .collect()
            .set{ salmon_log_ch }

        SALMON_TRANSCRIPT_MATRIX(salmon_transcript_results_ch, umi_only)

        COMPUTE_GENELEVEL_COUNTS(
            SALMON_TRANSCRIPT_MATRIX.out.transcript_matrix,
            gene_transcript_map_file_ch,
            gene_info_file_ch,
            transcript_info_file_ch,
            umi_only
        )

        salmon_transcript_matrix_ch = COMPUTE_GENELEVEL_COUNTS.out.transcript_matrix
        salmon_matrix_ch = COMPUTE_GENELEVEL_COUNTS.out.gene_matrix
        salmon_intron_matrix_ch = COMPUTE_GENELEVEL_COUNTS.out.intron_gene_matrix


        umitools_umionly_log_ch = Channel.empty()
        salmon_umionly_matrix_ch = Channel.empty()
        salmon_umionly_transcript_matrix_ch = Channel.empty()
        salmon_intron_matrix_umionly_ch = Channel.empty()

        if ( params.type_of_experiment == 'smartseq_fla_umi'){
            umi_only = true

            DEDUP_UMITOOLS_UMIONLY(STAR_ALIGN.out.transcript_bam)

            DEDUP_UMITOOLS_UMIONLY.out.dedup_bam
                .combine(transcript_fa_ch)
                .set { salmon_umionly_input_ch }

            DEDUP_UMITOOLS_UMIONLY.out.umitools_log
                .collect()
                .set{ umitools_umionly_log_ch }

            SALMON_QUANT_UMIONLY(salmon_umionly_input_ch, umi_only)

            SALMON_QUANT_UMIONLY.out.salmon_transcript_results
                .collect()
                .set{ salmon_umionly_transcript_results_ch }

            SALMON_QUANT_UMIONLY.out.salmon_logs
                .collect()
                .set{ salmon_umionly_log_ch }

            SALMON_TRANSCRIPT_MATRIX_UMIONLY(salmon_umionly_transcript_results_ch, umi_only)

            COMPUTE_GENELEVEL_COUNTS_UMIONLY(
                SALMON_TRANSCRIPT_MATRIX_UMIONLY.out.transcript_matrix,
                gene_transcript_map_file_ch,
                gene_info_file_ch,
                transcript_info_file_ch,
                umi_only
            )

            salmon_umionly_transcript_matrix_ch = COMPUTE_GENELEVEL_COUNTS_UMIONLY.out.transcript_matrix
            salmon_umionly_matrix_ch = COMPUTE_GENELEVEL_COUNTS_UMIONLY.out.gene_matrix
            salmon_intron_matrix_umionly_ch = COMPUTE_GENELEVEL_COUNTS_UMIONLY.out.intron_gene_matrix
        }

        mito_file_ch = Channel.empty()

        Channel.fromPath("${params.genome_dir}/${params.genome}/mito_reads/*.txt")
                .ifEmpty { error "Could not find mito read files!" }
                .collect()
                .set { mito_file_ch }

        CREATE_ANALYZE_STATS(
            cutadapt_log_ch.ifEmpty([]),
            sortmerna_log_ch.ifEmpty([]),
            salmon_log_ch.ifEmpty([]),
            star_log_ch.ifEmpty([]),
            star_reads_ch.ifEmpty([]),
            mito_file_ch.ifEmpty([]),
            salmon_intron_matrix_ch.ifEmpty([]),
            salmon_matrix_ch.ifEmpty([]),
            salmon_transcript_matrix_ch.ifEmpty([]),
            umitools_log_ch.ifEmpty([]),
            umitools_umionly_log_ch.ifEmpty([]),
            salmon_umionly_matrix_ch.ifEmpty([]),
            salmon_umionly_transcript_matrix_ch.ifEmpty([]),
            salmon_intron_matrix_umionly_ch.ifEmpty([]),
        )

        COGENTDS_REPORT(
            salmon_matrix_ch.ifEmpty([]),
            salmon_transcript_matrix_ch.ifEmpty([]),
            gene_info_file_ch,
            transcript_info_file_ch,
            clonotype_matrix_ch.ifEmpty([]),
            clonotype_summary_ch.ifEmpty([]),
            clonotype_metadata_ch.ifEmpty([]),
            junction_matrix_ch.ifEmpty([]),
            spanning_matrix_ch.ifEmpty([]),
            CREATE_ANALYZE_STATS.out.stats_file
        )
}
