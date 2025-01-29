#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { GENOME_DOWNLOAD } from '../modules/genome_download'
include { ANNOTATION_DOWNLOAD } from '../modules/annotation_download'
include { ADD_INTRONS_TO_GTF } from '../modules/add_introns_to_gtf'
include { GFFREAD_INTRON } from '../modules/gffread_intron'
include { STAR_INDEX } from '../modules/star_index'
include { STAR_FUSION_INDEX } from '../modules/star_fusion_index'
include { SORTMERNA_DOWNLOAD } from '../modules/sortmerna_download'
include { SORTMERNA_INDEX } from '../modules/sortmerna_index'
include { MITOCHONDRIAL_GENE_LIST } from '../modules/mitochondrial_gene_list'
include { CREATE_GENE_INFO } from '../modules/create_gene_info'
include { CREATE_TRANSCRIPT_INFO } from '../modules/create_transcript_info'


workflow ADD_GENOME_WORKFLOW {

    if(!params.genome_fasta_url) {
        error("no genome fasta URL specified --genome_fasta_url URL")
    }

    if(!params.annotation_gtf_url) {
        error("no gtf URL specified --annotation_gtf_url URL")
    }

    GENOME_DOWNLOAD(params.genome_fasta_url)
    ANNOTATION_DOWNLOAD(params.annotation_gtf_url)

    CREATE_GENE_INFO(ANNOTATION_DOWNLOAD.out.gtf)
    CREATE_TRANSCRIPT_INFO(ANNOTATION_DOWNLOAD.out.gtf)

    ADD_INTRONS_TO_GTF(ANNOTATION_DOWNLOAD.out.gtf)

    if (params.fusion) {
        STAR_FUSION_INDEX()
    }

    GFFREAD_INTRON(
        GENOME_DOWNLOAD.out.genome_fa,
        ADD_INTRONS_TO_GTF.out.intron_gtf)

    STAR_INDEX(
        GENOME_DOWNLOAD.out.genome_fa,
        ADD_INTRONS_TO_GTF.out.intron_gtf)
        
    Channel.fromPath("${params.sortmerna_fastas}")
        .ifEmpty { error "Could not find ribosomal database fasta files!" }
        .collect()
        .set { sortmerna_fastas_ch }

    SORTMERNA_INDEX(sortmerna_fastas_ch)

    Channel.fromPath("${params.mito_reads}")
        .ifEmpty { error "Could not find mitochondrial gene name file!" }
        .collect()
        .set { mito_reads_ch }

    MITOCHONDRIAL_GENE_LIST(mito_reads_ch)

}
