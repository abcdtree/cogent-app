nextflow.enable.dsl=2

conda_env_path = "${params.cogent_root}" + '/config/env/star.yaml'

process STAR_ALIGN {
    conda conda_env_path
    publishDir = [
      [
        path: { "${params.output_dir}" },
        mode: 'copy',
        pattern: "04-star_align_bams/*.bam",
        enabled: params.keep_bam
      ],
      [
        path: { "${params.output_dir}" },
        mode: 'copy',
        pattern: "logs/star_align/*.out",
      ],
      [
        path: { "${params.output_dir}" },
        mode: 'copy',
        pattern: "04-star_align_bams/*.out.junction",
      ],
      [
        path: { "${params.output_dir}" },
        mode: 'copy',
        pattern: "logs/star_align/*.out.tab",
      ]
    ]

    input:
        tuple val(sample_id), path(reads), path(star_index)

    output:
        tuple val(sample_id), path('04-star_align_bams/*Aligned.out.bam'), emit: bam
        tuple val(sample_id), path('04-star_align_bams/*toTranscriptome.out.bam'), emit: transcript_bam
        tuple val(sample_id), path('04-star_align_bams/*Chimeric.out.junction'), emit: junction_file
        path('logs/star_align/*Log.out'), emit: log_file
        path('logs/star_align/*Log.final.out'), emit: log_final_file
        path('logs/star_align/*Log.progress.out'), emit: log_progress_file
        path('logs/star_align/*ReadsPerGene.out.tab'), emit: read_counts_file

    script:

        def read1 = "${reads[0]}"
        def read2 = "${reads[1]}"
        """
        STAR \
            --genomeDir "${star_index}" \
            --readFilesIn $read1 $read2 \
            --outFileNamePrefix "04-star_align_bams/${sample_id}_" \
            --readFilesCommand 'zcat' \
            --runThreadN ${task.cpus} \
            --quantTranscriptomeBan ${params.quantTranscriptomeBan} \
            --quantMode ${params.quantMode} GeneCounts \
            --outSAMtype ${params.outSAMtype} \
            --genomeLoad ${params.genomeLoad} \
            --outReadsUnmapped ${params.outReadsUnmapped} \
            --outSAMstrandField ${params.outSAMstrandField} \
            --chimSegmentMin ${params.chimSegmentMin} \
            --chimJunctionOverhangMin ${params.chimJunctionOverhangMin} \
            --chimOutJunctionFormat ${params.chimOutJunctionFormat} \
            --alignSJDBoverhangMin ${params.alignSJDBoverhangMin} \
            --alignMatesGapMax ${params.alignMatesGapMax} \
            --alignIntronMax ${params.alignIntronMax} \
            --alignSJstitchMismatchNmax ${params.alignSJstitchMismatchNmax} \
            --outSAMattrRGline ${params.outSAMattrRGline} \
            --chimMultimapScoreRange ${params.chimMultimapScoreRange} \
            --chimScoreJunctionNonGTAG ${params.chimScoreJunctionNonGTAG} \
            --chimMultimapNmax ${params.chimMultimapNmax} \
            --chimNonchimScoreDropMin ${params.chimNonchimScoreDropMin} \
            --peOverlapNbasesMin ${params.peOverlapNbasesMin} \
            --peOverlapMMp ${params.peOverlapMMp} \
            --alignInsertionFlush ${params.alignInsertionFlush} \
            --alignSplicedMateMapLminOverLmate ${params.alignSplicedMateMapLminOverLmate} \
            --alignSplicedMateMapLmin ${params.alignSplicedMateMapLmin} 
            
        mkdir -p logs/star_align
        mv 04-star_align_bams/*Log*out logs/star_align
        mv 04-star_align_bams/*.out.tab logs/star_align
        """
}
