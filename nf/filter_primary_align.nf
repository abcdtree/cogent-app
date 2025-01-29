nextflow.enable.dsl=2

conda_env_path = "${params.cogent_root}" + '/config/env/salmon.yaml'

process FILTER_PRIMARY_ALIGN {
    conda conda_env_path
    publishDir = [
      [
        path: { "${params.output_dir}" },
        mode: 'copy',
        pattern: "04b-primary_alignment*bam/*.bam",
        enabled: params.keep_bam
      ]
    ]

    input:
        tuple val(sample_id), path(bam)
        val (is_genome)

    output:
        tuple val(sample_id), path("04b-primary_alignment*bam/*.toTranscriptome.primary.out.bam"), emit: primary_transcript_bam, optional: true
        tuple val(sample_id), path("04b-primary_alignment*bam/*.genome.primary.out.bam"), emit: primary_genome_bam, optional: true
        path ("*counts*.txt"), emit: primary_align_counts

    script:
        def dir_name = "04b-primary_alignment_bam" 
        def output_bam_name = is_genome ? "${sample_id}_Aligned.genome.primary.out.bam" : "${sample_id}_Aligned.toTranscriptome.primary.out.bam"
        def count_file_name = is_genome ? "${sample_id}_counts_genome.txt" : "${sample_id}_counts_transcriptome.txt"

        """
        mkdir -p ${dir_name}

        samtools view -h -F 256 ${bam} \
        | tee \
            >(LINES=\$(grep -v -e "^@" | wc -l); echo \$((LINES/2)) > ${count_file_name}) \
            | (samtools view -b - > ${dir_name}/${output_bam_name})
        
        """
}
