nextflow.enable.dsl=2

conda_env_path = "${params.cogent_root}" + '/config/env/umi_tools.yaml'

process DEDUP_UMITOOLS {
    conda conda_env_path
    publishDir = [
      [
        path: { "${params.output_dir}" },
        mode: 'copy',
        pattern: "04a-dedup_bam/*dedup.out.bam",
        enabled: params.keep_bam
      ],
      [
        path: { "${params.output_dir}" },
        mode: 'copy',
        pattern: "logs/umi_tools/",
      ]
    ]

    input:
        tuple val(sample_id), path(bam)

    output:
        tuple val (sample_id), path("04a-dedup_bam/*dedup.out.bam"), emit: dedup_bam
        path("logs/umi_tools/*umitools.log"), emit: umitools_log

    script:


        def tag_file = 'UM'
        def modified_bam = 'temp_with_UM.bam'
        def paired_flag = '--paired'
        def unpaired_flag = '--unpaired-reads discard'

        def perl_cmd = "perl -pe 's/^([^\\t]+)_([^\\t]+)(.+)/\\1_\\2\\3\\tUM:Z:\\2/'"

        if ( params.type_of_experiment == "smartseq_fla_umi"){
            perl_cmd = "perl -ne 's/^([^\\t]+)_([^\\t]+)(.+)/\\1_\\2\\3\\tUM:Z:\\2/; print if \$2 ne 'NNNNNNNN''"
        }

        """
        mkdir -p 04a-dedup_bam
        mkdir -p logs/umi_tools

        samtools view -h ${bam} | ${perl_cmd} | samtools sort | samtools view -b > ${modified_bam}

        samtools index ${modified_bam}

        umi_tools dedup --stdin=${modified_bam} --log=logs/umi_tools/${sample_id}_umitools.log \
            ${paired_flag} --no-sort-output --extract-umi-method=tag ${unpaired_flag} \
            --umi-tag=UM | samtools sort -n > 04a-dedup_bam/${sample_id}_Aligned.toTranscriptome.dedup.out.bam

        """
}
