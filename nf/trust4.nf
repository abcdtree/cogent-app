nextflow.enable.dsl=2

conda_env_path = "${params.cogent_root}" + '/config/env/trust4.yaml'

process TRUST4 {
    conda conda_env_path
    publishDir "${params.output_dir}", mode: 'copy'

    input:
        tuple val(sample_id), path(reads)
        val(genome_ch)

    output:
        path('10-immune_profiling/*_report.tsv'), emit: trust4_report

    script:
        def read1 = "${reads[0]}"
        def read2 = "${reads[1]}"

        def hg38_bcr_path = "${params.cogent_root}" + '/lib/immune/hg38_bcrtcr.fa'
        def mm39_bcr_path = "${params.cogent_root}" + '/lib/immune/GRCm38_bcrtcr.fa'
        def hg38_imgtc_path = "${params.cogent_root}" + '/lib/immune/human_IMGT+C.fa'
        def mm39_imgtc_path = "${params.cogent_root}" + '/lib/immune/mouse_IMGT+C.fa'

        if (genome_ch == 'hg38')

            """
            run-trust4 -t ${task.cpus} -1 $read1 -2 $read2 -f ${hg38_bcr_path} --ref ${hg38_imgtc_path} \\
            -o ${sample_id} --od '10-immune_profiling'
            """
        else if(genome_ch == 'mm39')

            """
            run-trust4 -t ${task.cpus} -1 $read1 -2 $read2 -f ${mm39_bcr_path} --ref ${mm39_imgtc_path} \\
            -o ${sample_id} --od '10-immune_profiling'
            """
        else
            error "Invalid genome name: $genome_ch, must be one of hg38 or mm39"
}

