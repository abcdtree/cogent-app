conda_env_path = "${params.cogent_root}" + '/config/env/cutadapt.yaml'
adapters = "${params.cogent_root}" + '/config/adapters.fasta'

process CUTADAPT {
    conda conda_env_path
    publishDir "${params.output_dir}", mode: 'copy'

    input:
        tuple val(sample_id), path(reads)

    output:
        path("02-cutadapt_trimmed_fastqs/${sample_id}_trimmed_R*.fastq.gz"), emit: trimmed_reads
        path("logs/cutadapt/${sample_id}_cutadapt_log.txt"), emit: log
        path("logs/cutadapt/${sample_id}.cutadapt.json"), emit: json
        //path("logs"), emit: log
        path("00-gzipped_fastqs/${sample_id}_R*.fastq.gz"), optional:true

    script:

        def read1 = "${reads[0]}"
        def read2 = "${reads[1]}"

        def outr1 = "-o 02-cutadapt_trimmed_fastqs/${sample_id}_trimmed_R1.fastq.gz"
        def outr2 = "-p 02-cutadapt_trimmed_fastqs/${sample_id}_trimmed_R2.fastq.gz"

        def adapter_r1 = "-a file:${adapters}"
        def adapter_r2 = "-A file:${adapters}"

        def min_read_length = 20
        def n_bases_removed = 0.7
        def quality_score = 20

        def gzip_read1 = ''
        def gzip_read2 = ''

        if (! read1.endsWith('.gz')) {
            gzip_read1 = "mkdir -p 00-gzipped_fastqs && gzip -c ${read1} > 00-gzipped_fastqs/${sample_id}_R1.fastq.gz"
            read1 = "00-gzipped_fastqs/${sample_id}_R1.fastq.gz"
        }
        if (! read2.endsWith('.gz')) {
            gzip_read2 = "mkdir -p 00-gzipped_fastqs && gzip -c ${read2} > 00-gzipped_fastqs/${sample_id}_R2.fastq.gz"
            read2 = "00-gzipped_fastqs/${sample_id}_R2.fastq.gz"
        }


        """
        ${gzip_read1}
        ${gzip_read2}
        mkdir 02-cutadapt_trimmed_fastqs
        mkdir -p logs/cutadapt
        cutadapt -j ${task.cpus} -m ${min_read_length} --trim-n --max-n ${n_bases_removed} \
        -q ${quality_score} --json=logs/cutadapt/${sample_id}.cutadapt.json ${adapter_r1} ${adapter_r2} ${outr1} ${outr2} \
        ${read1} ${read2} > logs/cutadapt/${sample_id}_cutadapt_log.txt
        """
}
