nextflow.enable.dsl=2

conda_env_path = "${params.cogent_root}" + '/config/env/sortmerna.yaml'

process SORTMERNA {
    conda conda_env_path
    publishDir "${params.output_dir}", mode: 'copy'

    input:
        tuple val(sample_id), path(reads), path(ribo_db)
        path(fastas)

    output:
        tuple val(sample_id), path("03a-ribodepletion/*_non_rRNA_R*.fastq.gz"), emit: ribodepleted_fastqs
        path("logs/sortmerna/*.sortmerna.log"), emit: log
        path('logs/sortmerna/*_sortmerna.stdout'), emit: stdout

    script:

        def read1 = "--reads ${reads[0]}"
        def read2 = "--reads ${reads[1]} --paired_in --out2"

        def mv_r1 = "mv 03a-ribodepletion/non_rRNA_reads_fwd.f*q.gz 03a-ribodepletion/${sample_id}_non_rRNA_R1.fastq.gz"
        def mv_r2 = "mv 03a-ribodepletion/non_rRNA_reads_rev.f*q.gz 03a-ribodepletion/${sample_id}_non_rRNA_R2.fastq.gz"

        """
        mkdir 03a-ribodepletion
        mkdir -p logs/sortmerna

        (sortmerna \\
            ${'--ref '+fastas.join(' --ref ')} \\
            --idx-dir ${ribo_db} \\
            ${read1} \\
            ${read2} \\
            --threads ${task.cpus} \\
            --workdir 03a-ribodepletion/ \\
            --kvdb 03a-ribodepletion/kvdb \\
            --aligned 03a-ribodepletion/rRNA_reads \\
            --other 03a-ribodepletion/non_rRNA_reads \\
            --fastx \\
        ) 2>&1 | tee logs/sortmerna/${sample_id}_sortmerna.stdout

        ${mv_r1}
        ${mv_r2}
        mv 03a-ribodepletion/rRNA_reads.log logs/sortmerna/${sample_id}.sortmerna.log

        """
}
