nextflow.enable.dsl=2

conda_env_path = "${params.cogent_root}" + '/config/env/salmon.yaml'

process GET_INTERGENIC_COUNTS {
    conda conda_env_path

    input:
        tuple val(sample_id), path(bam), path(intergenic_bed)

    output:
        path ("*_counts_intergenic.txt"), emit: intergenic_counts

    script:

        """

        samtools view -b ${bam} -L ${intergenic_bed} \
        | samtools view -c -F 256 > ${sample_id}_counts_intergenic.txt
                
        """
}
