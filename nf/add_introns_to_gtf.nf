nextflow.enable.dsl=2

intron_to_gtf = "${params.cogent_root}" + '/lib/add_genome/convert_intronic_to_exonic_gtf.R'
generate_gene_transcript_map = "${params.cogent_root}" + '/lib/add_genome/generate_gene_transcript_map.R'
conda_env_path = "${params.cogent_root}" + '/config/env/add_introns_to_gtf.yaml'

process ADD_INTRONS_TO_GTF {
    storeDir "${params.genome_dir}"
    conda conda_env_path

    input:
        path (gtf)

    output:
        path("${params.genome_name}/${params.genome_name}.intron.annotation.gtf"), emit: intron_gtf
        path("${params.genome_name}/${params.genome_name}.gene_transcript_map.txt"), emit: gene_transcript_map

    script:
        """
        mkdir -p ${params.genome_name}

        awk 'BEGIN{FS=OFS="\\t"} \$3 == "gene"' ${gtf} | convert2bed -i gtf - > temp.gene.bed
        awk 'BEGIN{FS=OFS="\\t"} \$3 == "exon"' ${gtf} | convert2bed -i gtf - > temp.exons.bed
        bedtools subtract -a temp.gene.bed -b temp.exons.bed >  temp.introns.bed

        awk '{print \$1"\\t"\$7"\\t"\$8"\\t"(\$2+1)"\\t"\$3"\\t"\$5"\\t"\$6"\\t"\$9"\\t"(substr(\$0, index(\$0,\$10)))}' temp.introns.bed \
            > temp.introns.gtf

        Rscript ${intron_to_gtf} temp.introns.gtf temp.introns_to_exons.gtf
        sed 1,3d temp.introns_to_exons.gtf > temp.introns_to_exons.final.gtf
        cat ${gtf} temp.introns_to_exons.final.gtf > ${params.genome_name}/${params.genome_name}.intron.annotation.gtf

        Rscript ${generate_gene_transcript_map} --gtf ${params.genome_name}/${params.genome_name}.intron.annotation.gtf \
            --outPath ${params.genome_name} --outFile ${params.genome_name}.gene_transcript_map.txt

        """
}
