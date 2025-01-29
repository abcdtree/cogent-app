transcript_list_script = "${params.cogent_root}" + '/lib/salmon/create_gene_list.py'
conda_env_path = "${params.cogent_root}" + '/config/env/salmon_report.yaml'

process CREATE_GENE_INFO {
    conda conda_env_path
    storeDir "${params.genome_dir}"

    input:
        path (gtf_file)

    output:
        path("${params.genome_name}/gene_info.csv"), emit: gene_info

    script:
        """
        python3 ${transcript_list_script} --type 'genes' -r ${params.genome_name} -o gene_info.csv -g ${gtf_file}
        """
}
