gene_list_script = "${params.cogent_root}" + '/lib/salmon/create_gene_list.py'
conda_env_path = "${params.cogent_root}" + '/config/env/salmon_report.yaml'

process CREATE_TRANSCRIPT_INFO {
    conda conda_env_path
    storeDir "${params.genome_dir}"

    input:
        path (gtf_file)

    output:
        path("${params.genome_name}/transcript_info.csv"), emit: transcript_info

    script:
        """
        python3 ${gene_list_script} --type 'transcripts' -r ${params.genome_name} -o transcript_info.csv -g ${gtf_file}
        """
}
