nextflow.enable.dsl=2

cogentds_script = "${params.cogent_root}" + '/lib/report/CogentDS_DataHandler_and_Reporter.R'
cogentds_reporter = "${params.cogent_root}" + '/lib/report/'
conda_env_path = "${params.cogent_root}" + '/config/env/cogentds_report.yaml'

process COGENTDS_REPORT {
    publishDir "${params.output_dir}", mode: 'move'
    conda conda_env_path 

    input:
        path (salmon_gene_matrix)
        path (salmon_transcript_matrix)
        path (gene_info)
        path (transcript_info)
        path (clonotype_counts)
        path (clonotype_summary)
        path (clonotype_metadata)
        path (junction_matrix)
        path (spanning_matrix)
        path (analyze_stats)


    output:
        path "report/CogentDS_analysis.rds", emit: gene_object
        path "report/CogentDS_preliminary-analysis_report.html", emit: html_report
        path "report/CogentDS_isoform_only_analysis.rds", emit: isoform_object, optional: true

    script:

        def fusion_command = ''
        def immune_command = ''
        def isoforms_report_cmd = '--analyzeIsoformforReport FALSE'
        
        if (params.fusion){
            fusion_command = "--fusionjunctionCounts ${junction_matrix} --fusionspanFragments ${spanning_matrix}"
        }

        if (params.immune){
            immune_command = "--clonotypeCounts ${clonotype_counts} --clonotypeSummary ${clonotype_summary} --clonotypeMetadata ${clonotype_metadata}"
        }

        if (params.isoform_report){
            isoforms_report_cmd = '--analyzeIsoformforReport TRUE'
        }

        
        """
        export R_LIBS="./lib"
        
        Rscript ${cogentds_script} --pathtoCogentDSreporter ${cogentds_reporter} --analyzeStats ${analyze_stats} \
        --geneCounts ${salmon_gene_matrix} --isoformCounts ${salmon_transcript_matrix} \
        --geneInfo ${gene_info} --isoformInfo ${transcript_info} --outPath \$PWD/report \
        ${isoforms_report_cmd} \
        ${fusion_command} \
        ${immune_command} 
        """
}
