nextflow.enable.dsl=2

cogent_path = "${params.cogent_root}"  + '/lib/cogent.py'
conda_env_path = "${params.cogent_root}" + '/config/env/demux_env.yaml'

process DEMUX {
    conda conda_env_path
    stageOutMode 'rsync'
    publishDir "${params.output_dir}", mode: "${params.demux_mode}"

    input:
        path read1
        path read2
        path bcs_file

    output:
        path("**.f*q*"), emit: fastqs, optional:true
        path("**.log"), emit: log
        path("**.csv"), emit: stats

    script:

        def no_gz_flag = ''
        if ( params.no_gz ) {
            no_gz_flag = ' --no_gz'
        }

        def undetermined_fq_flag = ''
        if ( params.undetermined_fq ){
            undetermined_fq_flag = ' --undetermined_fq'
        }

        def hpc_flag = ''
        if ( params.hpc ){
            hpc_flag = ' --hpc'
        }

        def rand_pick_flag = ''
        if ( params.random_pick ){
            rand_pick_flag = '--random_pick'
        }

        def dry_run_flag = ''
        if ( params.dry_run ){
            dry_run_flag =  '--dry_run'
        }

        def umi_length_flag = ''
        if ( params.umi_length ){
            umi_length_flag = "--umi_length ${params.umi_length}"
        }
        
        def bcs_map_flag = ''
        if ( params.bcs_map ){
            bcs_map_flag = "--bcs_map ${params.bcs_map}"
        }

        """
        python3 ${cogent_path} demux -i ${read1} -p ${read2} -b ${bcs_file} -t ${params.type_of_experiment} -o 00b-demultiplexed_fastqs ${no_gz_flag} \
        -m ${params.mismatch} ${umi_length_flag} -n ${task.cpus} --n_writers ${params.n_writers} \
        --i7_rc ${params.i7_rc} --i5_rc ${params.i5_rc} --read_buffer ${params.read_buffer}  --prog ${params.prog} ${bcs_map_flag} \
        ${undetermined_fq_flag} ${hpc_flag} ${rand_pick_flag} ${dry_run_flag} --use_barcodes ${params.use_barcodes} \
        --check_reads ${params.check_reads} --min_reads ${params.min_reads}

        """
}

