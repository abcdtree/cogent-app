nextflow.enable.dsl=2

conda_env_path = "${params.cogent_root}" + '/config/env/star.yaml'

// This extra LoadGenome step is not strictly required since STAR
// align uses --genomeLoad LoadAndRemove. However, if STAR align
// is interrupted, the shared-memory genome can end up in a bad
// state and cause STAR align to freeze indefinitely on the next run.
// Explicitly running LoadAndExit in this step will clean up a
// corrupted genome so that STAR align does not hang.

process STAR_LOAD_GENOME {
    conda conda_env_path

    input:
        path (star_index)
        path (log)

    script:

        """
        STAR \
            --genomeDir "${star_index}" \
            --genomeLoad LoadAndExit \
        """
}