process INTEGRON_FINDER {
    tag "$meta.id"
    publishDir "${params.outdir}/06_mge_analysis/integron_finder/${meta.id}", mode: 'copy'

    input:
    tuple val(meta), path(assembly)

    output:
    path "Results_Integron_Finder_${meta.id}/*",  emit: results

    script:
    """
    integron_finder ${assembly} \
        --local-max \
        --func-annot \
        --outdir Results_Integron_Finder_${meta.id} \
        --cpu ${task.cpus}
    """
}
