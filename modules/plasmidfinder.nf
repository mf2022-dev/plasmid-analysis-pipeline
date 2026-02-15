process PLASMIDFINDER {
    tag "$meta.id"
    publishDir "${params.outdir}/05_plasmid_typing/plasmidfinder/${meta.id}", mode: 'copy'

    input:
    tuple val(meta), path(assembly)

    output:
    path "${meta.id}_plasmidfinder_results/*",  emit: results

    script:
    """
    mkdir -p ${meta.id}_plasmidfinder_results
    plasmidfinder.py \
        -i ${assembly} \
        -o ${meta.id}_plasmidfinder_results \
        -p ${params.plasmidfinder_db} \
        -l ${params.plasmidfinder_cov} \
        -t ${params.plasmidfinder_id}
    """
}
