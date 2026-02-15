process MOBILE_ELEMENT_FINDER {
    tag "$meta.id"
    publishDir "${params.outdir}/06_mge_analysis/mobile_element_finder/${meta.id}", mode: 'copy'

    input:
    tuple val(meta), path(assembly)

    output:
    path "${meta.id}_mge_results/*",  emit: results

    script:
    """
    mkdir -p ${meta.id}_mge_results
    mefinder find \
        --contig ${assembly} \
        --outdir ${meta.id}_mge_results \
    || echo "No mobile elements found" > ${meta.id}_mge_results/summary.txt
    """
}
