process QUAST {
    tag "$meta.id"
    publishDir "${params.outdir}/02_assembly/quast/${meta.id}", mode: 'copy'

    input:
    tuple val(meta), path(assembly)

    output:
    path "${meta.id}_quast/*",  emit: results

    script:
    """
    quast ${assembly} -o ${meta.id}_quast --threads ${task.cpus}
    """
}
