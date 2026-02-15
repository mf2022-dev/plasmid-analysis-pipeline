process AMRFINDERPLUS {
    tag "$meta.id"
    publishDir "${params.outdir}/04_amr_virulence/amrfinderplus/${meta.id}", mode: 'copy'

    input:
    tuple val(meta), path(assembly)
    tuple val(meta), path(proteins)
    tuple val(meta), path(gff)

    output:
    path "${meta.id}_amrfinder.tsv",  emit: results

    script:
    """
    amrfinder \
        --nucleotide ${assembly} \
        --protein ${proteins} \
        --gff ${gff} \
        --organism "${params.amrfinder_organism}" \
        --plus \
        --threads ${task.cpus} \
        --output ${meta.id}_amrfinder.tsv
    """
}
