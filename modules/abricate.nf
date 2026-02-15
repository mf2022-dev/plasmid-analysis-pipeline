process ABRICATE {
    tag "$meta.id"
    publishDir "${params.outdir}/04_amr_virulence/abricate/${meta.id}", mode: 'copy'

    input:
    tuple val(meta), path(assembly)

    output:
    path "${meta.id}_card.tsv",            emit: card
    path "${meta.id}_vfdb.tsv",            emit: vfdb
    path "${meta.id}_plasmidfinder.tsv",   emit: plasmidfinder
    path "${meta.id}_abricate_*.tsv",      emit: results

    script:
    """
    abricate --db card           ${assembly} > ${meta.id}_abricate_card.tsv
    abricate --db vfdb           ${assembly} > ${meta.id}_abricate_vfdb.tsv
    abricate --db plasmidfinder  ${assembly} > ${meta.id}_abricate_plasmidfinder.tsv

    cp ${meta.id}_abricate_card.tsv           ${meta.id}_card.tsv
    cp ${meta.id}_abricate_vfdb.tsv           ${meta.id}_vfdb.tsv
    cp ${meta.id}_abricate_plasmidfinder.tsv  ${meta.id}_plasmidfinder.tsv
    """
}
