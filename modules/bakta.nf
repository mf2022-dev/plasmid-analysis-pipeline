process BAKTA {
    tag "$meta.id"
    publishDir "${params.outdir}/03_annotation/bakta/${meta.id}", mode: 'copy'

    input:
    tuple val(meta), path(assembly)

    output:
    tuple val(meta), path("${meta.id}.gbk"),    emit: gbk
    tuple val(meta), path("${meta.id}.gff3"),   emit: gff
    tuple val(meta), path("${meta.id}.faa"),    emit: proteins
    path "${meta.id}.tsv",                       emit: summary
    path "${meta.id}.txt",                       emit: txt

    script:
    def db_arg = params.bakta_db ? "--db ${params.bakta_db}" : ""
    """
    bakta ${db_arg} \
        --output . \
        --prefix ${meta.id} \
        --genus ${params.genus} --species "${params.species}" \
        --strain "${meta.id}" \
        --threads ${task.cpus} \
        ${assembly}
    """
}
