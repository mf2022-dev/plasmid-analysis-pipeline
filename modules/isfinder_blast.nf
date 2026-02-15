process ISFINDER_BLAST {
    tag "$meta.id"
    publishDir "${params.outdir}/06_mge_analysis/isfinder/${meta.id}", mode: 'copy'

    input:
    tuple val(meta), path(assembly)

    output:
    path "${meta.id}_is_elements.tsv",  emit: results

    script:
    """
    blastn \
        -query ${assembly} \
        -db \${ISFINDER_DB:-/databases/isfinder/ISfinder_db} \
        -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle" \
        -evalue ${params.isfinder_evalue} \
        -perc_identity ${params.isfinder_identity} \
        -num_threads ${task.cpus} \
        -out ${meta.id}_is_elements.tsv \
    || echo "No IS elements found" > ${meta.id}_is_elements.tsv
    """
}
