process FILTLONG {
    tag "$meta.id"
    publishDir "${params.outdir}/01_qc/filtlong", mode: 'copy'

    input:
    tuple val(meta), path(long_reads)

    output:
    tuple val(meta), path("${meta.id}_filtered.fastq.gz"), emit: reads

    script:
    """
    filtlong \
        --min_length 1000 \
        --keep_percent 95 \
        ${long_reads} | gzip > ${meta.id}_filtered.fastq.gz
    """
}
