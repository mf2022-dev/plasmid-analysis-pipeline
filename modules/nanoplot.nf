process NANOPLOT {
    tag "$meta.id"
    publishDir "${params.outdir}/01_qc/nanoplot/${meta.id}", mode: 'copy'

    input:
    tuple val(meta), path(long_reads)

    output:
    path "*.html",           emit: html
    path "NanoStats.txt",    emit: stats
    path "*.png",            emit: plots

    script:
    """
    NanoPlot --fastq ${long_reads} --outdir . --threads ${task.cpus} --plots dot
    """
}
