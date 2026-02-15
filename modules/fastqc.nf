process FASTQC {
    tag "$meta.id"
    publishDir "${params.outdir}/01_qc/fastqc", mode: 'copy'

    input:
    tuple val(meta), path(reads1), path(reads2)

    output:
    path "*.html",  emit: html
    path "*.zip",   emit: zip

    script:
    """
    fastqc --threads ${task.cpus} --outdir . ${reads1} ${reads2}
    """
}
