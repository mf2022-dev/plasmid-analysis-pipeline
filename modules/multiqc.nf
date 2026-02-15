process MULTIQC {
    tag "multiqc"
    publishDir "${params.outdir}/08_report", mode: 'copy'

    input:
    path fastqc_results
    path fastp_results
    path quast_results

    output:
    path "multiqc_report.html",  emit: report
    path "multiqc_data/*",       emit: data

    script:
    def config = params.multiqc_config ? "--config ${params.multiqc_config}" : ""
    """
    multiqc . ${config} --outdir . --filename multiqc_report
    """
}
