process CLINKER {
    tag "comparative"
    publishDir "${params.outdir}/07_comparative/clinker", mode: 'copy'

    input:
    path gbk_files

    output:
    path "clinker_output.html",  emit: html
    path "clinker_output.csv",   emit: csv

    script:
    """
    clinker ${gbk_files} \
        -p clinker_output.html \
        --identity 0.3 \
        -o clinker_output.csv
    """
}
