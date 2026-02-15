process PLASMID_REPORT {
    tag "report"
    publishDir "${params.outdir}/08_report", mode: 'copy'

    input:
    path bakta_summaries
    path abricate_results
    path amrfinder_results
    path plasmidfinder_results
    path mob_suite_results
    path isfinder_results
    path integron_results

    output:
    path "plasmid_analysis_report.html",  emit: report
    path "summary_tables/*.tsv",          emit: tables
    path "figures/*.png",                 emit: figures

    script:
    """
    generate_plasmid_report.py \
        --bakta ${bakta_summaries} \
        --abricate ${abricate_results} \
        --amrfinder ${amrfinder_results} \
        --plasmidfinder ${plasmidfinder_results} \
        --mob_suite ${mob_suite_results} \
        --isfinder ${isfinder_results} \
        --integrons ${integron_results} \
        --output_dir .
    """
}
