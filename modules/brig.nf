process BRIG {
    tag "comparative"
    publishDir "${params.outdir}/07_comparative/brig", mode: 'copy'

    input:
    path reference
    path plasmids

    output:
    path "brig_output/*",  emit: results

    script:
    """
    mkdir -p brig_output query_dir
    cp ${plasmids} query_dir/

    # Generate BRIG configuration and run
    generate_brig_config.py \
        --reference ${reference} \
        --query_dir query_dir/ \
        --output brig_output/
    """
}
