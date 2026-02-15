process MOB_SUITE {
    tag "$meta.id"
    publishDir "${params.outdir}/05_plasmid_typing/mob_suite/${meta.id}", mode: 'copy'

    input:
    tuple val(meta), path(assembly)

    output:
    path "${meta.id}_mob_recon/*",   emit: results
    path "${meta.id}_plasmids/*.fasta", emit: plasmids, optional: true

    script:
    """
    mob_recon \
        --infile ${assembly} \
        --outdir ${meta.id}_mob_recon \
        --num_threads ${task.cpus}

    # Collect all plasmid FASTAs
    mkdir -p ${meta.id}_plasmids
    for f in ${meta.id}_mob_recon/plasmid_*.fasta; do
        if [ -f "\$f" ]; then
            base=\$(basename \$f .fasta)
            cp \$f ${meta.id}_plasmids/${meta.id}_\${base}.fasta
        fi
    done

    # Also run mob_typer on each plasmid
    for f in ${meta.id}_plasmids/*.fasta; do
        if [ -f "\$f" ]; then
            base=\$(basename \$f .fasta)
            mob_typer --infile \$f --outdir ${meta.id}_mob_recon/\${base}_typed/ || true
        fi
    done
    """
}
