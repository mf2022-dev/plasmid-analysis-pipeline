process HYBRACTER {
    tag "$meta.id"
    publishDir "${params.outdir}/02_assembly/hybracter/${meta.id}", mode: 'copy'

    input:
    tuple val(meta), path(short_r1), path(short_r2), path(long_reads)

    output:
    tuple val(meta), path("${meta.id}_assembly.fasta"),  emit: assembly

    script:
    """
    hybracter hybrid \
        -l ${long_reads} \
        -1 ${short_r1} -2 ${short_r2} \
        --output hybracter_out \
        --threads ${task.cpus}

    cp hybracter_out/final_assembly.fasta ${meta.id}_assembly.fasta
    """
}
