process FLYE {
    tag "$meta.id"
    publishDir "${params.outdir}/02_assembly/flye/${meta.id}", mode: 'copy'

    input:
    tuple val(meta), path(long_reads)

    output:
    tuple val(meta), path("${meta.id}_assembly.fasta"),  emit: assembly
    path "${meta.id}_assembly_info.txt",                  emit: info

    script:
    """
    flye \
        --nano-raw ${long_reads} \
        --out-dir flye_out \
        --threads ${task.cpus}

    cp flye_out/assembly.fasta       ${meta.id}_assembly.fasta
    cp flye_out/assembly_info.txt    ${meta.id}_assembly_info.txt
    """
}
