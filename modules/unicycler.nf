process UNICYCLER {
    tag "$meta.id"
    publishDir "${params.outdir}/02_assembly/unicycler/${meta.id}", mode: 'copy'

    input:
    tuple val(meta), path(short_r1), path(short_r2), path(long_reads)

    output:
    tuple val(meta), path("${meta.id}_assembly.fasta"),  emit: assembly
    path "${meta.id}_assembly.gfa",                       emit: graph
    path "${meta.id}_unicycler.log",                      emit: log

    script:
    """
    unicycler \
        -1 ${short_r1} -2 ${short_r2} -l ${long_reads} \
        -o unicycler_out \
        --threads ${task.cpus} \
        --mode normal

    cp unicycler_out/assembly.fasta ${meta.id}_assembly.fasta
    cp unicycler_out/assembly.gfa   ${meta.id}_assembly.gfa
    cp unicycler_out/unicycler.log  ${meta.id}_unicycler.log
    """
}
