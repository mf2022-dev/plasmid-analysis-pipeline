process FASTP {
    tag "$meta.id"
    publishDir "${params.outdir}/01_qc/fastp", mode: 'copy'

    input:
    tuple val(meta), path(reads1), path(reads2)

    output:
    tuple val(meta), path("${meta.id}_trimmed_R1.fastq.gz"), path("${meta.id}_trimmed_R2.fastq.gz"), emit: reads
    path "${meta.id}_fastp.json",  emit: json
    path "${meta.id}_fastp.html",  emit: html

    script:
    """
    fastp \
        -i ${reads1} -I ${reads2} \
        -o ${meta.id}_trimmed_R1.fastq.gz -O ${meta.id}_trimmed_R2.fastq.gz \
        --detect_adapter_pe \
        --qualified_quality_phred 20 \
        --length_required 50 \
        --thread ${task.cpus} \
        --json ${meta.id}_fastp.json \
        --html ${meta.id}_fastp.html
    """
}
