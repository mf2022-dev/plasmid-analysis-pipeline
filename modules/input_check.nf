process INPUT_CHECK {
    tag "samplesheet"

    input:
    path samplesheet

    output:
    tuple val(meta), path(sr1), path(sr2), path(lr), emit: reads
    tuple val(meta), path(assembly),                  emit: assemblies, optional: true

    script:
    """
    check_samplesheet.py ${samplesheet}
    """
}
