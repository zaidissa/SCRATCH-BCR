process SAMPLESHEET_CHECK {
    tag "${samplesheet.name}"
    label 'process_single'
    container "${params.container_immcantation}"

    input:
      path samplesheet

    output:
      path "samplesheet.valid.csv", emit: csv

    script:
    """
    check_samplesheet.py ${samplesheet} samplesheet.valid.csv \\
        --input-type ${params.input_type} \\
        --subject-col '${params.clone_by}'
    """

    // Pure stdlib Python, so the stub runs the real check: a stub test then exercises the
    // same validation a real run does.
    stub:
    """
    check_samplesheet.py ${samplesheet} samplesheet.valid.csv \\
        --input-type ${params.input_type} \\
        --subject-col '${params.clone_by}'
    """
}
