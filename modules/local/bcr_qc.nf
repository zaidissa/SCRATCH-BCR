/*
 * BCR_QC - per-sample cell-level QC (see bin/bcr_qc.py for the rules), and
 * BCR_QC_COLLECT - cohort tables named exactly like TCRtoolkit's VDJ_QC outputs
 * (contigs_before_qc.tsv / contigs_after_qc.tsv), published to BCR_QC/tables/, which is
 * where SCRATCH-TCRmerge (receptor = BCR) derives them from.
 */
process BCR_QC {
    tag "${meta.id}"
    label 'process_low'
    container "${params.container_immcantation}"

    input:
      tuple val(meta), path(annotations), path(db)

    output:
      tuple val(meta), path("${meta.id}_airr_qc.tsv"), emit: airr
      path "${meta.id}_contigs_before_qc.tsv",         emit: before
      path "${meta.id}_contigs_after_qc.tsv",          emit: after
      path "${meta.id}_qc_summary.tsv",                emit: summary
      path "versions.yml",                             emit: versions

    script:
    """
    bcr_qc.py \\
        --sample '${meta.id}' \\
        --subject '${meta.subject}' \\
        --contigs ${annotations} \\
        --airr ${db} \\
        --min-umis ${params.bcr_min_umis} \\
        --productive-only ${params.bcr_productive_only} \\
        --drop-multi-heavy ${params.bcr_drop_multi_heavy} \\
        --require-heavy ${params.bcr_require_heavy} \\
        --multi-light ${params.bcr_multi_light}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 -c 'import platform; print(platform.python_version())')
        pandas: \$(python3 -c 'import pandas; print(pandas.__version__)')
    END_VERSIONS
    """

    stub:
    """
    for t in airr_qc contigs_before_qc contigs_after_qc qc_summary; do
        printf 'sample\\tsubject\\tstatus\\n${meta.id}\\t${meta.subject}\\tok\\n' > ${meta.id}_\$t.tsv
    done
    printf '"%s":\\n    pandas: stub\\n' "${task.process}" > versions.yml
    """
}

process BCR_QC_COLLECT {
    tag "${summaries.size()} sample(s)"
    label 'process_low'
    container "${params.container_immcantation}"

    input:
      path before,    stageAs: 'before/*'
      path after,     stageAs: 'after/*'
      path summaries, stageAs: 'summary/*'

    output:
      path "contigs_before_qc.tsv", emit: before
      path "contigs_after_qc.tsv",  emit: after
      path "bcr_qc_summary.tsv",    emit: summary

    script:
    """
    concat_tables.py contigs_before_qc.tsv before/*
    concat_tables.py contigs_after_qc.tsv  after/*
    concat_tables.py bcr_qc_summary.tsv    summary/*
    """

    stub:
    """
    concat_tables.py contigs_before_qc.tsv before/* || touch contigs_before_qc.tsv
    concat_tables.py contigs_after_qc.tsv  after/*  || touch contigs_after_qc.tsv
    concat_tables.py bcr_qc_summary.tsv    summary/* || touch bcr_qc_summary.tsv
    """
}
