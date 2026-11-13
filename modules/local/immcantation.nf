/*
 * Per-subject Immcantation steps. A "subject" is the samplesheet's --clone_by column
 * (default patient_id): every sample of one subject is pooled, so clones shared across a
 * patient's timepoints/tissues get one clone_id.
 *
 *   DEFINE_CLONES  SHazaM threshold + SCOPer hierarchical clones        (always)
 *   GERMLINE_SHM   dowser germlines + SHazaM observed mutations         (run_shm || run_lineage)
 *   LINEAGE        dowser lineage trees                                 (run_lineage)
 */

process DEFINE_CLONES {
    tag "${subject}"
    label 'process_medium'
    container "${params.container_immcantation}"

    input:
      tuple val(subject), path(airr, stageAs: 'in/*')

    output:
      tuple val(subject), path("${subject}_clone-pass.tsv"), emit: clones
      path "${subject}_threshold.tsv",                      emit: threshold
      path "${subject}_dist_to_nearest.png",                emit: plot, optional: true
      path "versions.yml",                                  emit: versions

    script:
    """
    define_clones.R \\
        --subject '${subject}' \\
        --threshold ${params.clonal_threshold} \\
        --fallback-threshold ${params.clonal_threshold_fallback} \\
        --method ${params.clone_method} \\
        --warn-above ${params.clonal_threshold_warn_above} \\
        --warn-below ${params.clonal_threshold_warn_below} \\
        --nproc ${task.cpus} \\
        ${(airr instanceof List ? airr : [airr]).join(' ')}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        shazam: \$(Rscript -e 'cat(as.character(packageVersion("shazam")))' 2>/dev/null)
        scoper: \$(Rscript -e 'cat(as.character(packageVersion("scoper")))' 2>/dev/null)
        dowser: \$(Rscript -e 'cat(as.character(packageVersion("dowser")))' 2>/dev/null)
        airr: \$(Rscript -e 'cat(as.character(packageVersion("airr")))' 2>/dev/null)
    END_VERSIONS
    """

    stub:
    """
    printf 'sequence_id\\tcell_id\\tclone_id\\n' > ${subject}_clone-pass.tsv
    printf 'subject\\tthreshold\\n${subject}\\t0.15\\n' > ${subject}_threshold.tsv
    printf '"%s":\\n    scoper: stub\\n' "${task.process}" > versions.yml
    """
}

process GERMLINE_SHM {
    tag "${subject}"
    label 'process_medium'
    container "${params.container_immcantation}"

    input:
      tuple val(subject), path(clones)
      val species
      val run_shm

    output:
      tuple val(subject), path("${subject}_germ-pass.tsv"), emit: db
      path "versions.yml",                                  emit: versions

    script:
    """
    germline_shm.R \\
        --subject '${subject}' \\
        --input ${clones} \\
        --germline-dir ${params.imgt_germlines} \\
        --species ${species} \\
        --run-shm ${run_shm} \\
        --nproc ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dowser: \$(Rscript -e 'cat(as.character(packageVersion("dowser")))' 2>/dev/null)
        shazam: \$(Rscript -e 'cat(as.character(packageVersion("shazam")))' 2>/dev/null)
    END_VERSIONS
    """

    stub:
    """
    cp ${clones} ${subject}_germ-pass.tsv
    printf '"%s":\\n    dowser: stub\\n' "${task.process}" > versions.yml
    """
}

process LINEAGE {
    tag "${subject}"
    label 'process_lineage'
    container "${params.container_immcantation}"

    input:
      tuple val(subject), path(db)

    output:
      path "${subject}_lineage_summary.tsv", emit: summary
      path "${subject}_trees.rds",           emit: trees, optional: true
      path "${subject}_top_trees.pdf",       emit: plots, optional: true
      path "versions.yml",                   emit: versions

    script:
    """
    build_lineage.R \\
        --subject '${subject}' \\
        --input ${db} \\
        --build ${params.lineage_build} \\
        --min-seqs ${params.lineage_min_seqs} \\
        --max-clones ${params.lineage_max_clones} \\
        --chain ${params.lineage_chain} \\
        --plot-top ${params.lineage_plot_top} \\
        --igphyml ${params.igphyml_exec} \\
        --nproc ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dowser: \$(Rscript -e 'cat(as.character(packageVersion("dowser")))' 2>/dev/null)
        phangorn: \$(Rscript -e 'cat(as.character(packageVersion("phangorn")))' 2>/dev/null)
    END_VERSIONS
    """

    stub:
    """
    printf 'subject\\tclone_id\\n' > ${subject}_lineage_summary.tsv
    printf '"%s":\\n    dowser: stub\\n' "${task.process}" > versions.yml
    """
}
