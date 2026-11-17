/*
 * BCR_OBJECT - the per-cell BCR bundle, published to bridge/merged_bcr_object/.
 *
 * Same file names as TCRtoolkit's bridge/merged_vdj_object/, in a folder of its own, so
 * SCRATCH-TCRmerge finds it with receptor = BCR and never confuses it with a TCR bundle
 * when both datasets are selected. See bin/merge_bcr_object.R.
 *
 * ATTACH_GEX - optional: add the post-QC BCR columns to an existing GEX Seurat object
 * (--input_gex_object). For large cohorts use SCRATCH-TCRmerge on the .h5ad instead.
 */
process BCR_OBJECT {
    tag "${params.project_name}"
    label 'process_medium'
    container "${params.container_sc}"

    input:
      path contigs_before
      path contigs_after
      path airr,    stageAs: 'airr/*'
      path lineage, stageAs: 'lineage/*'      // NO_FILE_lineage when lineage trees are off
      path samplesheet                        // samplesheet.valid.csv: per-sample metadata

    output:
      path "post_qc_cells.tsv", emit: post_cells
      path "pre_qc_cells.tsv",  emit: pre_cells
      path "*_summary.tsv",     emit: summary
      path "*.rds",             emit: rds, optional: true
      path "versions.yml",      emit: versions

    script:
    def lin = (lineage instanceof List ? lineage : [lineage]).findAll { f -> !f.name.startsWith('NO_FILE') }
    """
    merge_bcr_object.R \\
        --contigs-pre ${contigs_before} \\
        --contigs-post ${contigs_after} \\
        --airr '${(airr instanceof List ? airr : [airr]).join(',')}' \\
        --lineage '${lin.join(',')}' \\
        --unmutated-cutoff ${params.shm_unmutated_cutoff} \\
        --samplesheet ${samplesheet} \\
        --project '${params.project_name}'

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Seurat: \$(Rscript -e 'cat(as.character(packageVersion("Seurat")))' 2>/dev/null)
        scRepertoire: \$(Rscript -e 'cat(as.character(packageVersion("scRepertoire")))' 2>/dev/null)
        data.table: \$(Rscript -e 'cat(as.character(packageVersion("data.table")))' 2>/dev/null)
    END_VERSIONS
    """

    stub:
    """
    for s in pre_qc post_qc; do
        printf 'cell_id\\tsample\\tbarcode\\n' > \${s}_cells.tsv
        printf 'sample\\tcells\\n' > \${s}_summary.tsv
    done
    printf '"%s":\\n    Seurat: stub\\n' "${task.process}" > versions.yml
    """
}

process ATTACH_GEX {
    tag "${gex.name}"
    label 'process_high_memory'
    container "${params.container_sc}"

    input:
      path gex
      path cells

    output:
      path "bcr_gex_seurat.rds",           emit: seurat
      path "bcr_barcode_join_report.tsv",  emit: report
      path "versions.yml",                 emit: versions

    script:
    """
    attach_bcr_to_seurat.R \\
        --seurat ${gex} \\
        --cells ${cells} \\
        --prefix '${params.obs_prefix}' \\
        --library-tokens '${params.library_tokens}' \\
        --min-match-rate ${params.min_match_rate} \\
        --out bcr_gex_seurat.rds

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Seurat: \$(Rscript -e 'cat(as.character(packageVersion("Seurat")))' 2>/dev/null)
    END_VERSIONS
    """

    stub:
    """
    touch bcr_gex_seurat.rds bcr_barcode_join_report.tsv
    printf '"%s":\\n    Seurat: stub\\n' "${task.process}" > versions.yml
    """
}

// Comma-joined staged files, without the NO_FILE placeholder.
def realFiles(x) {
    (x instanceof List ? x : [x]).findAll { f -> !f.name.startsWith('NO_FILE') }.join(',')
}

/*
 * BCR_REPORT - one self-contained HTML report (bin/bcr_report.R): sample status, cell QC,
 * isotypes, clonal expansion, clones shared across samples/time points, SHM, IGHV usage,
 * the clonal threshold and lineage trees. Opens offline; light and dark page.
 */
process BCR_REPORT {
    tag "${params.project_name}"
    label 'process_low'
    container "${params.container_immcantation}"

    input:
      path cells
      path qc,         stageAs: 'qc/*'
      path status,     stageAs: 'status/*'
      path thresholds, stageAs: 'threshold/*'
      path dist_plots, stageAs: 'dist/*'
      path lineage,    stageAs: 'lineage/*'
      path trees,      stageAs: 'trees/*'

    output:
      path "bcr_report.html", emit: html

    script:
    """
    bcr_report.R \\
        --cells ${cells} \\
        --qc '${realFiles(qc)}' \\
        --status '${realFiles(status)}' \\
        --thresholds '${realFiles(thresholds)}' \\
        --dist-plots '${realFiles(dist_plots)}' \\
        --lineage '${realFiles(lineage)}' \\
        --trees '${realFiles(trees)}' \\
        --project '${params.project_name}' \\
        --version '${workflow.manifest.version}' \\
        --out bcr_report.html
    """

    stub:
    """
    echo '<!doctype html><title>SCRATCH-BCR report</title>' > bcr_report.html
    """
}
