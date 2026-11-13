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
