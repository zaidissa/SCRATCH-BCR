/*
 * CELLRANGER_VDJ - assemble and annotate BCR contigs from 10x 5' VDJ-B FASTQs.
 *
 * Differs from SCRATCH-E2E's module only in --chain: that one is fixed to TR, this one
 * runs params.vdj_chain (IG). The 10x VDJ reference holds both TR and IG segments, so the
 * SAME reference directory serves TCR and BCR runs.
 *
 * FASTQs are staged into fastq/ and renamed to Cell Ranger's
 * <sample>_S<n>_L<lane>_R<1|2>_001.fastq.gz by SCRATCH-E2E's cellranger_renaming.py, so
 * any FASTQ naming works. Cell Ranger runs under meta.cr_id (the sample name with every
 * character outside [A-Za-z0-9_-] replaced, which its --id requires); the run folder is
 * then renamed to the samplesheet's sample name, which is the prefix of every cell_id.
 */
process CELLRANGER_VDJ {
    tag "${meta.id}"
    label 'process_high'
    container "${params.container_cellranger}"

    input:
      tuple val(meta), path(reads, stageAs: 'fastq/*')
      path reference

    output:
      tuple val(meta), path("${meta.id}/outs/filtered_contig_annotations.csv"),
                       path("${meta.id}/outs/filtered_contig.fasta"),           emit: contigs
      path "${meta.id}/outs/*",                                                emit: outs
      path "versions.yml",                                                     emit: versions

    script:
    def args = task.ext.args ?: ''
    def cr   = meta.cr_id ?: meta.id
    """
    cellranger_renaming.py "${cr}" fastq

    cellranger vdj \\
        --id="${cr}" \\
        --fastqs=fastq \\
        --sample="${cr}" \\
        --reference="${reference}" \\
        --chain=${params.vdj_chain} \\
        --localcores=${task.cpus} \\
        --localmem=${task.memory.toGiga()} \\
        ${args}

    if [ "${cr}" != "${meta.id}" ]; then mv "${cr}" "${meta.id}"; fi

    # Cell Ranger leaves ~GBs of intermediate files in <id>/SC_VDJ_ASSEMBLER_CS; only
    # outs/ is declared, so the rest stays in the work dir and is never published.
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cellranger: \$(cellranger --version 2>&1 | grep -oE '[0-9]+\\.[0-9]+\\.[0-9]+' | head -1)
    END_VERSIONS
    """

    stub:
    def cr = meta.cr_id ?: meta.id
    """
    cellranger_renaming.py "${cr}" fastq
    mkdir -p ${meta.id}/outs
    touch ${meta.id}/outs/filtered_contig_annotations.csv ${meta.id}/outs/filtered_contig.fasta \\
          ${meta.id}/outs/airr_rearrangement.tsv ${meta.id}/outs/metrics_summary.csv \\
          ${meta.id}/outs/clonotypes.csv
    printf '"%s":\\n    cellranger: stub\\n' "${task.process}" > versions.yml
    """
}
