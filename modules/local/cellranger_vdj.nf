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
 *
 * Every run writes <sample>_cellranger_status.tsv (bin/sample_status.py): the failure
 * reason read from Cell Ranger's message (corrupt_fastq, chemistry_not_detected,
 * too_few_reads, ...) or, on success, key metrics and low-quality flags. A data failure
 * with --skip_failed_samples true exits 0 with no contigs, so the sample drops out with
 * its reason recorded; an out-of-memory failure exits 137 so it is retried with more.
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
                       path("${meta.id}/outs/filtered_contig.fasta"),           emit: contigs, optional: true
      path "${meta.id}/outs/*",                                                emit: outs,    optional: true
      path "${meta.id}_cellranger_status.tsv",                                 emit: status
      path "cellranger.log",                                                   emit: log
      path "versions.yml",                                                     emit: versions

    script:
    def args = task.ext.args ?: ''
    def cr   = meta.cr_id ?: meta.id
    def skip = params.skip_failed_samples.toString() == 'true'
    def last = task.attempt > 2
    """
    cellranger_renaming.py "${cr}" fastq

    set +e
    cellranger vdj \\
        --id="${cr}" \\
        --fastqs=fastq \\
        --sample="${cr}" \\
        --reference="${reference}" \\
        --chain=${params.vdj_chain} \\
        --localcores=${task.cpus} \\
        --localmem=${task.memory.toGiga()} \\
        ${args} 2>&1 | tee cellranger.log
    rc=\${PIPESTATUS[0]}
    set -e

    if [ \$rc -eq 0 ] && [ "${cr}" != "${meta.id}" ]; then mv "${cr}" "${meta.id}"; fi
    sample_status.py cellranger --sample '${meta.id}' --subject '${meta.subject}' \\
        --exit-code \$rc --log cellranger.log --metrics "${meta.id}/outs/metrics_summary.csv" \\
        --pipestance "${cr}"

    if [ \$rc -ne 0 ]; then
        reason=\$(awk -F'\\t' 'NR==1{for(i=1;i<=NF;i++) if(\$i=="reason") c=i} NR==2{print \$c}' ${meta.id}_cellranger_status.tsv)
        echo "Cell Ranger failed for ${meta.id}: \$reason" >&2
        # More memory may fix it: fail the task so Nextflow retries with a bigger one.
        if [ "\$reason" = out_of_memory ] && [ "${last}" = false ]; then exit 137; fi
        # Data problem: keep the status row, drop the sample (or stop the run).
        rm -rf "${cr}"
        if [ "${skip}" = false ]; then exit \$rc; fi
    fi

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
    touch cellranger.log
    printf 'sample\\tsubject\\tstage\\tstatus\\treason\\n${meta.id}\\t${meta.subject}\\tcellranger\\tok\\t\\n' > ${meta.id}_cellranger_status.tsv
    mkdir -p ${meta.id}/outs
    touch ${meta.id}/outs/filtered_contig_annotations.csv ${meta.id}/outs/filtered_contig.fasta \\
          ${meta.id}/outs/airr_rearrangement.tsv ${meta.id}/outs/metrics_summary.csv \\
          ${meta.id}/outs/clonotypes.csv
    printf '"%s":\\n    cellranger: stub\\n' "${task.process}" > versions.yml
    """
}
