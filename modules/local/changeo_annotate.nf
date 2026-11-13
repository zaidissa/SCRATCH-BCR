/*
 * CHANGEO_ANNOTATE - re-annotate Cell Ranger's BCR contigs against IMGT with IgBLAST.
 *
 * Cell Ranger's own V(D)J calls use the 10x reference, which lacks IMGT gapping. Every
 * Immcantation step after this (SHazaM distances, SCOPer clones, germline reconstruction,
 * SHM, Dowser trees) needs IMGT-gapped alignments (sequence_alignment, v_germline_start
 * ...), so contigs are re-annotated exactly as in the Immcantation 10x tutorial and
 * nf-core/airrflow:
 *
 *   AssignGenes.py igblast   -> <id>_igblast.fmt7
 *   MakeDb.py igblast --10x  -> <id>_db-pass.tsv (AIRR; Cell Ranger's c_call, umi_count,
 *                               cell_id merged in from filtered_contig_annotations.csv)
 *
 * The IgBLAST database and IMGT references ship inside the immcantation/suite image
 * (params.igblast_db, params.imgt_germlines). Pin the image tag for reproducible calls:
 * a new tag can carry a new IMGT release.
 */
process CHANGEO_ANNOTATE {
    tag "${meta.id}"
    label 'process_medium'
    container "${params.container_immcantation}"

    input:
      tuple val(meta), path(annotations), path(fasta)
      val species

    output:
      tuple val(meta), path(annotations), path("${meta.id}_db-pass.tsv"), emit: db
      path "${meta.id}_db-fail.tsv",                                      emit: failed, optional: true
      path "versions.yml",                                                emit: versions

    script:
    """
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        changeo: \$(MakeDb.py --version 2>&1 | grep -oE '[0-9]+\\.[0-9]+\\.[0-9]+' | head -1)
        igblast: \$(igblastn -version 2>&1 | grep -oE '[0-9]+\\.[0-9]+\\.[0-9]+' | head -1)
    END_VERSIONS

    # A sample with no contigs (failed library, TCR library run as BCR) must not stop the
    # cohort: pass a header-only AIRR table on, so BCR_QC records it as status=no_contigs
    # and the run warns about it at the end.
    if [ ! -s "${fasta}" ]; then
        echo "WARNING: ${fasta} is empty for ${meta.id} - Cell Ranger called no contigs." >&2
        echo "Check that the library is VDJ-B and that cellranger ran with --chain IG." >&2
        printf 'sequence_id\\tsequence\\trev_comp\\tproductive\\tv_call\\td_call\\tj_call\\tc_call\\tsequence_alignment\\tgermline_alignment\\tjunction\\tjunction_aa\\tlocus\\tcell_id\\tumi_count\\n' \\
            > ${meta.id}_db-pass.tsv
        exit 0
    fi

    AssignGenes.py igblast \\
        -s ${fasta} \\
        -b ${params.igblast_db} \\
        --organism ${species} \\
        --loci ig \\
        --format blast \\
        --nproc ${task.cpus} \\
        --outname ${meta.id}

    MakeDb.py igblast \\
        -i ${meta.id}_igblast.fmt7 \\
        -s ${fasta} \\
        -r ${params.imgt_germlines}/${species}/vdj \\
        --10x ${annotations} \\
        --extended \\
        --failed \\
        --outname ${meta.id}
    """

    stub:
    """
    touch ${meta.id}_db-pass.tsv
    printf '"%s":\\n    changeo: stub\\n' "${task.process}" > versions.yml
    """
}
