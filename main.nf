#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SCRATCH-BCR - single-cell B-cell receptor processing
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    10x 5' VDJ-B FASTQ (or existing Cell Ranger contigs)
      -> CELLRANGER_VDJ --chain IG                         per sample
      -> CHANGEO_ANNOTATE  IgBLAST vs IMGT                 per sample
      -> BCR_QC            productive, 1 heavy / cell ...  per sample, then cohort tables
      -> DEFINE_CLONES     SHazaM threshold + SCOPer       per subject (patient)
      -> GERMLINE_SHM      dowser germlines + SHazaM       per subject   [run_shm]
      -> LINEAGE           dowser trees                    per subject   [run_lineage]
      -> BCR_OBJECT        per-cell bundle (cells.tsv, combineBCR, Seurat)
      -> ATTACH_GEX        BCR columns on a GEX Seurat     [input_gex_object]

    The science follows nf-core/airrflow (Immcantation). The output contract follows
    TCRtoolkit + SCRATCH-TCRmerge, so GEX, TCR and BCR end up in one object:
    see README "Merging with scRNA-seq and scTCR".
----------------------------------------------------------------------------------------
*/

include { SAMPLESHEET_CHECK            } from './modules/local/samplesheet_check.nf'
include { CELLRANGER_VDJ               } from './modules/local/cellranger_vdj.nf'
include { CHANGEO_ANNOTATE             } from './modules/local/changeo_annotate.nf'
include { BCR_QC; BCR_QC_COLLECT       } from './modules/local/bcr_qc.nf'
include { DEFINE_CLONES; GERMLINE_SHM; LINEAGE } from './modules/local/immcantation.nf'
include { BCR_OBJECT; ATTACH_GEX       } from './modules/local/bcr_object.nf'

// Functions, not `def x = { }` closures: Nextflow 26.x's strict parser rejects calling a
// closure variable like a function (the same trap SCRATCH-E2E documents).

// Cirro sends JSON booleans; the command line sends the STRING "true"/"false".
def truthy(value) {
    value?.toString()?.toLowerCase() in ['true', '1', 'yes', 'y', 'on']
}

// Absolute paths and URIs pass through; a relative path resolves against the
// SAMPLESHEET's directory, so a sheet works wherever the run is launched from.
def resolveSheetPath(p, sheetDir) {
    def t = p.toString()
    return (t.startsWith('/') || t =~ /^[a-zA-Z0-9]+:\/\//) ? t : sheetDir.resolve(t).toString()
}

workflow {

    // ── parameters ───────────────────────────────────────────────────────────
    if (!params.input) {
        error "No samplesheet. Set --input samplesheet.csv (see README: Samplesheet)."
    }
    if (!(params.input_type in ['fastq', 'cellranger'])) {
        error "--input_type must be 'fastq' or 'cellranger', got '${params.input_type}'."
    }
    def species = params.species ?: (params.genome.toString().startsWith('GRCm') ? 'mouse' : 'human')
    if (!(species in ['human', 'mouse'])) {
        error "--species must be human or mouse (IMGT/IgBLAST references in the container), got '${species}'."
    }
    def run_shm     = truthy(params.run_shm)
    def run_lineage = truthy(params.run_lineage)
    if (!(params.lineage_build in ['pratchet', 'pml', 'igphyml'])) {
        error "--lineage_build must be pratchet, pml or igphyml, got '${params.lineage_build}'."
    }
    log.info """
    SCRATCH-BCR  ${workflow.manifest.version}
      input       : ${params.input}  (${params.input_type})
      genome      : ${params.genome}  -> species ${species}
      clones by   : ${params.clone_by} (or patient / subject_id / donor if absent; see samplesheet log), threshold ${params.clonal_threshold} (${params.clone_method})
      SHM         : ${run_shm ? 'on' : 'off'}
      lineage     : ${run_lineage ? "on (${params.lineage_build}, clones >= ${params.lineage_min_seqs} seqs)" : 'off'}
      GEX object  : ${params.input_gex_object ?: '-'}
      outdir      : ${params.outdir}
    """.stripIndent()

    // ── samplesheet ──────────────────────────────────────────────────────────
    def sheetDir = file(params.input).parent
    ch_sheet = SAMPLESHEET_CHECK(channel.fromPath(params.input, checkIfExists: true)).csv
    ch_rows  = ch_sheet.splitCsv(header: true)
    ch_versions = channel.empty()

    // ── contigs: Cell Ranger, or supplied ────────────────────────────────────
    if (params.input_type == 'fastq') {
        def ref = params.vdj_reference ?: params.genomes?.get(params.genome)?.vdj
        if (!ref) {
            error "No Cell Ranger VDJ reference for --genome ${params.genome}. Set --vdj_reference " +
                  "or --genomes.${params.genome}.vdj (10x refdata-cellranger-vdj-*)."
        }
        // Local paths are checked now, before a 16-core task is queued for nothing;
        // remote URIs are left to the staging step (credentials live there).
        if (!(ref.toString() ==~ /^[a-z0-9+.-]+:\/\/.*/) && !file(ref).exists()) {
            error "Cell Ranger VDJ reference not found: ${ref}\n" +
                  "  Download refdata-cellranger-vdj-* from 10x Genomics and set --vdj_reference,\n" +
                  "  or start from existing contigs with --input_type cellranger."
        }
        // One row per lane -> one task per sample with every lane's R1/R2.
        ch_reads = ch_rows
            .map { r ->
                tuple(r.sample, r.subject, r.cr_id,
                      [ file(resolveSheetPath(r.fastq_1, sheetDir), checkIfExists: true),
                        file(resolveSheetPath(r.fastq_2, sheetDir), checkIfExists: true) ])
            }
            .groupTuple(by: [0, 1, 2])
            .map { sample, subject, cr_id, reads ->
                tuple([id: sample, subject: subject, cr_id: cr_id], reads.flatten()) }

        CELLRANGER_VDJ(ch_reads, file(ref))
        ch_contigs  = CELLRANGER_VDJ.out.contigs
        ch_versions = ch_versions.mix(CELLRANGER_VDJ.out.versions)
    } else {
        ch_contigs = ch_rows.map { r ->
            tuple([id: r.sample, subject: r.subject],
                  file(resolveSheetPath(r.contig_annotations, sheetDir), checkIfExists: true),
                  file(resolveSheetPath(r.contig_fasta, sheetDir),       checkIfExists: true))
        }
    }

    // ── annotation + QC ──────────────────────────────────────────────────────
    CHANGEO_ANNOTATE(ch_contigs, species)
    BCR_QC(CHANGEO_ANNOTATE.out.db)
    BCR_QC_COLLECT(BCR_QC.out.before.collect(), BCR_QC.out.after.collect(), BCR_QC.out.summary.collect())
    ch_versions = ch_versions.mix(CHANGEO_ANNOTATE.out.versions, BCR_QC.out.versions)

    // A sample with no usable cells stays in the tables (status column) instead of
    // stopping the cohort; say so where it will be seen.
    BCR_QC.out.summary
        .flatMap { f -> f.splitCsv(sep: '\t', header: true) }
        .filter { r -> r.status && r.status != 'ok' }
        .subscribe { r ->
            log.warn "Sample ${r.sample}: ${r.status} - it contributes no BCR cells. " +
                     "See BCR_QC/tables/bcr_qc_summary.tsv."
        }

    // ── clones, per subject ──────────────────────────────────────────────────
    ch_subject = BCR_QC.out.airr
        .map { meta, f -> tuple(meta.subject, f) }
        .groupTuple()
    DEFINE_CLONES(ch_subject)
    ch_final = DEFINE_CLONES.out.clones
    ch_versions = ch_versions.mix(DEFINE_CLONES.out.versions)

    DEFINE_CLONES.out.threshold
        .flatMap { f -> f.splitCsv(sep: '\t', header: true) }
        .filter { r -> r.threshold_note && r.threshold_note != 'NA' }
        .subscribe { r ->
            log.warn "Patient ${r.subject}: clonal threshold ${r.threshold_note}. Check " +
                     "immcantation/clones/${r.subject}_dist_to_nearest.png; set --clonal_threshold to override."
        }

    // Germlines are needed by SHM and by lineage trees (trees root on the germline).
    if (run_shm || run_lineage) {
        GERMLINE_SHM(ch_final, species, run_shm)
        ch_final = GERMLINE_SHM.out.db
        ch_versions = ch_versions.mix(GERMLINE_SHM.out.versions)
    }

    def no_lineage = file("${projectDir}/assets/NO_FILE_lineage", checkIfExists: true)
    ch_lineage = channel.value([no_lineage])
    if (run_lineage) {
        LINEAGE(GERMLINE_SHM.out.db)
        ch_lineage = LINEAGE.out.summary.collect().ifEmpty([no_lineage])
        ch_versions = ch_versions.mix(LINEAGE.out.versions)
    }

    // ── per-cell object ──────────────────────────────────────────────────────
    BCR_OBJECT(
        BCR_QC_COLLECT.out.before,
        BCR_QC_COLLECT.out.after,
        ch_final.map { _subject, f -> f }.collect(),
        ch_lineage,
        ch_sheet
    )
    ch_versions = ch_versions.mix(BCR_OBJECT.out.versions)

    if (params.input_gex_object) {
        ATTACH_GEX(file(params.input_gex_object, checkIfExists: true), BCR_OBJECT.out.post_cells)
        ch_versions = ch_versions.mix(ATTACH_GEX.out.versions)
    }

    // ── software versions: one file for the methods section ─────────────────
    // Every task writes the same block for its process, so de-duplicate on content.
    ch_versions
        .map { f -> f.text.trim() }
        .unique()
        .mix(channel.of("Workflow:\n    SCRATCH-BCR: ${workflow.manifest.version}\n    Nextflow: ${nextflow.version}".toString()))
        .collectFile(name: 'software_versions.yml', storeDir: "${params.outdir}/pipeline_info",
                     newLine: true, sort: true)

    workflow.onComplete = {
        log.info(workflow.success
            ? "SCRATCH-BCR finished. Per-cell BCR table: ${params.outdir}/bridge/merged_bcr_object/post_qc_cells.tsv"
            : "SCRATCH-BCR failed - see the error above.")
    }
}
