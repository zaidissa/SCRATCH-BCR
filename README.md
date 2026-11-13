# BCR-Pipeline

Single-cell **B-cell receptor** processing for 10x 5′ VDJ-B data, built to sit beside
[SCRATCH-E2E](https://github.com/zaidissa/SCRATCH-E2E) (scRNA-seq) and
[TCRtoolkit](https://github.com/KarchinLab/TCRtoolkit/tree/singlecell-master-summary) (scTCR),
and to merge with them through [SCRATCH-TCRmerge](https://github.com/zaidissa/SCRATCH-TCRmerge).

The science follows [nf-core/airrflow](https://nf-co.re/airrflow) (the Immcantation framework);
the inputs, cell naming, output bundle and Cirro packaging follow your existing pipelines, so a
GEX, a TCR and a BCR result land in **one object**.

```
10x VDJ-B FASTQ ─┐                                   (or existing Cell Ranger contigs)
                 ▼
  CELLRANGER_VDJ --chain IG ........................ per sample
  CHANGEO_ANNOTATE  IgBLAST vs IMGT (AssignGenes / MakeDb --10x)
  BCR_QC            productive, ≥1 UMI, one heavy chain per cell, top light chain
  DEFINE_CLONES     SHazaM threshold → SCOPer hierarchical clones ... per patient
  GERMLINE_SHM      dowser germlines → SHazaM mutation frequency .... per patient  [on]
  LINEAGE           dowser trees (parsimony by default) ............. per patient  [on]
  BCR_OBJECT        post_qc_cells.tsv · combineBCR.rds · Seurat .rds
  ATTACH_GEX        bcr_* columns on your GEX Seurat ................ [optional]
```

## Quick start

```bash
# 1. wiring check - no containers, no data (seconds)
nextflow run . -profile test -stub

# 2. unit tests of the QC and object code
python3 tests/run_unit_tests.py --with-r

# 3. real run (Docker)
nextflow run . \
    --input samplesheet.csv \
    --genome GRCh38 \
    --vdj_reference /refs/refdata-cellranger-vdj-GRCh38-alts-ensembl-7.1.0 \
    --outdir results
```

Mouse: `--genome GRCm39` (Cell Ranger VDJ reference GRCm38-alts-ensembl-7.0.0; IMGT/IgBLAST
`mouse`). HPC without Docker: `-profile singularity`. Apple Silicon Mac: `-profile arm64`
(the images are x86_64 and run under Rosetta).

**Resource ceilings.** Set them on the command line or in a `-params-file`:
`--max_cpus 14 --max_memory 52.GB`. A `params { max_cpus = ... }` block inside a `-c` config
file is read too late to change the ceiling (conf/base.config has already computed
`resourceLimits`); in a `-c` file set `process.resourceLimits = [cpus: 14, memory: 52.GB, time: 240.h]`
instead. A task asking for more CPUs than the machine has fails with "Process requirement
exceeds available CPUs".

## Samplesheet

**FASTQ** (`--input_type fastq`, default) - SCRATCH-E2E's layout, one row per lane:

```csv
sample,fastq_1,fastq_2,modality,patient_id,timepoint
PT1-BCR-D0,PT1-BCR-D0_S1_L001_R1_001.fastq.gz,PT1-BCR-D0_S1_L001_R2_001.fastq.gz,BCR,PT1,D0
PT1-BCR-D0,PT1-BCR-D0_S1_L002_R1_001.fastq.gz,PT1-BCR-D0_S1_L002_R2_001.fastq.gz,BCR,PT1,D0
PT1-BCR-D30,D30_R1.fastq.gz,D30_R2.fastq.gz,BCR,PT1,D30
PT1-GEX-D0,...,...,GEX,PT1,D0        <- ignored: only BCR / VDJ-B rows are used
```

- `sample` is the prefix of every `cell_id`, kept exactly as written (dots and `+` allowed) so
  it matches the GEX object's names. Cell Ranger's `--id` only takes letters, digits, `_` and
  `-`, so Cell Ranger runs under a cleaned copy (`GBM1.P7.d30+C3` -> `GBM1_P7_d30_C3`, column
  `cr_id`); output folders and cell IDs use the original name.
- `patient_id` (or any column via `--clone_by`) groups samples for **clone calling**: all of a
  patient's samples are pooled, so a clone seen at D0 and D30 gets one `clone_id`. If that column
  is missing, the first of `patient`, `subject_id`, `donor`, `donor_id` is used; with none of
  them, each sample is its own group.
- `modality` is optional; when present only `BCR`/`VDJ-B`/`VDJB`/`IG` rows run, so one
  SCRATCH samplesheet can drive both pipelines.
- **No `modality` column** (BTC/Cirro FASTQ datasets list a sample's GEX, TCR and BCR libraries
  under one sample name): the library is read from each FASTQ path. Rows whose path names
  `BCR` / `VDJ-B` / `IG` as a whole word (`BTC-P7-S3-BCR-LIB/...`) are kept; rows naming
  `GEX` / `TCR` / `VDJ-T` / `ADT` / `HTO` / ... are dropped, and the log says how many. If no
  path names any library, every row is treated as BCR.
- FASTQ names are normalised to Cell Ranger's convention, so `X_R1.fastq.gz` works. Relative
  paths resolve against the samplesheet's folder.

**Existing Cell Ranger output** (`--input_type cellranger`) - skips alignment:

```csv
sample,contig_annotations,contig_fasta,patient_id
PT1-BCR-D0,/runs/PT1-BCR-D0/outs/filtered_contig_annotations.csv,/runs/PT1-BCR-D0/outs/filtered_contig.fasta,PT1
```

For `cellranger multi` runs use `outs/per_sample_outs/<sample>/vdj_b/`.

## Outputs

```
results/
├── bridge/merged_bcr_object/      ← the per-cell BCR bundle (TCRtoolkit's layout)
│   ├── post_qc_cells.tsv            one row per QC-passing cell  ← join THIS onto GEX
│   ├── pre_qc_cells.tsv             one row per Cell Ranger cell, before QC
│   ├── {pre,post}_qc_summary.tsv    one row per sample
│   ├── {pre,post}_qc_combineBCR.rds scRepertoire::combineBCR() objects
│   └── {pre,post}_qc_seurat.rds     Seurat, BCR in @meta.data, placeholder assay
├── BCR_QC/
│   ├── tables/contigs_{before,after}_qc.tsv  one row per contig, with qc_fail_reason
│   ├── tables/bcr_qc_summary.tsv             per-sample QC counts + status (ok / no_contigs / no_cells_passed)
│   └── per_sample/
├── immcantation/
│   ├── annotate/<sample>_db-pass.tsv      IMGT-annotated AIRR rearrangements
│   ├── clones/<patient>_clone-pass.tsv    + _threshold.tsv, _dist_to_nearest.png
│   ├── germline_shm/<patient>_germ-pass.tsv   final AIRR table (germlines, mu_freq)
│   └── lineage/<patient>_trees.rds, _lineage_summary.tsv, _top_trees.pdf
├── cellranger/<sample>/outs/              Cell Ranger VDJ outputs (incl. airr_rearrangement.tsv)
├── seurat/bcr_gex_seurat.rds              only with --input_gex_object
└── pipeline_info/                         samplesheet.valid.csv, software_versions.yml, report, timeline, trace
```

### `post_qc_cells.tsv` columns

| Column | Meaning |
|---|---|
| `cell_id`, `sample`, `subject`, `barcode` | `cell_id = <sample>_<barcode>` (no `-1`), TCRtoolkit's convention |
| `patient_id`, `timepoint`, `tissue`, ... | every other samplesheet column (not file paths or `modality`), copied per sample |
| `cdr3_heavy`, `cdr3_light`, `light_locus` | CDR3 amino acids (IMGT junction without the conserved C/W·F); IGK or IGL |
| `ighv`, `ighd`, `ighj`, `iglv`, `iglj`, `ighv_allele` | IMGT gene calls (alleles stripped, except `ighv_allele`) |
| `isotype`, `isotype_class` | constant gene (e.g. IGHG1) and class (IgM/IgD/IgG/IgA/IgE) |
| `clone_id`, `clone_subgroup`, `clone_size`, `clone_size_sample`, `clone_n_samples`, `clone_size_category` | SCOPer clone (light-chain subgroup within it), cells in the clone (patient-wide / this sample), samples it spans, Single…Hyperexpanded |
| `mu_freq_heavy`, `mu_count_heavy`, `mu_freq_light`, `shm_status` | SHM; `shm_status` = unmutated if heavy `mu_freq` < 1% (`--shm_unmutated_cutoff`) |
| `in_lineage_tree` | the cell's clone got a dowser tree |
| `paired_bcr`, `multi_heavy_pre_qc`, `multi_light_pre_qc` | heavy+light both present; chain multiplicity **before** QC trimmed it |
| `CTaa`, `CTgene` | `H:<cdr3>\|L:<cdr3>` (TCRtoolkit's labelled style) and gene string |

## Failed samples

A sample whose library gives no contigs (failed library, a TCR library labelled BCR) does
not stop the run. It stays in `BCR_QC/tables/bcr_qc_summary.tsv` with `status = no_contigs`
(or `no_cells_passed` when contigs exist but every cell fails QC), the log ends with a
warning naming it, and the other samples finish. If **every** sample fails, `BCR_OBJECT`
stops with an error. A Cell Ranger crash is still a task failure: Nextflow retries it on
out-of-memory signals, then stops.

## Merging with scRNA-seq and scTCR

BCR data is per-cell **metadata**, not a second matrix: it is added to the GEX object's cells.

### Recommended: SCRATCH-TCRmerge on the `.h5ad` (any cohort size)

Use the `bcr-support` branch of SCRATCH-TCRmerge. Run it twice - TCR first, then BCR on the
result:

```bash
# 1. GEX + TCR
nextflow run SCRATCH-TCRmerge --gex_h5ad gex_annotated.h5ad \
    --tcr_tables tcrtoolkit/bridge/merged_vdj_object/post_qc_cells.tsv --outdir m_tcr
# 2. + BCR
nextflow run SCRATCH-TCRmerge --gex_h5ad m_tcr/merged.h5ad --receptor bcr \
    --tcr_tables bcr/bridge/merged_bcr_object/post_qc_cells.tsv \
    --contigs bcr/BCR_QC/tables/contigs_before_qc.tsv \
    --contigs_passed_qc bcr/BCR_QC/tables/contigs_after_qc.tsv --outdir m_tcr_bcr
```

On Cirro: the second run uses **Receptor = BCR**, with the first run's `merged.h5ad` picked as
the GEX object. The contig tables are found automatically.

Result: one `.h5ad` with `tcr_*` and `bcr_*` columns, `tcr_has_tcr_match`, `bcr_has_bcr_match`,
and **`vdj_dual_receptor`** for cells carrying both (almost always T+B doublets; flagged, never
removed). Convert to Seurat with SCRATCH-E2E's `bin/sceasy_converter.R` if the cohort fits in R.

Paired libraries named `...-GEX` / `...-TCR` / `...-BCR` match automatically: those words are
dropped from sample names when matching barcodes, and every barcode convention is scored. The
winner is reported in `tables/barcode_join_report.tsv`.

### Seurat directly (small/medium cohorts)

`--input_gex_object gex_seurat.rds` writes `seurat/bcr_gex_seurat.rds`: your object plus
`bcr_*` metadata. Barcodes are matched the same way as TCRmerge, the join is reported in
`bcr_barcode_join_report.tsv`, and it adds `vdj_dual_receptor` if TCR columns are already
present. R's sparse matrices cap at about 2.1e9 non-zeros, so use the `.h5ad` route for large
cohorts.

## SHM and lineage trees

Both are **on by default**.

| Step | Typical cost (5-20k B cells / patient) | Why |
|---|---|---|
| Clones | minutes | always on; the core output |
| SHM (`--run_shm`) | ~2-10 min per patient | naive vs antigen-experienced cells, germinal-centre activity, plasmablast maturity |
| Trees (`--run_lineage`) | `pratchet` 10-60 min · `pml` slower · `igphyml` hours | affinity maturation, clone evolution across timepoints |

Trees are built for the 500 largest clones per patient (`--lineage_max_clones`, `0` = every
clone) with at least `--lineage_min_seqs` sequences (default 3). The cap keeps an expanded
repertoire (post-vaccination, plasmablast-rich tissue) from running for many hours; clones
past it get `in_lineage_tree = FALSE`. Germlines are rebuilt whenever either step is on, since trees
are rooted on them.

## Key parameters

| Param | Default | |
|---|---|---|
| `--input_type` | `fastq` | `fastq` or `cellranger` |
| `--genome` / `--species` | `GRCh38` / derived | `GRCm39` → mouse |
| `--clone_by` | `patient_id` | samplesheet column grouping samples for clones |
| `--clonal_threshold` | `auto` | SHazaM fit per patient, else `--clonal_threshold_fallback` (0.15) |
| `--clonal_threshold_warn_above` / `_below` | `0.25` / `0.05` | a fitted threshold outside this range is flagged in `<patient>_threshold.tsv` (`threshold_note`) and the run log; it is still used |
| `--clone_method` | `hierarchical` | or `identical` |
| `--bcr_drop_multi_heavy` | `true` | drop cells with more than one productive heavy chain |
| `--bcr_multi_light` | `keep_top` | `drop_cell`, `keep_all` |
| `--lineage_build` | `pratchet` | `pml`, `igphyml` |
| `--lineage_max_clones` | `500` | trees for the N largest clones per patient; `0` = all |
| `--max_cpus` / `--max_memory` / `--max_time` | `32` / `240.GB` / `240.h` | per-task ceilings (see Quick start for `-c` files) |
| `--container_immcantation` | `immcantation/suite:4.5.0` | pin it: a new tag can change IMGT calls |

## Cirro

Register the repo with `.cirro/` as the process definition (see `.cirro/README.md`). Use a
FASTQ dataset (keep `uses_sample_sheet` on, and include `patient_id` or `patient` in the sample
metadata) or a Cell Ranger output dataset ("Start from: cellranger"). A mixed GEX/TCR/BCR
FASTQ dataset can be selected as is: only the BCR lanes run. The VDJ reference is read from
`s3://btc-refdata/scratch/refdata-cellranger-vdj-*`, the same directory SCRATCH-E2E uses for
TCR, because the 10x VDJ reference holds IG segments too. **Confirm that path before the first
FASTQ run.**

## Containers

| Image | Used by |
|---|---|
| `syedsazaidi/cellranger:v9.0.1` | CELLRANGER_VDJ (SCRATCH-E2E's image) |
| `immcantation/suite:4.5.0` | annotation, QC, clones, germlines, SHM, trees (bundles IgBLAST + IMGT refs) |
| `syedsazaidi/scratch-tcr:latest` | BCR_OBJECT, ATTACH_GEX (TCRtoolkit's single-cell image: Seurat, scRepertoire) |

## Validation

- `tests/run_unit_tests.py --with-r`: QC edge cases (one synthetic cell per failure reason),
  empty samples, the BTC samplesheet layout, metadata carry-over, per-sample summaries.
- `-profile test -stub`: wiring of every process.
- Real data (Docker, Apple Silicon via `-profile arm64`): two 10x 5' BCR libraries from tumour
  tissue, two patients, taken from a mixed GEX/TCR/BCR Cirro dataset. All 15 steps
  completed: Cell Ranger 9.0.1, IgBLAST/IMGT annotation of every contig, QC, clone calling, SHM,
  7 lineage trees, the per-cell bundle. Not yet run: on Cirro itself, the automatic threshold
  fit on a cohort large enough for it, igphyml/pml trees, the merge with a real GEX object.

## Citations

Immcantation: Gupta et al. 2015 (Change-O), Gupta et al. 2017 (SHazaM/SCOPer), Hoehn et al. 2022
(dowser); IgBLAST: Ye et al. 2013; IMGT: Lefranc et al.; scRepertoire: Borcherding et al. 2020;
nf-core/airrflow: Gabernet et al. 2024.
