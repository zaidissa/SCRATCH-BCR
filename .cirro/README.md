# Cirro process definition

The files in this folder register SCRATCH-BCR as a Cirro custom pipeline, the same way
SCRATCH-E2E, TCRtoolkit and SCRATCH-TCRmerge are registered.

| File | Role |
|---|---|
| `process-form.json` | the launch form (QC, clones, SHM/lineage, optional GEX object, VDJ reference override) |
| `process-input.json` | maps form values to Nextflow params; `outdir` is the new dataset's S3 folder |
| `preprocess.py` | builds `samplesheet.csv` from the selected dataset, drops blank optional params |
| `process-compute.config` | AWS Batch executor, S3 reference layout, resource ceilings |

## Register

Cirro → Pipelines → Add custom pipeline:

- **Repository**: this GitHub repo, branch `main`
- **Workflow type**: Nextflow; **entrypoint** `main.nf`
- **Configuration folder**: `.cirro`
- **Uses sample sheet**: ON (sample names and `patient`/`patient_id` come from it)
- **Input datasets**: paired FASTQ datasets (10x 5' VDJ-B, alone or mixed with GEX/TCR lanes),
  or Cell Ranger vdj/multi output datasets (form: *Start from = cellranger*)

## Inputs it handles

- **Mixed GEX/TCR/BCR FASTQ datasets** (BTC layout: one sample name, libraries told apart by
  path, e.g. `.../BTC-P7-S3-BCR-LIB/..._S21_L001_R1_001.fastq.gz`): only the BCR lanes run.
  Lanes are pivoted on (sampleIndex, sample, lane), so libraries of one sample with different
  S-numbers stay separate.
- **Patient column** `patient_id`, `patient`, `subject_id` or `donor` groups samples for clone
  calling. Without one, each sample is its own patient (a warning is logged).
- **Sample names with `.` or `+`** are kept for cell IDs; Cell Ranger runs under a cleaned copy.

## Before the first FASTQ run

The default VDJ reference is
`s3://btc-refdata/scratch/refdata-cellranger-vdj-GRCh38-alts-ensembl-7.1.0`, the directory
SCRATCH-E2E uses for TCR (the 10x VDJ reference holds IG segments too). It has not been checked
from this repo. If the run fails staging it, put the right S3 folder in the form's
**Cell Ranger VDJ reference** field; no code change is needed.

Starting from Cell Ranger output skips the reference entirely.

## Outputs (new Cirro dataset)

`bridge/merged_bcr_object/post_qc_cells.tsv` is the per-cell table SCRATCH-TCRmerge reads with
`--receptor bcr`; see the main README for the full tree.
