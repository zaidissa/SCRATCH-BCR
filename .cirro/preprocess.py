#!/usr/bin/env python3
"""
Cirro preprocessing for BCR-Pipeline.

1. Build the samplesheet main.nf reads, in one of two shapes:

   fastq       Cirro annotates a FASTQ dataset one row per FILE; Cell Ranger wants one row
               per sample/lane with R1 and R2 in columns. Pivoted exactly as SCRATCH-E2E's
               preprocess does (the version known to work on real Cirro FASTQ datasets).
               If the dataset's sample sheet has a `modality` column, check_samplesheet.py
               keeps only BCR rows; without one it reads the library from the FASTQ path
               (BTC-...-BCR-LIB/...), so a mixed GEX/TCR/BCR dataset - BTC datasets list
               all three under one sample name - can be selected as is.

   cellranger  An existing Cell Ranger vdj/multi output dataset: filtered_contig_annotations
               .csv + filtered_contig.fasta per sample. ds.files is Cirro's INDEXED listing,
               which often shows only some files (SCRATCH-align registers just
               airr_rearrangement.tsv) - so the two files are derived as siblings in the
               same outs/ folder of whatever file IS listed for the sample, the approach
               TCRtoolkit's singlecell preprocess uses.

2. Drop optional inputs left blank in the form, so "" never reaches Nextflow as a path.
3. Log enough to diagnose a failed launch without re-running it.

Process definition: keep `uses_sample_sheet` ENABLED for fastq input (sample names and
patient_id come from it).
"""

import pandas as pd
from cirro.helpers.preprocess_dataset import PreprocessDataset

BLANK_MEANS_UNSET = ["input_gex_object", "vdj_reference", "species"]
OPTIONAL_META = ["patient_id", "patient", "subject_id", "donor", "modality", "timepoint",
                 "tissue", "tissue_type", "sex", "batch", "condition"]
PATIENT_COLS = ("patient_id", "patient", "subject_id", "donor")


def _meta_columns(ds):
    try:
        return [c for c in OPTIONAL_META if c in ds.samplesheet.columns]
    except Exception:
        return []


def samplesheet_fastq(ds) -> pd.DataFrame:
    if "sample" not in ds.samplesheet.columns:
        raise ValueError(f"Dataset samplesheet has no 'sample' column. Present: {list(ds.samplesheet.columns)}")
    table = ds.pivot_samplesheet(
        index=["sampleIndex", "sample", "lane"],
        pivot_columns="read",
        metadata_columns=_meta_columns(ds),
        column_prefix="fastq_",
    ).sort_values(by="sample")
    missing = [c for c in ("fastq_1", "fastq_2") if c not in table.columns]
    if missing:
        raise ValueError(f"Pivoted samplesheet lacks {missing} - is this a paired FASTQ dataset? "
                         f"Columns: {list(table.columns)}")
    if not any(c in table.columns for c in PATIENT_COLS):
        ds.logger.warning("No patient_id (or patient) in the dataset sample sheet: each sample will be its own "
                          "clone-calling subject, so clones shared across a patient's samples "
                          "will NOT get one clone_id. Add patient_id to the sample metadata.")
    return table


def samplesheet_cellranger(ds) -> pd.DataFrame:
    files = ds.files.copy()
    want = {"contig_annotations": "filtered_contig_annotations.csv",
            "contig_fasta": "filtered_contig.fasta"}
    rows = []
    for sample, grp in files.groupby("sample"):
        listed = grp["file"].tolist()
        row = {"sample": sample}
        for col, name in want.items():
            exact = [f for f in listed if f.rsplit("/", 1)[-1] == name]
            if exact:
                row[col] = exact[0]
                continue
            # Derive from any file under the sample's outs/ (indexed listing is partial).
            outs = [f for f in listed if "/outs/" in f]
            if not outs:
                raise ValueError(f"{sample}: no {name} and no file under an outs/ folder to "
                                 f"derive it from. Listed: {listed[:5]}")
            row[col] = outs[0].rsplit("/outs/", 1)[0] + "/outs/" + name
            ds.logger.info(f"{sample}: {name} derived as sibling -> {row[col]}")
        rows.append(row)
    table = pd.DataFrame(rows)
    meta = _meta_columns(ds)
    if meta:
        table = table.merge(ds.samplesheet[["sample", *meta]].drop_duplicates("sample"),
                            on="sample", how="left")
    return table


def main():
    ds = PreprocessDataset.from_running()
    ds.logger.info("Files annotated in the dataset:")
    ds.logger.info(ds.files)

    input_type = str(ds.params.get("input_type", "fastq")).strip().lower()
    table = samplesheet_fastq(ds) if input_type == "fastq" else samplesheet_cellranger(ds)
    table.to_csv("samplesheet.csv", index=None)
    ds.logger.info(f"Samplesheet ({input_type}):")
    ds.logger.info(table.to_csv(index=None))
    # ${launchDir}: the head job and the tasks do not share a filesystem view on Batch.
    ds.add_param("input", "${launchDir}/samplesheet.csv", overwrite=True)

    dropped = [k for k in BLANK_MEANS_UNSET if k in ds.params and not str(ds.params[k]).strip()]
    for k in dropped:
        ds.remove_param(k)
    if dropped:
        ds.logger.info(f"Dropped blank optional params: {', '.join(dropped)}")

    ds.logger.info("Final parameters:")
    ds.logger.info(ds.params)


if __name__ == "__main__":
    main()
