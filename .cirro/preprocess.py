#!/usr/bin/env python3
"""
Cirro preprocessing for SCRATCH-BCR.

1. Build the samplesheet main.nf reads, in one of two shapes:

   fastq       Cirro annotates a FASTQ dataset one row per FILE; Cell Ranger wants one row
               per R1/R2 pair. Pairs are matched by file name (see samplesheet_fastq for why
               not pivot_samplesheet), with one metadata row per sample.
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

import re

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


# Read 1 / read 2 token in a FASTQ name: X_S1_L001_R1_001.fastq.gz, X_R1.fastq.gz, X.R1.fq.gz.
READ_RE = re.compile(r"([._])R([12])((?:_\d{3})?\.f(?:ast)?q\.gz)$")


def _sample_meta(ds) -> pd.DataFrame:
    """One metadata row per sample. A sample present in two selected datasets is listed twice
    in ds.samplesheet (seen on Cirro); merging on that would duplicate every file row."""
    meta_cols = _meta_columns(ds)
    sheet = ds.samplesheet[["sample", *meta_cols]]
    dup = sorted(sheet.loc[sheet["sample"].duplicated(), "sample"].unique())
    if dup:
        ds.logger.info(f"Samples listed more than once in the sample metadata (kept the first "
                       f"non-blank value per column): {dup}")
    return sheet.groupby("sample", as_index=False, sort=False).first()


def samplesheet_fastq(ds) -> pd.DataFrame:
    """One row per R1/R2 pair, paired by FILE NAME.

    Not ds.pivot_samplesheet(): Cirro's file table has no `lane` column and its
    `sampleIndex` is the lane number, so pivoting on (sampleIndex, sample, lane) would pair
    lane-1 reads of different libraries (GEX/TCR/BCR) listed under one sample name.
    """
    if "sample" not in ds.samplesheet.columns:
        raise ValueError(f"Dataset samplesheet has no 'sample' column. Present: {list(ds.samplesheet.columns)}")
    files = ds.files.drop_duplicates(subset=["sample", "file"])
    pairs, unpaired = {}, []
    for sample, f in zip(files["sample"], files["file"]):
        m = READ_RE.search(f)
        if not m:
            unpaired.append(f)          # index reads (I1/I2) and anything that is not R1/R2
            continue
        key = (sample, f[:m.start()] + m.group(1) + "R?" + m.group(3))
        pairs.setdefault(key, {})[f"fastq_{m.group(2)}"] = f
    if unpaired:
        ds.logger.info(f"Ignored {len(unpaired)} file(s) that are not R1/R2 reads, e.g. {unpaired[:3]}")
    rows, missing = [], []
    for (sample, _), r in sorted(pairs.items()):
        if "fastq_1" in r and "fastq_2" in r:
            rows.append({"sample": sample, **r})
        else:
            missing.append(next(iter(r.values())))
    if missing:
        raise ValueError(f"{len(missing)} FASTQ file(s) have no R1/R2 partner, e.g. {missing[:3]}. "
                         f"Is this a paired FASTQ dataset?")
    if not rows:
        raise ValueError("No R1/R2 FASTQ pairs found in the dataset.")
    table = pd.DataFrame(rows).merge(_sample_meta(ds), on="sample", how="left")
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
    if _meta_columns(ds):
        table = table.merge(_sample_meta(ds), on="sample", how="left")
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
