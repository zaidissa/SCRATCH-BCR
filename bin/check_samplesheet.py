#!/usr/bin/env python3
"""
check_samplesheet.py SHEET OUT --input-type fastq|cellranger [--subject-col patient_id]

Validates the BCR samplesheet and writes a normalised copy (samplesheet.valid.csv).

fastq       sample,fastq_1,fastq_2[,modality][,patient_id][,...]
            SCRATCH-E2E's layout: one row per lane, several rows per sample allowed. With a
            `modality` column only BCR rows are kept (BCR / VDJ-B / VDJB / IG); GEX and TCR
            rows are dropped with a note, so one SCRATCH samplesheet can drive both.
cellranger  sample,contig_annotations,contig_fasta[,patient_id][,...]
            existing Cell Ranger vdj/multi outputs (filtered_contig_annotations.csv and
            filtered_contig.fasta), skipping alignment.

Library type without a `modality` column: Cirro FASTQ datasets often list a sample's GEX,
TCR and BCR libraries under ONE sample name, telling them apart only by the path
(BTC-...-BCR-LIB/..., ..._GEX_...). Each row's FASTQ path is then read for a library word:
rows marked BCR / VDJ-B / IG are kept, rows marked GEX / TCR / VDJ-T / ADT / HTO / ATAC /
CRISPR are dropped. If no row carries a library word, every row is assumed to be BCR.

A `subject` column is added from --subject-col (clones are called per subject). If that
column is missing, the first of patient_id, patient, subject_id, donor, donor_id that the
sheet has is used; when none is present, or a value is blank, the sample is its own subject.

`cr_id` is the sample name made safe for Cell Ranger, whose --id accepts only letters,
digits, '_' and '-' ("GBM1.DFCI7.d30+C3" -> "GBM1_DFCI7_d30_C3"). It is used only for the
Cell Ranger run; `sample` - and so every cell_id - keeps the original name, which is the
name the GEX object (SCRATCH-E2E) uses. Every other column is carried
through as sample-level metadata. File existence is NOT checked here - paths may be S3
URIs; Nextflow checks them when it stages the files.
"""

import argparse
import csv
import re
import sys

BCR_MODALITIES = {"BCR", "VDJ-B", "VDJB", "VDJ_B", "IG"}
REQUIRED = {"fastq": ["sample", "fastq_1", "fastq_2"],
            "cellranger": ["sample", "contig_annotations", "contig_fasta"]}
SAMPLE_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._+-]*$")
SUBJECT_FALLBACKS = ("patient_id", "patient", "subject_id", "donor", "donor_id")
# Library words in a FASTQ path, matched as whole tokens of the upper-cased path.
_TOK = r"(?<![A-Z0-9])({})(?![A-Z0-9])"
BCR_PATH_RE = re.compile(_TOK.format(r"BCR|VDJ[-_]?B|IG|BCELL"))
OTHER_PATH_RE = re.compile(_TOK.format(r"GEX|TCR|VDJ[-_]?T|ADT|HTO|CITE|ATAC|CRISPR|CMO|FB|5GEX|3GEX"))


def library_from_path(*paths):
    """'BCR', 'other', or '' from the library words in the FASTQ paths."""
    s = " ".join(paths).upper()
    if BCR_PATH_RE.search(s):
        return "BCR"
    return "other" if OTHER_PATH_RE.search(s) else ""


def cellranger_id(sample):
    return re.sub(r"_+", "_", re.sub(r"[^A-Za-z0-9_-]", "_", sample)).strip("_")


def die(msg):
    sys.exit(f"ERROR (samplesheet): {msg}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("sheet")
    ap.add_argument("out")
    ap.add_argument("--input-type", default="fastq", choices=sorted(REQUIRED))
    ap.add_argument("--subject-col", default="patient_id")
    a = ap.parse_args()

    with open(a.sheet, newline="") as fh:
        rows = list(csv.DictReader(fh))
    if not rows:
        die(f"{a.sheet} has no rows")
    cols = [c.strip() for c in rows[0].keys()]
    rows = [{k.strip(): (v or "").strip() for k, v in r.items() if k is not None} for r in rows]

    # Exact duplicate rows (a launcher listing one sample's metadata twice) are dropped;
    # the same file in two DIFFERENT rows is still an error below.
    uniq = list({tuple(sorted(r.items())): r for r in rows}.values())
    if len(uniq) < len(rows):
        print(f"[samplesheet] dropped {len(rows) - len(uniq)} exact duplicate row(s)")
        rows = uniq

    missing = [c for c in REQUIRED[a.input_type] if c not in cols]
    if missing:
        die(f"missing column(s) {missing} for --input_type {a.input_type}. Found: {cols}")

    if "modality" not in cols and a.input_type == "fastq":
        libs = [library_from_path(r["fastq_1"], r["fastq_2"]) for r in rows]
        if any(libs):
            kept = [dict(r, modality="BCR") for r, lib in zip(rows, libs) if lib == "BCR"]
            n_other = sum(lib == "other" for lib in libs)
            n_unknown = sum(lib == "" for lib in libs)
            print(f"[samplesheet] no modality column: library read from FASTQ paths - "
                  f"{len(kept)} BCR row(s) kept, {n_other} GEX/TCR/other row(s) dropped"
                  + (f", {n_unknown} row(s) with no library word dropped" if n_unknown else ""))
            if not kept:
                die("no FASTQ path names a BCR library (BCR / VDJ-B / IG) and some name another "
                    "library (GEX / TCR / ...). Add a modality column marking the BCR rows.")
            rows = kept
            cols = cols + ["modality"]

    if "modality" in cols:
        kept = [r for r in rows if r["modality"].upper() in BCR_MODALITIES]
        dropped = sorted({r["modality"] or "(blank)" for r in rows} - {r["modality"] for r in kept})
        if dropped:
            print(f"[samplesheet] ignoring rows with modality {dropped} (not BCR)")
        rows = kept
        if not rows:
            die(f"no row has a BCR modality ({sorted(BCR_MODALITIES)}). "
                f"Is this the right samplesheet, or should the modality column be BCR?")

    subject_col = a.subject_col
    if subject_col not in cols:
        subject_col = next((c for c in SUBJECT_FALLBACKS if c in cols), a.subject_col)
        if subject_col != a.subject_col:
            print(f"[samplesheet] no '{a.subject_col}' column: grouping clones by '{subject_col}'")

    files = REQUIRED[a.input_type][1:]
    seen_files = set()
    for i, r in enumerate(rows, start=2):
        s = r["sample"]
        if not SAMPLE_RE.match(s):
            die(f"row {i}: sample '{s}' - use letters, digits, '.', '_', '+' or '-', no spaces. "
                f"It is the prefix of every cell_id.")
        r["cr_id"] = cellranger_id(s)
        for f in files:
            if not r[f]:
                die(f"row {i} ({s}): empty {f}")
            if r[f] in seen_files:
                die(f"row {i} ({s}): {r[f]} is listed twice")
            seen_files.add(r[f])
        if a.input_type == "fastq":
            for f in files:
                if not re.search(r"\.f(ast)?q\.gz$", r[f]):
                    die(f"row {i} ({s}): {f} '{r[f]}' is not a .fastq.gz")
        r["subject"] = r.get(subject_col, "") or s
        if not SAMPLE_RE.match(r["subject"]):
            die(f"row {i} ({s}): {subject_col} '{r['subject']}' - use letters, digits, "
                f"'.', '_', '+' or '-'. It names the per-subject clone files and prefixes clone_id.")

    if a.input_type == "cellranger":
        dup = sorted({r["sample"] for r in rows if [x["sample"] for x in rows].count(r["sample"]) > 1})
        if dup:
            die(f"samples listed more than once with --input_type cellranger: {dup}")

    # Two sample names that only differ in characters Cell Ranger rejects would collide.
    by_id = {}
    for r in rows:
        if by_id.setdefault(r["cr_id"], r["sample"]) != r["sample"]:
            die(f"samples '{by_id[r['cr_id']]}' and '{r['sample']}' both become Cell Ranger id "
                f"'{r['cr_id']}'. Rename one of them.")

    # One subject per sample - a sample split across two subjects would split its cells
    # into two clone-calling groups.
    subj = {}
    for r in rows:
        if subj.setdefault(r["sample"], r["subject"]) != r["subject"]:
            die(f"sample {r['sample']} has more than one {subject_col}: "
                f"{subj[r['sample']]} and {r['subject']}")

    out_cols = REQUIRED[a.input_type] + ["subject", "cr_id"] + \
        [c for c in cols if c not in REQUIRED[a.input_type] + ["subject", "cr_id"]]
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=out_cols, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)

    n_samples = len({r["sample"] for r in rows})
    n_subjects = len(set(subj.values()))
    print(f"[samplesheet] {len(rows)} row(s), {n_samples} sample(s), {n_subjects} subject(s) "
          f"(clones grouped by {subject_col if subject_col in cols else 'sample'})")


if __name__ == "__main__":
    main()
