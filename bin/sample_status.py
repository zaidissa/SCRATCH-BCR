#!/usr/bin/env python3
"""
sample_status.py - why each sample did or did not make it into the results.

  cellranger  one sample, run inside CELLRANGER_VDJ after `cellranger vdj`:
              sample_status.py cellranger --sample S --subject P --exit-code N \
                  --log cellranger.log [--metrics outs/metrics_summary.csv]
              -> S_cellranger_status.tsv
              A failed run gets a reason read from Cell Ranger's own messages; a
              successful one gets its key metrics and quality flags.

  collect     the cohort, one row per samplesheet sample:
              sample_status.py collect --samplesheet samplesheet.valid.csv \
                  [--qc bcr_qc_summary.tsv] status/*.tsv
              -> sample_status.tsv

Standard library only: it runs in the Cell Ranger image as well as the Immcantation one.
"""
import argparse, csv, os, re, sys

# Cell Ranger message -> reason. First match wins; matched case-insensitively.
FAILURES = [
    ("corrupt_fastq",         r"corrupt deflate stream|invalid gzip|unexpected end of file|"
                              r"not in gzip format|crc32 mismatch|invalid block type"),
    ("no_fastqs_found",       r"no input fastqs were found|could not find fastq"),
    ("out_of_memory",         r"out of memory|cannot allocate memory|memoryerror|signal: killed"),
    ("low_valid_barcodes",    r"extremely low rate of correct barcodes|barcode.{0,40}(low|insufficient)"),
    ("reads_too_short",       r"read length.{0,80}(short|less than|insufficient)|reads? (are|is) too short"),
    ("chemistry_not_detected", r"(could not|unable to|failed to) (auto-?)?detect.{0,40}chemistry|"
                              r"chemistry.{0,60}(could not|unable|failed|incompatible|mismatch)"),
    ("too_few_reads",         r"(insufficient|too few|not enough|zero|no) (reads|read pairs|data)|"
                              r"no reads? (were )?(found|mapped)"),
    ("no_vdj_reads",          r"no (vdj|v\(d\)j)|insufficient.{0,40}v\(?d\)?j"),
    ("disk_full",             r"no space left on device|disk quota"),
]

# A run that finished but is weak. (metric, flag, test, threshold): flagged, never dropped.
FLAGS = [
    ("Estimated Number of Cells",       "few_cells",           "lt", 50),
    ("Valid Barcodes",                  "low_valid_barcodes",  "lt", 75),
    ("Reads Mapped to Any V(D)J Gene",  "low_vdj_mapping",     "lt", 50),
    ("Fraction Reads in Cells",         "low_reads_in_cells",  "lt", 50),
    ("Mean Read Pairs per Cell",        "low_read_depth",      "lt", 5000),
    ("Cells With Productive V-J Spanning Pair", "low_pairing", "lt", 30),
]
KEEP = ["Estimated Number of Cells", "Number of Read Pairs", "Mean Read Pairs per Cell",
        "Valid Barcodes", "Reads Mapped to Any V(D)J Gene", "Fraction Reads in Cells",
        "Cells With Productive V-J Spanning Pair", "Q30 Bases in RNA Read 1"]
col = lambda m: "cr_" + re.sub(r"[^a-z0-9]+", "_", m.lower()).strip("_")


def num(v):
    try:
        return float(str(v).replace(",", "").rstrip("%"))
    except ValueError:
        return None


def write(path, rows, cols):
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, cols, delimiter="\t", extrasaction="ignore", restval="")
        w.writeheader()
        w.writerows(rows)


def cellranger(a):
    row = {"sample": a.sample, "subject": a.subject, "stage": "cellranger",
           "exit_code": a.exit_code, "status": "ok", "reason": "", "flags": "", "detail": ""}
    if a.exit_code != 0:
        text = open(a.log, errors="replace").read() if os.path.exists(a.log) else ""
        # Search Cell Ranger's error block (it starts at "[error]"), not the progress lines
        # of a run that got far before failing; fall back to the log's tail.
        i = text.lower().find("[error]")
        text = text[i:] if i >= 0 else "\n".join(text.splitlines()[-60:])
        hit = next(((r, pat) for r, pat in FAILURES if re.search(pat, text, re.I)), None)
        reason = hit[0] if hit else "cellranger_error"
        # The most specific line: the one that gave the reason, else Cell Ranger's own error
        # text, else the last non-blank line.
        lines = [l.strip() for l in text.splitlines() if l.strip()]
        err = ([l for l in lines if re.search(hit[1], l, re.I)] if hit else []) or \
              [l for l in lines if re.search(r"error|fail|corrupt|invalid", l, re.I)]
        row.update(status="failed", reason=reason,
                   detail=((err or lines or ["no log"])[-1])[:300].replace("\t", " "))
    elif a.metrics and os.path.exists(a.metrics):
        with open(a.metrics, newline="") as f:
            m = next(csv.DictReader(f), {}) or {}
        row.update({col(k): m.get(k, "") for k in KEEP})
        flags = [flag for metric, flag, _, thr in FLAGS
                 if num(m.get(metric)) is not None and num(m.get(metric)) < thr]
        row["flags"] = ",".join(flags)
        if flags:
            row["status"] = "ok_low_quality"
    write(f"{a.sample}_cellranger_status.tsv", [row], list(row))


COLS = ["sample", "subject", "status", "stage", "reason", "flags", "detail", "exit_code",
        "cells_passed", "pct_cells_passed", "pct_paired"] + [col(k) for k in KEEP]


def collect(a):
    def read(p):
        with open(p, newline="") as f:
            return list(csv.DictReader(f, delimiter="\t"))
    sheet = []
    with open(a.samplesheet, newline="") as f:
        for r in csv.DictReader(f):
            if r["sample"] not in [s for s, _ in sheet]:
                sheet.append((r["sample"], r.get("subject", "")))
    cr = {r["sample"]: r for p in a.status if os.path.getsize(p) for r in read(p)}
    qc = {r["sample"]: r for r in read(a.qc)} if a.qc and os.path.getsize(a.qc) else {}

    rows = []
    for s, subject in sheet:
        row = {"sample": s, "subject": subject, "status": "ok", "stage": "", "reason": ""}
        row.update({k: v for k, v in cr.get(s, {}).items() if v != ""})
        q = qc.get(s)
        if q:
            row.update({k: q.get(k, "") for k in ("cells_passed", "pct_cells_passed", "pct_paired")})
        if cr.get(s, {}).get("status") == "failed":
            pass                                         # reason from Cell Ranger stands
        elif q and q.get("status") not in (None, "", "ok"):
            row.update(status="failed", stage="bcr_qc", reason=q["status"])
        elif not q:
            # No QC row and no recorded Cell Ranger failure: the task crashed some other
            # way (killed, retries exhausted). Its Nextflow log has the error.
            # Cell Ranger ok (or not run) means it was lost after it: annotation or QC.
            row.update(status="failed",
                       stage="annotate_or_qc" if cr.get(s, {}).get("status", "ok") != "failed" and s in cr
                             else "unknown",
                       reason=row.get("reason") or "task_failed",
                       detail=row.get("detail") or "see the pipeline log for this sample's task")
        else:
            row["stage"] = "done"
            if row.get("flags"):
                row["status"] = "ok_low_quality"
        rows.append(row)
    write("sample_status.tsv", rows, COLS)
    for r in rows:
        if r["status"] != "ok":
            print(f"{r['sample']}\t{r['status']}\t{r.get('reason') or r.get('flags')}")


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    sub = p.add_subparsers(dest="mode", required=True)
    c = sub.add_parser("cellranger")
    c.add_argument("--sample", required=True); c.add_argument("--subject", default="")
    c.add_argument("--exit-code", type=int, required=True)
    c.add_argument("--log", required=True); c.add_argument("--metrics")
    k = sub.add_parser("collect")
    k.add_argument("--samplesheet", required=True); k.add_argument("--qc")
    k.add_argument("status", nargs="*")
    a = p.parse_args()
    (cellranger if a.mode == "cellranger" else collect)(a)
