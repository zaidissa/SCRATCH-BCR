#!/usr/bin/env python3
"""
bcr_qc.py - per-sample BCR contig QC.

Joins Cell Ranger's filtered_contig_annotations.csv to the IgBLAST/IMGT re-annotation
(Change-O MakeDb.py --10x output, one AIRR row per contig) and applies single-cell BCR
QC at the CELL level, because clone calling and every downstream readout is per cell.

Rules, applied in order (each contig keeps the FIRST reason it failed):

  not_annotated      IgBLAST produced no AIRR row for the contig (MakeDb fail)
  non_productive     IMGT productive != T                  (--productive-only)
  low_umi            umi_count < --min-umis
  multi_heavy        the cell has >1 passing IGH contig     (--drop-multi-heavy)
                     Two productive heavy chains in one barcode is a doublet signal far
                     more often than allelic inclusion; Immcantation drops these cells.
  no_heavy           the cell has no passing IGH contig     (--require-heavy)
                     Clones are defined on the heavy chain, so a light-only cell cannot
                     be placed in a clone.
  extra_light        the cell has >1 passing light contig and --multi-light keep_top:
                     the highest-UMI light is kept, the rest fail with this reason.
                     drop_cell removes the whole cell instead (reason multi_light).

Outputs, named so SCRATCH-TCRmerge reads them exactly like TCRtoolkit's VDJ_QC tables:

  <sample>_contigs_before_qc.tsv  every Cell Ranger contig + sample + IMGT calls +
                                  passed_bcr_qc / qc_fail_reason. Cell Ranger column
                                  names are kept, because scRepertoire::combineBCR() and
                                  TCRmerge's per-contig table read them.
  <sample>_contigs_after_qc.tsv   the passing subset of the above.
  <sample>_airr_qc.tsv            passing AIRR rows (Immcantation input for clone calling),
                                  with cell_id = <sample>_<barcode> so cells stay unique
                                  when samples of one subject are pooled.
  <sample>_qc_summary.tsv         one row of counts.
"""

import argparse
import os
import sys

import pandas as pd

LIGHT = ("IGK", "IGL")


def read_table(path, sep, empty_cols=()):
    """A 0-byte file (no contigs at all) reads as an empty table with `empty_cols`."""
    if os.path.getsize(path) == 0:
        return pd.DataFrame(columns=list(empty_cols), dtype=str)
    return pd.read_csv(path, sep=sep, dtype=str, keep_default_na=False, na_values=[""])


def truthy(series):
    return series.astype(str).str.upper().isin(["T", "TRUE", "1", "YES"])


def bare(barcode):
    """AAACCTG...-1 -> AAACCTG... (TCRtoolkit's per-cell key drops the GEM suffix)."""
    return barcode.astype(str).str.replace(r"-\d+$", "", regex=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--subject", default="", help="patient/subject id (clone-calling group)")
    ap.add_argument("--contigs", required=True, help="Cell Ranger filtered_contig_annotations.csv")
    ap.add_argument("--airr", required=True, help="MakeDb.py igblast --10x db-pass.tsv")
    ap.add_argument("--min-umis", type=int, default=1)
    ap.add_argument("--productive-only", default="true")
    ap.add_argument("--drop-multi-heavy", default="true")
    ap.add_argument("--require-heavy", default="true")
    ap.add_argument("--multi-light", default="keep_top", choices=("keep_top", "drop_cell", "keep_all"))
    ap.add_argument("--outdir", default=".")
    a = ap.parse_args()

    flag = lambda v: str(v).lower() in ("true", "1", "yes", "y", "on")
    s, subject = a.sample, (a.subject or a.sample)

    cr = read_table(a.contigs, ",", ("barcode", "contig_id", "chain", "umis"))
    for col in ("barcode", "contig_id", "chain"):
        if col not in cr.columns:
            sys.exit(f"ERROR: {a.contigs} has no '{col}' column - is it filtered_contig_annotations.csv?")
    # Cell Ranger also reports TR contigs in an IG run (and vice versa) when the library
    # carries some; this pipeline is about B cells, so only IG chains are considered.
    n_non_ig = int((~cr["chain"].isin(("IGH",) + LIGHT)).sum())
    cr = cr[cr["chain"].isin(("IGH",) + LIGHT)].copy()

    airr = read_table(a.airr, "\t", ("sequence_id", "locus", "productive", "v_call", "j_call", "junction"))
    need = {"sequence_id", "locus", "productive", "v_call", "j_call", "junction"}
    missing = need - set(airr.columns)
    if missing:
        sys.exit(f"ERROR: {a.airr} lacks AIRR columns {sorted(missing)}")
    umi_col = next((c for c in ("umi_count", "duplicate_count", "consensus_count") if c in airr.columns), None)

    # ── contig-level table: Cell Ranger rows + IMGT calls ───────────────────
    imgt_cols = {c: f"imgt_{c}" for c in
                 ("locus", "productive", "v_call", "d_call", "j_call", "c_call", "junction_aa")
                 if c in airr.columns}
    ann = airr[["sequence_id", *imgt_cols]].rename(columns=imgt_cols)
    ct = cr.merge(ann, left_on="contig_id", right_on="sequence_id", how="left").drop(columns="sequence_id")
    # Cell Ranger 9 writes its own `sample` column (the Cell Ranger --id, made safe for it);
    # the samplesheet name is the one every cell_id uses, so it replaces that column.
    ct = ct.drop(columns=[c for c in ("sample", "subject") if c in ct.columns])
    ct.insert(0, "sample", s)
    ct.insert(1, "subject", subject)
    ct["cell_id"] = s + "_" + bare(ct["barcode"])
    ct["umis_num"] = pd.to_numeric(ct.get("umis", pd.Series(index=ct.index, dtype=float)), errors="coerce").fillna(0)

    reason = pd.Series("", index=ct.index, dtype=object)

    def fail(mask, why):
        reason.loc[mask & (reason == "")] = why

    fail(ct["imgt_locus"].isna(), "not_annotated")
    if flag(a.productive_only):
        fail(~truthy(ct["imgt_productive"].fillna("")), "non_productive")
    fail(ct["umis_num"] < a.min_umis, "low_umi")

    # The locus that counts is IMGT's when present (it re-derives it), Cell Ranger's otherwise.
    locus = ct["imgt_locus"].fillna(ct["chain"])
    alive = reason == ""
    n_heavy = (alive & (locus == "IGH")).groupby(ct["cell_id"]).transform("sum")
    n_light = (alive & locus.isin(LIGHT)).groupby(ct["cell_id"]).transform("sum")

    multi_heavy_cell = n_heavy > 1
    if flag(a.drop_multi_heavy):
        fail(multi_heavy_cell, "multi_heavy")
    if flag(a.require_heavy):
        fail(n_heavy == 0, "no_heavy")

    if a.multi_light == "drop_cell":
        fail(n_light > 1, "multi_light")
    elif a.multi_light == "keep_top":
        alive = reason == ""
        lights = ct[alive & locus.isin(LIGHT)].sort_values("umis_num", ascending=False)
        extra = lights.index[lights.duplicated(subset="cell_id", keep="first")]
        fail(ct.index.isin(extra), "extra_light")

    ct["passed_bcr_qc"] = reason == ""
    ct["qc_fail_reason"] = reason
    ct["multi_heavy_pre_qc"] = multi_heavy_cell.values
    ct["multi_light_pre_qc"] = (n_light > 1).values
    ct = ct.drop(columns="umis_num")

    passed = ct[ct["passed_bcr_qc"]]

    # ── AIRR rows for clone calling ──────────────────────────────────────────
    aq = airr[airr["sequence_id"].isin(set(passed["contig_id"]))].copy()
    bc_of = dict(zip(passed["contig_id"], passed["barcode"]))
    aq["barcode"] = aq["sequence_id"].map(bc_of)
    aq["cell_id"] = s + "_" + bare(aq["barcode"])          # was the bare 10x barcode
    aq["sample"], aq["subject"] = s, subject
    # Cell Ranger names contigs <barcode>-1_contig_N, so two libraries of one patient can
    # share a sequence_id; clones, germlines and tree tips are per patient and need it unique.
    aq["sequence_id"] = s + "_" + aq["sequence_id"]
    if umi_col and umi_col != "umi_count":
        aq["umi_count"] = aq[umi_col]

    # ── summary ──────────────────────────────────────────────────────────────
    cells_all = ct["cell_id"].nunique()
    cells_pass = passed["cell_id"].nunique()
    per_cell_pass = passed.assign(L=passed["imgt_locus"]).groupby("cell_id")["L"].agg(set)
    paired = int(per_cell_pass.map(lambda x: "IGH" in x and bool(x & set(LIGHT))).sum())
    # status lets the run (and a reader of bcr_qc_summary.tsv) tell a failed library from
    # a small one: the sample stays in every table instead of stopping the cohort.
    status = "no_contigs" if len(ct) == 0 else ("no_cells_passed" if cells_pass == 0 else "ok")
    summary = pd.DataFrame([{
        "sample": s, "subject": subject, "status": status,
        "contigs_ig": len(ct), "contigs_non_ig_ignored": n_non_ig,
        "contigs_annotated": int(ct["imgt_locus"].notna().sum()),
        "contigs_passed": len(passed),
        "cells_total": cells_all, "cells_passed": cells_pass,
        "pct_cells_passed": round(100 * cells_pass / max(cells_all, 1), 2),
        "cells_paired_heavy_light": paired,
        "pct_paired": round(100 * paired / max(cells_pass, 1), 2),
        "cells_multi_heavy": int(ct.loc[ct["multi_heavy_pre_qc"], "cell_id"].nunique()),
        "cells_multi_light": int(ct.loc[ct["multi_light_pre_qc"], "cell_id"].nunique()),
        **{f"contigs_failed_{r}": int((reason == r).sum())
           for r in ("not_annotated", "non_productive", "low_umi", "multi_heavy",
                     "no_heavy", "multi_light", "extra_light")},
    }])

    o = a.outdir.rstrip("/")
    ct.to_csv(f"{o}/{s}_contigs_before_qc.tsv", sep="\t", index=False)
    passed.to_csv(f"{o}/{s}_contigs_after_qc.tsv", sep="\t", index=False)
    aq.to_csv(f"{o}/{s}_airr_qc.tsv", sep="\t", index=False)
    summary.to_csv(f"{o}/{s}_qc_summary.tsv", sep="\t", index=False)

    r = summary.iloc[0]
    print(f"[bcr_qc] {s}: {r.cells_passed}/{r.cells_total} cells pass ({r.pct_cells_passed}%), "
          f"{r.cells_paired_heavy_light} paired H+L, {r.cells_multi_heavy} multi-heavy, "
          f"{len(passed)}/{len(ct)} IG contigs kept")
    if status != "ok":
        print(f"[bcr_qc] WARNING {s}: {status} - this sample contributes no cells. Check the "
              f"chain type of this library (was Cell Ranger run with --chain IG?)", file=sys.stderr)


if __name__ == "__main__":
    main()
