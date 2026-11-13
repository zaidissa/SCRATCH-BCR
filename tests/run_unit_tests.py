#!/usr/bin/env python3
"""
Unit tests for the analysis scripts that run without Immcantation.

    python3 tests/run_unit_tests.py            # BCR QC (needs pandas)
    python3 tests/run_unit_tests.py --with-r   # + merge_bcr_object.R (needs R: Seurat,
                                               #   scRepertoire, data.table, dplyr)

Clone calling, germlines, SHM and trees need the immcantation/suite container and are
covered by a real-data run, not here.
"""
import os, subprocess, sys, tempfile
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BIN, failures = os.path.join(ROOT, "bin"), []

def run(*cmd, cwd):
    r = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True)
    if r.returncode:
        sys.exit(f"FAILED: {' '.join(cmd)}\n{r.stdout}\n{r.stderr}")
    return r.stdout + r.stderr

def check(name, got, want):
    ok = got == want
    print(f"  {'ok  ' if ok else 'FAIL'} {name}: {got!r}" + ("" if ok else f" (want {want!r})"))
    if not ok: failures.append(name)

with tempfile.TemporaryDirectory() as d:
    run("python3", f"{ROOT}/tests/unit/make_fixture.py", d, "S1", cwd=d)
    run("python3", f"{BIN}/bcr_qc.py", "--sample", "S1", "--subject", "P1",
        "--contigs", "S1_filtered_contig_annotations.csv", "--airr", "S1_db-pass.tsv", cwd=d)
    b = pd.read_csv(f"{d}/S1_contigs_before_qc.tsv", sep="\t", keep_default_na=False)
    reason = lambda bc, chain: sorted(b[(b.barcode.str.startswith(bc + "G")) & (b.chain == chain)].qc_fail_reason)
    print("bcr_qc.py")
    check("c1 paired passes",          reason("AAAC1", "IGH") + reason("AAAC1", "IGK"), ["", ""])
    check("c2 extra light dropped",    reason("AAAC2", "IGL"), ["extra_light"])
    check("c3 multi-heavy cell",       sorted(set(reason("AAAC3", "IGH") + reason("AAAC3", "IGK"))), ["multi_heavy"])
    check("c4 light only",             reason("AAAC4", "IGK"), ["no_heavy"])
    check("c5 non-productive heavy",   reason("AAAC5", "IGH"), ["non_productive"])
    check("c6 unannotated contig",     reason("AAAC6", "IGK"), ["not_annotated"])
    check("c7 TRA ignored",            len(b[b.chain == "TRA"]), 0)
    check("c8 low-UMI heavy",          reason("AAAC8", "IGH"), ["low_umi"])
    s = pd.read_csv(f"{d}/S1_qc_summary.tsv", sep="\t").iloc[0]
    check("cells passed",              int(s.cells_passed), 5)
    check("cells paired",              int(s.cells_paired_heavy_light), 4)
    a = pd.read_csv(f"{d}/S1_airr_qc.tsv", sep="\t")
    check("airr cell_id prefixed",     a.cell_id.str.startswith("S1_").all(), True)
    check("airr cell_id no -1",        a.cell_id.str.endswith("-1").any(), False)
    check("status ok",                 s.status, "ok")

    # A sample with no contigs at all (0-byte CSV, header-only AIRR, as CHANGEO_ANNOTATE
    # writes for an empty FASTA) is recorded, not fatal.
    open(f"{d}/E1_contigs.csv", "w").close()
    with open(f"{d}/E1_db-pass.tsv", "w") as f:
        f.write("sequence_id\tlocus\tproductive\tv_call\tj_call\tjunction\n")
    run("python3", f"{BIN}/bcr_qc.py", "--sample", "E1", "--subject", "P1",
        "--contigs", "E1_contigs.csv", "--airr", "E1_db-pass.tsv", cwd=d)
    e = pd.read_csv(f"{d}/E1_qc_summary.tsv", sep="\t").iloc[0]
    check("empty sample status",       e.status, "no_contigs")
    check("empty sample cells",        int(e.cells_total), 0)

    # BTC/Cirro layout: GEX, TCR and BCR lanes under ONE sample name, no modality column,
    # a `patient` column, and a sample name Cell Ranger's --id would reject.
    print("check_samplesheet.py")
    with open(f"{d}/btc.csv", "w") as f:
        f.write("sample,fastq_1,fastq_2,patient,timepoint\n")
        for lib, n in (("BCR", 21), ("GEX", 3), ("TCR", 12)):
            for lane in (1, 2):
                stem = f"X-7-S3-{lib}-LIB/X-7-S3-{lib}-LIB_S{n}_L00{lane}"
                f.write(f"G.D7.d30+C3.v1,{stem}_R1_001.fastq.gz,{stem}_R2_001.fastq.gz,D7,30\n")
    out = run("python3", f"{BIN}/check_samplesheet.py", "btc.csv", "btc.valid.csv", cwd=d)
    v = pd.read_csv(f"{d}/btc.valid.csv", dtype=str)
    check("only BCR lanes kept",       len(v), 2)
    check("all kept lanes are BCR",    v.fastq_1.str.contains("-BCR-LIB").all(), True)
    check("patient used as subject",   set(v.subject), {"D7"})
    check("sample name kept",          set(v["sample"]), {"G.D7.d30+C3.v1"})
    check("Cell Ranger-safe id",       set(v.cr_id), {"G_D7_d30_C3_v1"})

    # A sample's rows repeated verbatim (Cirro listed its metadata twice) are dropped, not fatal.
    with open(f"{d}/dup.csv", "w") as f:
        f.write("sample,fastq_1,fastq_2,patient\n")
        for _ in range(2):
            f.write("A.d63+C3,A-BCR_S5_L001_R1_001.fastq.gz,A-BCR_S5_L001_R2_001.fastq.gz,P1\n")
        f.write("A.d63+C3,A-BCR_S5_L002_R1_001.fastq.gz,A-BCR_S5_L002_R2_001.fastq.gz,P1\n")
    run("python3", f"{BIN}/check_samplesheet.py", "dup.csv", "dup.valid.csv", cwd=d)
    check("exact duplicate rows dropped", len(pd.read_csv(f"{d}/dup.valid.csv")), 2)

    if "--with-r" in sys.argv:
        print("merge_bcr_object.R")
        a["clone_id"] = "P1_" + a.groupby("cell_id").ngroup().astype(str)
        a.loc[a.barcode.str.startswith(("AAAC1G", "AAAC10")), "clone_id"] = "P1_shared"
        a["mu_freq"] = 0.03
        a.to_csv(f"{d}/P1_germ-pass.tsv", sep="\t", index=False)
        # Two lanes of S1 in the sheet (metadata must collapse to one row per sample),
        # plus the empty sample E1, whose header-only AIRR table must not break the merge.
        with open(f"{d}/samplesheet.valid.csv", "w") as f:
            f.write("sample,fastq_1,fastq_2,modality,patient_id,timepoint,tissue,subject\n"
                    "S1,a_R1.fastq.gz,a_R2.fastq.gz,BCR,P1,D0,blood,P1\n"
                    "S1,b_R1.fastq.gz,b_R2.fastq.gz,BCR,P1,D0,blood,P1\n"
                    "E1,c_R1.fastq.gz,c_R2.fastq.gz,BCR,P1,D30,blood,P1\n")
        a.iloc[0:0].to_csv(f"{d}/E1_germ-pass.tsv", sep="\t", index=False)
        run("Rscript", f"{BIN}/merge_bcr_object.R", "--contigs-pre", "S1_contigs_before_qc.tsv",
            "--contigs-post", "S1_contigs_after_qc.tsv", "--airr", "P1_germ-pass.tsv,E1_germ-pass.tsv",
            "--samplesheet", "samplesheet.valid.csv", cwd=d)
        c = pd.read_csv(f"{d}/post_qc_cells.tsv", sep="\t")
        check("post-QC cells",             len(c), 5)
        check("shared clone size",         sorted(c.loc[c.clone_id == "P1_shared", "clone_size"]), [2, 2])
        check("isotype class",             set(c.isotype_class), {"IgM"})
        check("CTaa format",               c.CTaa.iloc[0].startswith("H:"), True)
        check("timepoint carried",         set(c.timepoint), {"D0"})
        check("tissue carried",            set(c.tissue), {"blood"})
        check("file columns not carried",  "fastq_1" in c.columns or "modality" in c.columns, False)
        check("rows not duplicated by lanes", c.cell_id.is_unique, True)
        ps = pd.read_csv(f"{d}/post_qc_summary.tsv", sep="\t").set_index("sample")
        check("summary unique_clones",     int(ps.loc["S1", "unique_clones"]), c.clone_id.nunique())
        check("summary expanded cells",    int(ps.loc["S1", "cells_in_expanded_clones"]), 2)
        check("summary median SHM",        float(ps.loc["S1", "median_mu_freq_heavy"]), 0.03)
        check("cr_id not carried",         "cr_id" in c.columns, False)
        for f in ("pre_qc_cells.tsv", "post_qc_seurat.rds", "post_qc_combineBCR.rds"):
            check(f"{f} written", os.path.exists(f"{d}/{f}"), True)

print(f"\n{'ALL PASSED' if not failures else f'{len(failures)} FAILED: {failures}'}")
sys.exit(1 if failures else 0)
