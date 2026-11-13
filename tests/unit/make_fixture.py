#!/usr/bin/env python3
"""
Synthetic Cell Ranger + IgBLAST/MakeDb outputs for one sample, one cell per QC edge case.

  c1  IGH + IGK                       -> passes, paired
  c2  IGH + IGK(10 umi) + IGL(2 umi)  -> IGL fails extra_light, cell passes
  c3  IGH + IGH + IGK                 -> whole cell fails multi_heavy
  c4  IGK only                        -> no_heavy
  c5  IGH non-productive + IGK        -> IGH non_productive, IGK no_heavy
  c6  IGH + IGL + IGK not in IgBLAST  -> IGK not_annotated, cell passes (IGH+IGL)
  c7  TRA contig in an IG library     -> ignored (non-IG)
  c8  IGH umi 0 + IGK                 -> IGH low_umi, IGK no_heavy
  c9  IGH only                        -> passes, unpaired
  c10 IGH + IGK, same clone as c1     -> passes (clone mate, for the object tests)

Usage: make_fixture.py OUTDIR SAMPLE
"""
import csv
import os
import sys

out, sample = sys.argv[1], sys.argv[2]
os.makedirs(out, exist_ok=True)

V = {"IGH": ("IGHV3-23*01", "IGHD3-10*01", "IGHJ4*02", "IGHM"),
     "IGK": ("IGKV1-39*01", "", "IGKJ1*01", "IGKC"),
     "IGL": ("IGLV2-14*01", "", "IGLJ2*01", "IGLC2"),
     "TRA": ("TRAV12-2*01", "", "TRAJ37*01", "TRAC")}
cells = {   # barcode: [(chain, umis, productive, annotated, cdr3)]
    "AAAC1": [("IGH", 20, True, True, "ARDYW"), ("IGK", 8, True, True, "QQSYST")],
    "AAAC2": [("IGH", 15, True, True, "ARGGW"), ("IGK", 10, True, True, "QQYNS"), ("IGL", 2, True, True, "QSYDS")],
    "AAAC3": [("IGH", 12, True, True, "ARPPW"), ("IGH", 9, True, True, "ARTTW"), ("IGK", 5, True, True, "QQANS")],
    "AAAC4": [("IGK", 7, True, True, "QQLNS")],
    "AAAC5": [("IGH", 11, False, True, "ARKKW"), ("IGK", 6, True, True, "QQHHS")],
    "AAAC6": [("IGH", 14, True, True, "ARSSW"), ("IGL", 9, True, True, "QSWDS"), ("IGK", 3, True, False, "QQXXS")],
    "AAAC7": [("TRA", 9, True, False, "AVNRA")],
    "AAAC8": [("IGH", 0, True, True, "ARMMW"), ("IGK", 4, True, True, "QQMMS")],
    "AAAC9": [("IGH", 18, True, True, "ARQQW")],
    "AAAC10": [("IGH", 16, True, True, "ARDYW"), ("IGK", 7, True, True, "QQSYST")],
}

cr_rows, airr_rows = [], []
for bc, contigs in cells.items():
    barcode = f"{bc}GGGGGGGGGGG"[:16] + "-1"
    for i, (chain, umis, prod, ann, cdr3) in enumerate(contigs, start=1):
        cid = f"{barcode}_contig_{i}"
        v, d, j, c = V[chain]
        cr_rows.append(dict(sample=f"{sample}_cr", barcode=barcode, is_cell="true", contig_id=cid, high_confidence="true",
                            length=550, chain=chain, v_gene=v.split("*")[0], d_gene=d.split("*")[0],
                            j_gene=j.split("*")[0], c_gene=c, full_length="true",
                            productive=str(prod).lower(), fwr1="", fwr1_nt="", cdr1="", cdr1_nt="",
                            fwr2="", fwr2_nt="", cdr2="", cdr2_nt="", fwr3="", fwr3_nt="",
                            cdr3=cdr3, cdr3_nt="TGT" * len(cdr3), fwr4="", fwr4_nt="",
                            reads=umis * 10, umis=umis, raw_clonotype_id=f"clonotype{bc}",
                            raw_consensus_id=f"clonotype{bc}_consensus_{i}"))
        if ann and chain.startswith("IG"):
            airr_rows.append(dict(sequence_id=cid, sequence="ACGT" * 100, rev_comp="F",
                                  productive="T" if prod else "F", v_call=v, d_call=d, j_call=j,
                                  c_call=c, locus=chain, junction=f"TGT{'GCC' * len(cdr3)}TGG",
                                  junction_aa=f"C{cdr3}W", junction_length=3 * len(cdr3) + 6,
                                  sequence_alignment="ACGT" * 90, germline_alignment="ACGT" * 90,
                                  v_germline_start=1, cell_id=barcode, umi_count=umis,
                                  consensus_count=umis * 10))

with open(f"{out}/{sample}_filtered_contig_annotations.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(cr_rows[0])); w.writeheader(); w.writerows(cr_rows)
with open(f"{out}/{sample}_db-pass.tsv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(airr_rows[0]), delimiter="\t"); w.writeheader(); w.writerows(airr_rows)
print(f"fixture: {len(cr_rows)} contigs, {len(airr_rows)} AIRR rows, {len(cells)} cells -> {out}")
