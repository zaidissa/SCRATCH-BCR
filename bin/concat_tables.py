#!/usr/bin/env python3
"""
concat_tables.py OUT IN1 [IN2 ...]

Row-binds tab-separated tables with a header, taking the UNION of columns. Cell Ranger
and IgBLAST versions differ in the columns they emit, so a plain `cat`/collectFile
(which assumes identical headers) would misalign samples processed differently.
Empty/zero-byte inputs are skipped.
"""

import os
import sys

import pandas as pd

if len(sys.argv) < 3:
    sys.exit(__doc__)

out, ins = sys.argv[1], [p for p in sys.argv[2:] if os.path.getsize(p) > 0]
frames = [pd.read_csv(p, sep="\t", dtype=str, keep_default_na=False) for p in ins]
df = pd.concat(frames, ignore_index=True, sort=False) if frames else pd.DataFrame()
df.to_csv(out, sep="\t", index=False)
print(f"[concat] {out}: {len(df)} rows from {len(frames)} table(s)")
