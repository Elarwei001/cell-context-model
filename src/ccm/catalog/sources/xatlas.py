"""X-Atlas/Orion (Xaira Therapeutics): genome-wide CRISPRi Perturb-seq in HCT116 and HEK293T.

Reads only the ``gene_target`` column of each parquet shard directly from Hugging Face (column
projection, so the expression lists are never downloaded). The shard file name encodes the cell
line (``data/<LINE>_BatchN.parquet``).
"""

from __future__ import annotations

import re
from collections import Counter

import pandas as pd

from ..schema import CELLOSAURUS, normalise_cell_line

REPO = "datasets/Xaira-Therapeutics/X-Atlas-Orion"
REVISION = "53a5bc98d49247bcf967500292575c3d3602de31"
LICENCE = "CC-BY-NC-SA-4.0"
CONTROL = "Non-Targeting"


def xatlas_contexts(revision: str = REVISION, verbose: bool = True) -> pd.DataFrame:
    import pyarrow.parquet as pq
    from huggingface_hub import HfFileSystem

    fs = HfFileSystem()
    files = sorted(fs.glob(f"{REPO}@{revision}/data/*.parquet"))
    counts: dict[str, Counter] = {}
    for i, path in enumerate(files):
        line = normalise_cell_line(re.match(r"(.+?)_Batch", path.rsplit("/", 1)[-1]).group(1))
        with fs.open(path) as fh:
            col = pq.read_table(fh, columns=["gene_target"]).column(0).to_pylist()
        counts.setdefault(line, Counter()).update(col)
        if verbose and i % 25 == 0:
            print(f"[xatlas] {i + 1}/{len(files)} shards", flush=True)
    rows = []
    for line, c in counts.items():
        for target, n in c.items():
            ctrl = target == CONTROL
            rows.append(
                {
                    "source": f"xatlas-orion:{line.lower()}",
                    "species": "human",
                    "identity_kind": "cell_line",
                    "identity": line,
                    "cell_line_id": CELLOSAURUS.get(line, ""),
                    "intervention_kind": "none" if ctrl else "genetic",
                    "intervention_target": "" if ctrl else target,
                    "intervention_detail": "" if ctrl else "CRISPRi",
                    "assay": "10x Flex",
                    "n_cells": n,
                    "licence": LICENCE,
                }
            )
    return pd.DataFrame(rows)
