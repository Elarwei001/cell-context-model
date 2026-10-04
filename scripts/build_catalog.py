"""Build the context catalog from source metadata.

Usage:
    python scripts/build_catalog.py --data-dir data --sources census,h5ad,xatlas,tahoe,l1000

Per-source rows are cached under <data-dir>/sources/<name>.parquet; pass --refresh to rebuild.
Outputs go to catalog/: contexts.parquet, interventions.parquet, identities.parquet (large,
git-ignored) and summary.md, gaps.md (committed).

Source inputs:
- census: streamed from the CELLxGENE Census (no local files needed).
- h5ad:   local h5ad files of single-cell CRISPR screens in <data-dir>/h5ad/ (see H5AD_SCREENS).
- xatlas, tahoe: read from Hugging Face with column projection (no local files needed).
- l1000:  <data-dir>/l1000/siginfo_beta.txt and cellinfo_beta.txt, downloaded on first use.
"""

from __future__ import annotations

import argparse
import os
import sys
import urllib.request

import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from ccm.catalog.build import build_contexts, build_identities, build_interventions  # noqa: E402
from ccm.catalog.report import render_gaps, render_summary  # noqa: E402

H5AD_SCREENS = [
    # (source id, file name, cell line, licence)
    ("replogle2022:k562-gwps", "scperturb_replogle_k562_gwps.h5ad", "K562", "CC-BY-4.0"),
    ("replogle2022:k562-essential", "replogle_k562_ess_sc.h5ad", "K562", "CC-BY-4.0"),
    ("replogle2022:rpe1", "scperturb_replogle_rpe1.h5ad", "hTERT RPE-1", "CC-BY-4.0"),
    ("nadig2025:jurkat", "nadig_jurkat.h5ad", "Jurkat", "GEO GSE264667 (public)"),
    ("nadig2025:hepg2", "nadig_hepg2.h5ad", "HepG2", "GEO GSE264667 (public)"),
]
SOURCE_NOTES = {
    "census": "CELLxGENE Census LTS 2025-11-08 (human, mouse, macaque, marmoset, chimpanzee)",
    "h5ad": "Replogle 2022 (K562 genome-wide and essential, RPE1) and Nadig 2025 (Jurkat, HepG2)",
    "xatlas": "X-Atlas/Orion (HCT116, HEK293T)",
    "tahoe": "Tahoe-100M (50 cell lines x drugs x doses)",
    "l1000": "LINCS L1000 2020 release (bulk signatures)",
}


def _cached(path: str, refresh: bool, fn):
    if os.path.exists(path) and not refresh:
        return pd.read_parquet(path)
    df = fn()
    os.makedirs(os.path.dirname(path), exist_ok=True)
    df.to_parquet(path, index=False)
    return df


def load_source(name: str, data_dir: str, refresh: bool) -> list[pd.DataFrame]:
    cache = os.path.join(data_dir, "sources")
    if name == "census":
        from ccm.catalog.sources.census import SPECIES, census_contexts

        return [
            _cached(os.path.join(cache, f"census_{sp}.parquet"), refresh,
                    lambda sp=sp: census_contexts(sp))
            for sp in SPECIES
        ]
    if name == "h5ad":
        from ccm.catalog.sources.perturb_h5ad import H5adScreen, screen_contexts

        out = []
        for source, fname, line, licence in H5AD_SCREENS:
            path = os.path.join(data_dir, "h5ad", fname)
            if not os.path.exists(path):
                print(f"[h5ad] missing {path}, skipping {source}", flush=True)
                continue
            screen = H5adScreen(source=source, path=path, cell_line=line, licence=licence)
            out.append(_cached(os.path.join(cache, f"{source.replace(':', '_')}.parquet"),
                               refresh, lambda s=screen: screen_contexts(s)))
        return out
    if name == "xatlas":
        from ccm.catalog.sources.xatlas import xatlas_contexts

        return [_cached(os.path.join(cache, "xatlas_orion.parquet"), refresh, xatlas_contexts)]
    if name == "tahoe":
        from ccm.catalog.sources.tahoe import tahoe_contexts

        return [_cached(os.path.join(cache, "tahoe_100m.parquet"), refresh, tahoe_contexts)]
    if name == "l1000":
        from ccm.catalog.sources.l1000 import BASE, l1000_contexts

        d = os.path.join(data_dir, "l1000")
        os.makedirs(d, exist_ok=True)
        for f in ("siginfo_beta.txt", "cellinfo_beta.txt"):
            if not os.path.exists(os.path.join(d, f)):
                print(f"[l1000] downloading {f}", flush=True)
                urllib.request.urlretrieve(f"{BASE}/{f}", os.path.join(d, f))
        return [_cached(os.path.join(cache, "lincs_l1000.parquet"), refresh,
                        lambda: l1000_contexts(os.path.join(d, "siginfo_beta.txt"),
                                               os.path.join(d, "cellinfo_beta.txt")))]
    raise ValueError(f"unknown source {name}")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", default="data")
    ap.add_argument("--out-dir", default="catalog")
    ap.add_argument("--sources", default="census,h5ad,xatlas,tahoe,l1000")
    ap.add_argument("--refresh", action="store_true")
    args = ap.parse_args()
    names = [s for s in args.sources.split(",") if s]
    frames = []
    per_source = {}
    for name in names:
        loaded = load_source(name, args.data_dir, args.refresh)
        frames += loaded
        for f in loaded:
            label = str(f["source"].iloc[0]).split(":")[0] if len(f) else name
            per_source.setdefault(label, []).append(f)
        print(f"[{name}] loaded", flush=True)
    contexts = build_contexts(frames)
    per_source_contexts = {k: build_contexts(v) for k, v in per_source.items()}
    interventions = build_interventions(contexts)
    identities = build_identities(contexts)
    os.makedirs(args.out_dir, exist_ok=True)
    contexts.to_parquet(os.path.join(args.out_dir, "contexts.parquet"), index=False)
    interventions.to_parquet(os.path.join(args.out_dir, "interventions.parquet"), index=False)
    identities.to_parquet(os.path.join(args.out_dir, "identities.parquet"), index=False)
    note = "Sources: " + "; ".join(SOURCE_NOTES[n] for n in names) + "."
    with open(os.path.join(args.out_dir, "summary.md"), "w") as fh:
        fh.write(render_summary(contexts, interventions, identities, note, per_source_contexts))
    with open(os.path.join(args.out_dir, "gaps.md"), "w") as fh:
        fh.write(render_gaps(interventions, identities))
    print(f"contexts {len(contexts):,}; interventions {len(interventions):,}; "
          f"identities {len(identities):,}", flush=True)


if __name__ == "__main__":
    main()
