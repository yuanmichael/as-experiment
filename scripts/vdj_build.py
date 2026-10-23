"""Download the paired GEX + V(D)J manifest and build per-library training data.

One row of training data is ONE CELL: its receptor, and its own transcriptome.
The barcode is what links them, so this script's whole job is to perform that
join correctly and record what it had to throw away.

Measured on 7 samples beforehand: ~82% of T-cell-sorted GEX cells carry a
productive TRB. A library that comes back far below that is a mismatch, not
biology -- GSE167118's Bac17L paired against the CD3-NEGATIVE barcodes and
reported 0.8%, which reads exactly like poor recovery. Every library's pairing
rate is written out so that failure cannot hide.

Two antibody-capture traps are handled explicitly:
  * 10x .h5 files carry a feature_type axis; Antibody Capture rows are dropped.
    GSE277117 taught this the hard way -- 163 antibody features sat in the gene
    matrix, confounded with the outcome.
  * features.tsv third column is the same thing for the mtx format.

Outputs, per library, to /data/vdj/<gse>/<library>/:
    cells.parquet   barcode, cdr3_aa, v_gene, j_gene, clone_size, gse, library
    expr.npz        CSR counts (cells x genes), rows aligned to cells.parquet
    genes.txt       gene symbols for that matrix's columns
    stats.json      pairing rate and what was dropped

Gene harmonisation across studies is deliberately a SEPARATE later step: each
library keeps its own gene axis here, so a bad intersection cannot silently
destroy data at ingest.

    modal run scripts/vdj_build.py::plan                    # cost, no download
    modal run --detach scripts/vdj_build.py::build_all
    modal run scripts/vdj_build.py::report
"""

import json
import pathlib

import modal

app = modal.App("as-experiment")
vol = modal.Volume.from_name("as-data", create_if_missing=True)
OUT = "/data/vdj"

image = (
    modal.Image.debian_slim(python_version="3.12")
    .apt_install("wget")
    .pip_install("pandas", "numpy", "scipy", "h5py", "pyarrow")
)

MANIFEST = pathlib.Path(__file__).resolve().parents[1] / "data" / "airr" / "vdj_manifest.tsv"

# A library whose pairing rate falls below this is flagged, not silently used.
MIN_PAIRING = 0.25


def load_manifest() -> list[dict]:
    import csv
    with open(MANIFEST, newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    for r in rows:
        r["bytes"] = int(r.get("bytes") or 0)
    return rows


@app.function(image=image, volumes={"/data": vol}, timeout=600)
def plan():
    """What this will cost and how much it will move, before moving any of it."""
    rows = load_manifest.local() if hasattr(load_manifest, "local") else None
    print("plan() runs locally; see the caller output")


@app.function(image=image, volumes={"/data": vol}, cpu=2.0, memory=8192,
              timeout=3600, max_containers=20, retries=2)
def build_one(row: dict) -> dict:
    """Download one library's files, join barcodes, write training arrays."""
    import gzip
    import os
    import shutil
    import subprocess
    import time

    import h5py
    import numpy as np
    import pandas as pd
    import scipy.io
    import scipy.sparse as sp

    gse, lib = row["gse"], row["library"]
    safe = "".join(c if c.isalnum() or c in "-_." else "_" for c in lib) or "series"
    dest = f"{OUT}/{gse}/{safe}"
    if os.path.exists(f"{dest}/stats.json") and not row.get("_force"):
        with open(f"{dest}/stats.json") as fh:
            st = json.load(fh)
        # Recompute if the record predates the T-cell enrichment diagnostic --
        # without it a low pairing rate cannot be told from a bad join.
        if "tcell_enrichment" in st:
            return {**st, "status": "cached"}

    work = f"/tmp/{gse}_{safe}"
    os.makedirs(work, exist_ok=True)
    os.makedirs(dest, exist_ok=True)
    t0 = time.time()
    FTP = "https://ftp.ncbi.nlm.nih.gov/geo"

    def url_for(name: str) -> str | None:
        import re
        m = re.match(r"(GSM\d+)", name)
        if m:
            g = m.group(1)
            return f"{FTP}/samples/{g[:-3]}nnn/{g}/suppl/{name}"
        if name.startswith(gse):
            return f"{FTP}/series/{gse[:-3]}nnn/{gse}/suppl/{name}"
        return None

    def grab(name: str) -> str | None:
        if not name:
            return None
        u = url_for(name)
        if not u:
            return None
        p = f"{work}/{os.path.basename(name)}"
        r = subprocess.run(["wget", "-q", "-t", "3", "-O", p, u])
        if r.returncode != 0 or not os.path.exists(p) or os.path.getsize(p) == 0:
            return None
        return p

    # ---------------------------------------------------------------- contigs
    cpath = grab(row["contig"])
    if not cpath:
        return {"gse": gse, "library": lib, "status": "contig download failed"}
    op = gzip.open if cpath.endswith(".gz") else open
    with op(cpath, "rt", errors="replace") as fh:
        ctg = pd.read_csv(fh, low_memory=False)
    ctg.columns = [c.strip().lower() for c in ctg.columns]

    def col(*names):
        for n in names:
            if n in ctg.columns:
                return n
        return None

    bc_c = col("barcode", "cell_id", "cellbarcode")
    cdr_c = col("cdr3", "cdr3_aa", "cdr3aa")
    if not bc_c or not cdr_c:
        return {"gse": gse, "library": lib, "status": "contig schema",
                "columns": list(ctg.columns)[:12]}
    ch_c, pr_c = col("chain"), col("productive")
    v_c, j_c = col("v_gene"), col("j_gene")

    trb = ctg
    n_contigs = len(trb)
    if ch_c:
        trb = trb[trb[ch_c].astype(str).str.upper() == "TRB"]
    if pr_c:
        trb = trb[trb[pr_c].astype(str).str.lower().isin(
            ("true", "t", "yes", "productive"))]
    trb = trb[trb[cdr_c].notna()]
    trb = trb[~trb[cdr_c].astype(str).str.contains(r"[\*_]", na=False, regex=True)]
    # one receptor per cell: keep the most-supported contig
    umi_c = col("umis", "umi", "reads")
    if umi_c:
        trb = trb.sort_values(umi_c, ascending=False)
    trb = trb.drop_duplicates(subset=[bc_c], keep="first")

    # ------------------------------------------------------------- expression
    genes = None
    if row.get("h5"):
        hp = grab(row["h5"])
        if not hp:
            return {"gse": gse, "library": lib, "status": "h5 download failed"}
        with h5py.File(hp, "r") as f:
            grp = "matrix" if "matrix" in f else list(f.keys())[0]
            g = f[grp]
            data = g["data"][:]
            indices = g["indices"][:]
            indptr = g["indptr"][:]
            shape = tuple(g["shape"][:])
            bcs = [b.decode() for b in g["barcodes"][:]]
            feat = g["features"] if "features" in g else g
            names = feat["name"][:] if "name" in feat else feat["gene_names"][:]
            genes = [n.decode() for n in names]
            ftype = ([t.decode() for t in feat["feature_type"][:]]
                     if "feature_type" in feat else None)
        # 10x h5 is genes x cells in CSC; transpose to cells x genes
        M = sp.csc_matrix((data, indices, indptr), shape=shape).T.tocsr()
    else:
        bp, fp, mp = grab(row["barcodes"]), grab(row["features"]), grab(row["matrix"])
        if not (bp and mp):
            return {"gse": gse, "library": lib, "status": "mtx download failed"}
        op = gzip.open if bp.endswith(".gz") else open
        with op(bp, "rt", errors="replace") as fh:
            bcs = [ln.strip() for ln in fh if ln.strip()]
        ftype = None
        if fp:
            op = gzip.open if fp.endswith(".gz") else open
            with op(fp, "rt", errors="replace") as fh:
                parts = [ln.rstrip("\n").split("\t") for ln in fh if ln.strip()]
            genes = [p[1] if len(p) > 1 else p[0] for p in parts]
            if parts and len(parts[0]) > 2:
                ftype = [p[2] for p in parts]
        if mp.endswith(".gz"):
            with gzip.open(mp, "rb") as fin, open(mp[:-3], "wb") as fo:
                shutil.copyfileobj(fin, fo)
            mp = mp[:-3]
        M = sp.csr_matrix(scipy.io.mmread(mp).T)   # mtx is genes x cells

    if genes is None:
        genes = [f"g{i}" for i in range(M.shape[1])]

    # Drop Antibody Capture / CITE features. GSE277117: 163 antibody features
    # hid in the gene matrix and were confounded with the outcome.
    n_ab = 0
    if ftype is not None and len(ftype) == M.shape[1]:
        keep = np.array([t.strip().lower() in ("gene expression", "gene_expression")
                         for t in ftype])
        if keep.any() and not keep.all():
            n_ab = int((~keep).sum())
            M = M[:, keep]
            genes = [g for g, k in zip(genes, keep) if k]

    # ------------------------------------------------------------------- join
    gex_idx = {b: i for i, b in enumerate(bcs)}
    hit = trb[bc_c].astype(str).map(gex_idx)
    paired = trb[hit.notna()].copy()
    paired["_row"] = hit[hit.notna()].astype(int).values

    rate = len(paired) / max(len(bcs), 1)
    if len(paired) == 0:
        st = {"gse": gse, "library": lib, "status": "zero overlap",
              "n_gex": len(bcs), "n_trb": len(trb), "pairing_rate": 0.0}
        with open(f"{dest}/stats.json", "w") as fh:
            json.dump(st, fh)
        vol.commit()
        return st

    # Is a low pairing rate a mismatch, or just a sample that isn't T-sorted?
    #
    # Rate alone cannot tell these apart. An unsorted biopsy where 15% of cells
    # are T cells SHOULD pair at ~15%, and that is fine -- we only train on the
    # cells that paired. A wrong join (GSE167118's Bac17L against CD3-negative
    # barcodes) also gives a low rate, but the cells it pairs are arbitrary.
    #
    # The discriminating test: paired cells must actually be T cells. Compare
    # CD3/TRAC expression in paired versus unpaired cells. A correct join shows
    # strong enrichment at ANY rate; a mismatch shows none.
    tmark = [i for i, g in enumerate(genes)
             if g.upper() in ("CD3E", "CD3D", "CD3G", "TRAC", "TRBC2")]
    enrich = None
    if tmark and len(paired) and len(bcs) > len(paired):
        rows_p = paired["_row"].values
        mask = np.ones(M.shape[0], dtype=bool)
        mask[rows_p] = False
        sub = M[:, tmark]
        tot = np.asarray(sub.sum(axis=1)).ravel()
        lib_sz = np.asarray(M.sum(axis=1)).ravel()
        lib_sz[lib_sz == 0] = 1
        frac = tot / lib_sz
        mp_, mu_ = float(frac[rows_p].mean()), float(frac[mask].mean())
        enrich = round(mp_ / mu_, 2) if mu_ > 0 else None

    clone = paired[cdr_c].astype(str).value_counts()
    cells = pd.DataFrame({
        "barcode": paired[bc_c].astype(str).values,
        "cdr3_aa": paired[cdr_c].astype(str).values,
        "v_gene": paired[v_c].astype(str).values if v_c else "",
        "j_gene": paired[j_c].astype(str).values if j_c else "",
        "clone_size": paired[cdr_c].astype(str).map(clone).values,
        "gse": gse, "library": lib,
    })
    X = M[paired["_row"].values, :]

    import scipy.sparse
    scipy.sparse.save_npz(f"{dest}/expr.npz", X.tocsr())
    cells.to_parquet(f"{dest}/cells.parquet", index=False)
    with open(f"{dest}/genes.txt", "w") as fh:
        fh.write("\n".join(genes))

    st = {
        "gse": gse, "library": lib, "status": "ok",
        "n_gex": len(bcs), "n_contigs": int(n_contigs), "n_trb": int(len(trb)),
        "n_paired": int(len(paired)), "pairing_rate": round(rate, 4),
        "n_genes": int(X.shape[1]), "n_antibody_dropped": n_ab,
        "n_clonotypes": int(clone.size),
        "max_clone": int(clone.max()), "singletons": int((clone == 1).sum()),
        "tcell_enrichment": enrich,
        # The enrichment test only discriminates when the rate is LOW. In a
        # T-sorted library the unpaired cells are also T cells, so the ratio is
        # ~1 even for a perfect join -- GSE272993 paired at 62.5% with
        # enrichment ~1 and the first version of this flag called it a bad join.
        # A high rate is its own evidence; only a low rate needs the marker test.
        "flag": ("BAD JOIN" if (rate < MIN_PAIRING and enrich is not None
                                and enrich < 1.5)
                 else "low rate, T-cell enriched, join OK" if rate < MIN_PAIRING
                 else ""),
        "secs": round(time.time() - t0, 1),
    }
    with open(f"{dest}/stats.json", "w") as fh:
        json.dump(st, fh)
    vol.commit()
    shutil.rmtree(work, ignore_errors=True)
    print(f"{gse}/{lib}: {len(paired):,} cells  pairing {100*rate:.1f}%  "
          f"{X.shape[1]:,} genes  {n_ab} antibody dropped"
          + ("  *** LOW PAIRING" if rate < MIN_PAIRING else ""))
    return st


@app.function(image=image, volumes={"/data": vol}, timeout=86400)
def build_all(rows: list[dict]):
    import os
    print(f"dispatching {len(rows)} libraries")
    out = list(build_one.map(rows, return_exceptions=True))
    ok = [r for r in out if isinstance(r, dict) and r.get("status") in ("ok", "cached")]
    print(f"\n{len(ok)}/{len(out)} libraries built")
    # The children wrote into OUT from their own containers; this one has to
    # reload the volume before it can see or write alongside their output.
    vol.reload()
    os.makedirs(OUT, exist_ok=True)
    with open(f"{OUT}/build_stats.json", "w") as fh:
        json.dump([r for r in out if isinstance(r, dict)], fh, indent=1)
    vol.commit()
    cells = sum(r.get("n_paired", 0) for r in ok)
    print(f"total cells: {cells:,}")
    for r in out:
        if not isinstance(r, dict):
            print(f"  EXCEPTION {r}")
        elif r.get("status") not in ("ok", "cached"):
            print(f"  FAILED {r.get('gse')}/{r.get('library')}: {r.get('status')}")
        elif r.get("flag"):
            print(f"  {r['flag']} {r['gse']}/{r['library']} "
                  f"{100*r['pairing_rate']:.1f}%")
    return out


@app.function(image=image, volumes={"/data": vol}, timeout=1800)
def report():
    import glob

    import pandas as pd

    rows = []
    for p in glob.glob(f"{OUT}/*/*/stats.json"):
        with open(p) as fh:
            rows.append(json.load(fh))
    if not rows:
        print("nothing built yet")
        return
    df = pd.DataFrame(rows)
    ok = df[df.get("status") == "ok"] if "status" in df else df
    print("=" * 92)
    print(f"{len(df)} libraries, {len(ok)} usable")
    print("=" * 92)
    if len(ok):
        g = (ok.groupby("gse")
               .agg(libraries=("library", "size"), cells=("n_paired", "sum"),
                    median_pairing=("pairing_rate", "median"),
                    clonotypes=("n_clonotypes", "sum"),
                    singletons=("singletons", "sum"),
                    antibody_dropped=("n_antibody_dropped", "sum"))
               .sort_values("cells", ascending=False))
        print(g.to_string())
        print(f"\nTOTAL CELLS: {int(ok['n_paired'].sum()):,}")
        print(f"TOTAL CLONOTYPES: {int(ok['n_clonotypes'].sum()):,}")
        print(f"SINGLETONS: {int(ok['singletons'].sum()):,} "
              f"({100*ok['singletons'].sum()/max(ok['n_clonotypes'].sum(),1):.1f}% "
              f"of clonotypes) -- the honest evaluation stratum")
        low = ok[ok["pairing_rate"] < MIN_PAIRING]
        if len(low):
            print(f"\n{len(low)} libraries below {100*MIN_PAIRING:.0f}% pairing "
                  f"-- inspect before use, these are usually mismatches:")
            print(low[["gse", "library", "n_gex", "n_trb",
                       "pairing_rate"]].to_string(index=False))
    bad = df[df.get("status") != "ok"] if "status" in df else df.iloc[:0]
    if len(bad):
        print(f"\n{len(bad)} not usable:")
        print(bad["status"].value_counts().to_string())


@app.local_entrypoint()
def main(dry: bool = False, limit: int = 0, gse: str = ""):
    rows = load_manifest()
    if gse:
        want = set(gse.split(","))
        rows = [r for r in rows if r["gse"] in want]
    if limit:
        # One library from each study, so a smoke test exercises every naming
        # convention rather than N variants of the easiest one.
        seen, picked = set(), []
        for r in rows:
            if r["gse"] not in seen:
                seen.add(r["gse"])
                picked.append(r)
            if len(picked) >= limit:
                break
        rows = picked
    gb = sum(r["bytes"] for r in rows) / 1e9
    print(f"manifest: {len(rows)} libraries, {gb:.2f} GB to move")
    by = {}
    for r in rows:
        by.setdefault(r["gse"], [0, 0])
        by[r["gse"]][0] += 1
        by[r["gse"]][1] += r["bytes"]
    for k, (n, b) in sorted(by.items(), key=lambda kv: -kv[1][1]):
        print(f"  {k:12s} {n:4d} libs  {b/1e9:7.2f} GB")
    if dry:
        print("\ndry run; nothing dispatched")
        return
    build_all.remote(rows)
