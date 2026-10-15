"""Per-study preprocessing for the AS-vs-control meta-analysis.

One study in, one (genes x samples) matrix plus a sample sheet out. Nothing is
merged across studies -- that only happens at the effect-size stage.

Microarray: parse the GEO series matrix, map probes to gene symbols with the
platform annotation, collapse duplicate probes by keeping the one with the
highest mean signal (averaging a dead probe into a live one destroys signal).

RNA-seq: read the supplementary table, log2(x+1), keep protein-coding genes.

Outputs data/bulk/prepared/<GSE>_expr.tsv.gz and <GSE>_samples.tsv.

Usage:  python scripts/bulk_prep.py
"""

from __future__ import annotations

import gzip
import re
import sys
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
BULK = ROOT / "data" / "bulk"
PREP = BULK / "prepared"
ANNOT = BULK / "platform_annot"

ANNOT_URL = "https://ftp.ncbi.nlm.nih.gov/geo/platforms/{stub}nnn/{gpl}/annot/{gpl}.annot.gz"


# --------------------------------------------------------------------------
# study definitions: how to read each one, and how to read its labels
# --------------------------------------------------------------------------

STUDIES: dict[str, dict] = {
    "GSE73754": dict(
        kind="matrix",
        gpl="GPL10558",
        tissue="whole_blood",
        logged=True,  # quantile-normalised in GeneSpring; check and fix if linear
        group_from=("characteristic", "disease"),
        case_pat=r"ankylosing",
        covars={"sex": ("characteristic", "sex"), "age": ("characteristic", "age")},
    ),
    "GSE25101": dict(
        kind="matrix",
        gpl="GPL6947",
        tissue="whole_blood",  # PAXgene whole blood per series design; the
        # 'cell type: PBMC' characteristic is a mis-annotation
        logged=True,  # lumi VST + robust spline normalisation
        group_from=("characteristic", "disease status"),
        case_pat=r"ankylosing",
        covars={},  # gender/age matched by design, not reported per sample
    ),
    "GSE18781": dict(
        kind="matrix",
        gpl="GPL570",
        tissue="peripheral_blood",
        logged=False,  # gcrma + invariant norm + median polish, linear scale
        group_from=("characteristic", "disease state"),
        case_pat=r"spondyloarthropathy",
        drop_pat=r"sarcoid",  # third arm; not part of this contrast
        covars={
            "sex": ("characteristic", "gender"),
            "age": ("characteristic", "age"),
            "set": ("characteristic", "set"),  # sets normalised independently
        },
    ),
    "GSE11886": dict(
        kind="matrix",
        gpl="GPL570",
        tissue="monocyte_derived_macrophage",
        logged=False,
        group_from=("title",),
        case_pat=r"\bPatient\b",
        keep_pat=r"\bNone\b",  # untreated arm only; IFN/LPS arms are perturbations
        covars={},
    ),
    "GSE299639": dict(
        kind="seq",
        file="GSE299639_genes_TPM.anno.txt.gz",
        tissue="peripheral_blood",
        value="tpm",  # already depth-normalised
        case_pat=r"^AS_",
        ctrl_pat=r"^HC_",
        sex_from_name=r"_([FM])\d+$",
    ),
    "GSE181364": dict(
        kind="seq",
        file="GSE181364_expressed_gene_reads.txt.gz",
        tissue="whole_blood",
        value="counts",  # raw reads -- MUST be depth-normalised to CPM
        case_pat=r"^AS_patient",
        ctrl_pat=r"^Health",
    ),
    "GSE221786": dict(
        kind="seq",
        file="GSE221786_matrix.txt.gz",
        tissue="il17_enriched_stimulated",
        value="counts",  # submitter calls it a matrix; depth varies ~10x, so
        # treat it as count-like and renormalise rather than trust the label
        case_pat=r"^AS_",
        ctrl_pat=r"^Ctrl_",
        sex_from_name=r"^[A-Za-z]+_([FM])_",
        id_cols=2,  # gene_id + gene
    ),
}


# --------------------------------------------------------------------------
# platform annotation
# --------------------------------------------------------------------------

def gpl_stub(gpl: str) -> str:
    num = gpl.removeprefix("GPL")
    return "GPL" + (num[:-3] if len(num) > 3 else "")


def probe_to_symbol(gpl: str) -> pd.Series:
    """probe id -> gene symbol, from the GEO platform .annot file."""
    ANNOT.mkdir(parents=True, exist_ok=True)
    local = ANNOT / f"{gpl}.annot.gz"
    if not local.exists():
        url = ANNOT_URL.format(stub=gpl_stub(gpl), gpl=gpl)
        req = urllib.request.Request(url, headers={"User-Agent": "as-experiment/1.0"})
        with urllib.request.urlopen(req, timeout=600) as fh:
            local.write_bytes(fh.read())
        print(f"    fetched {gpl}.annot.gz ({local.stat().st_size / 1e6:.1f} MB)")

    # The .annot file has a '!platform_table_begin' sentinel before the table.
    with gzip.open(local, "rt", encoding="utf-8", errors="replace") as fh:
        lines = fh.read().splitlines()
    start = next(i for i, l in enumerate(lines) if l.startswith("!platform_table_begin"))
    end = next(i for i, l in enumerate(lines) if l.startswith("!platform_table_end"))
    table = pd.read_csv(
        pd.io.common.StringIO("\n".join(lines[start + 1 : end])), sep="\t", dtype=str
    )
    idcol = table.columns[0]
    symcol = next(c for c in table.columns if "Gene symbol" in c)
    s = table.set_index(idcol)[symcol].dropna()
    # multi-mapping probes list several symbols; a probe that cannot be assigned
    # to exactly one gene is dropped rather than arbitrarily assigned
    s = s[~s.str.contains("///", na=False)]
    return s[s.str.strip() != ""]


# --------------------------------------------------------------------------
# series matrix parsing
# --------------------------------------------------------------------------

def read_series_matrix(path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Returns (expression probes x samples, metadata samples x fields)."""
    with open(path, encoding="utf-8", errors="replace") as fh:
        lines = fh.read().splitlines()

    header: dict[str, list[list[str]]] = {}
    tbl_start = None
    for i, line in enumerate(lines):
        if line.startswith("!series_matrix_table_begin"):
            tbl_start = i + 1
            break
        if line.startswith("!Sample_"):
            key, _, rest = line.partition("\t")
            vals = [v.strip().strip('"') for v in rest.split("\t")]
            header.setdefault(key.removeprefix("!Sample_"), []).append(vals)

    tbl_end = next(
        i for i, l in enumerate(lines) if l.startswith("!series_matrix_table_end")
    )
    expr = pd.read_csv(
        pd.io.common.StringIO("\n".join(lines[tbl_start:tbl_end])),
        sep="\t",
        index_col=0,
    )
    expr.columns = [c.strip().strip('"') for c in expr.columns]

    gsms = header["geo_accession"][0]
    meta = pd.DataFrame(index=gsms)
    meta["title"] = header["title"][0]
    # characteristics_ch1 appears once per field; each is "key: value"
    for row in header.get("characteristics_ch1", []):
        keys = {v.split(":", 1)[0].strip().lower() for v in row if ":" in v}
        if len(keys) != 1:
            # inconsistent field across samples -- store positionally
            continue
        key = keys.pop()
        meta[key] = [v.split(":", 1)[1].strip() if ":" in v else np.nan for v in row]
    return expr, meta


def label_matrix_study(gse: str, cfg: dict, meta: pd.DataFrame) -> pd.DataFrame:
    src = cfg["group_from"]
    field = "title" if src[0] == "title" else src[1]
    if field not in meta.columns:
        raise KeyError(f"{gse}: no metadata field {field!r}; have {list(meta.columns)}")
    txt = meta[field].astype(str)

    keep = pd.Series(True, index=meta.index)
    if cfg.get("drop_pat"):
        keep &= ~txt.str.contains(cfg["drop_pat"], case=False, regex=True)
    if cfg.get("keep_pat"):
        keep &= meta["title"].astype(str).str.contains(
            cfg["keep_pat"], case=False, regex=True
        )
    meta = meta[keep].copy()
    txt = txt[keep]

    meta["group"] = np.where(
        txt.str.contains(cfg["case_pat"], case=False, regex=True), "AS", "control"
    )
    for out, (_, field) in cfg.get("covars", {}).items():
        meta[out] = meta[field] if field in meta.columns else np.nan
    return meta


# --------------------------------------------------------------------------
# collapse probes -> genes
# --------------------------------------------------------------------------

def collapse_to_genes(expr: pd.DataFrame, sym: pd.Series) -> pd.DataFrame:
    shared = expr.index.intersection(sym.index)
    sub = expr.loc[shared]
    genes = sym.loc[shared].str.strip()
    # keep, per gene, the probe with the highest mean signal
    order = sub.mean(axis=1).groupby(genes).idxmax()
    out = sub.loc[order.values]
    out.index = order.index
    out.index.name = "gene"
    return out.sort_index()


# --------------------------------------------------------------------------
# per-kind drivers
# --------------------------------------------------------------------------

def prep_matrix(gse: str, cfg: dict) -> tuple[pd.DataFrame, pd.DataFrame]:
    expr, meta = read_series_matrix(BULK / gse / f"{gse}_series_matrix.txt")
    meta = label_matrix_study(gse, cfg, meta)
    expr = expr[meta.index]

    expr = expr.apply(pd.to_numeric, errors="coerce").dropna(how="all")
    # Trust the data, not the submitter's note: a 99th percentile in the
    # thousands means linear scale regardless of what the record claims.
    p99 = float(np.nanpercentile(expr.values, 99))
    if p99 > 100:
        print(f"    linear scale detected (p99={p99:.0f}) -> log2(x+1)")
        expr = np.log2(expr.clip(lower=0) + 1)
    else:
        print(f"    already log scale (p99={p99:.1f})")

    sym = probe_to_symbol(cfg["gpl"])
    expr = collapse_to_genes(expr, sym)

    # Probes near background measure noise, and noise with a small residual
    # variance produces a large standardised effect. Drop the bottom fifth.
    # Skipped for matrices the submitter already median-centred per gene, where
    # the mean carries no expression-level information.
    if abs(float(np.nanmedian(expr.values))) > 0.5:
        thresh = float(np.nanpercentile(expr.mean(axis=1), 20))
        before = len(expr)
        expr = expr[expr.mean(axis=1) > thresh]
        print(f"    background filter (mean > {thresh:.2f}): {before} -> {len(expr)}")
    else:
        print("    median-centred matrix: background filter not applicable")
    return expr, meta


def prep_seq(gse: str, cfg: dict) -> tuple[pd.DataFrame, pd.DataFrame]:
    path = BULK / gse / cfg["file"]
    df = pd.read_csv(path, sep="\t", dtype=str)

    tail = [c for c in ("Symbol", "Chr", "GeneType", "Description") if c in df.columns]
    n_id = cfg.get("id_cols", 1)
    id_cols = list(df.columns[:n_id])
    sample_cols = [c for c in df.columns if c not in id_cols + tail]

    if "Symbol" in df.columns:
        gene_col = "Symbol"
    elif "gene" in id_cols:
        gene_col = "gene"
    else:
        gene_col = id_cols[0]

    if "GeneType" in df.columns:
        pc = df["GeneType"].eq("protein_coding")
        print(f"    protein_coding {int(pc.sum())}/{len(df)}")
        df = df[pc].reset_index(drop=True)

    df = df[df[gene_col].notna() & (df[gene_col].astype(str).str.strip() != "")]
    mat = df[sample_cols].apply(pd.to_numeric, errors="coerce")
    mat = mat.set_axis(df[gene_col].astype(str).str.strip().values, axis=0)
    mat = mat.dropna(how="all")

    mat = mat.clip(lower=0)

    # Depth-normalise anything count-like. Skipping this is the fastest way to
    # manufacture a genome-wide false signal: if library size happens to track
    # disease status, EVERY gene inherits the difference and FDR cannot see it.
    if cfg["value"] == "counts":
        depth = mat.sum(axis=0)
        print(
            f"    library size: min {depth.min() / 1e6:.2f}M  max {depth.max() / 1e6:.2f}M "
            f"({depth.max() / max(depth.min(), 1):.1f}x spread) -> CPM"
        )
        mat = mat / depth * 1e6
    mat = np.log2(mat + 1)

    # Drop genes that are not meaningfully expressed. A gene that is near zero
    # in most samples has a residual variance close to zero, which sends the
    # standardised effect size to infinity at these sample sizes.
    expressed = (mat > 1).mean(axis=1)
    before = len(mat)
    mat = mat[expressed >= 0.5]
    print(f"    expression filter (>1 in >=50% of samples): {before} -> {len(mat)}")

    # One row per gene: keep the highest-expressed duplicate, as for probes.
    # Sort-then-drop rather than .loc[idxmax()] -- gene labels repeat, and .loc
    # with a duplicated label returns every matching row instead of selecting.
    mat = mat.iloc[np.argsort(-mat.mean(axis=1).to_numpy(), kind="stable")]
    mat = mat[~mat.index.duplicated(keep="first")]
    mat.index.name = "gene"
    mat = mat.sort_index()

    meta = pd.DataFrame(index=mat.columns)
    is_case = meta.index.str.contains(cfg["case_pat"], regex=True)
    is_ctrl = meta.index.str.contains(cfg["ctrl_pat"], regex=True)
    if not (is_case ^ is_ctrl).all():
        bad = meta.index[~(is_case ^ is_ctrl)].tolist()
        raise ValueError(f"{gse}: unclassifiable sample names {bad}")
    meta["group"] = np.where(is_case, "AS", "control")
    meta["title"] = meta.index
    if cfg.get("sex_from_name"):
        meta["sex"] = meta.index.str.extract(cfg["sex_from_name"], expand=False).values
    return mat, meta


def main() -> None:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    PREP.mkdir(parents=True, exist_ok=True)

    summary = []
    for gse, cfg in STUDIES.items():
        print(f"\n{gse}  ({cfg['kind']}, {cfg['tissue']})")
        try:
            expr, meta = (prep_matrix if cfg["kind"] == "matrix" else prep_seq)(gse, cfg)
        except Exception as exc:
            print(f"    FAILED: {type(exc).__name__}: {exc}")
            continue

        meta["gse"] = gse
        meta["tissue"] = cfg["tissue"]
        counts = meta["group"].value_counts().to_dict()
        print(f"    {expr.shape[0]} genes x {expr.shape[1]} samples   {counts}")

        expr.to_csv(PREP / f"{gse}_expr.tsv.gz", sep="\t", compression="gzip")
        meta.to_csv(PREP / f"{gse}_samples.tsv", sep="\t")
        summary.append(
            dict(
                gse=gse,
                tissue=cfg["tissue"],
                genes=expr.shape[0],
                AS=counts.get("AS", 0),
                control=counts.get("control", 0),
                covars=",".join(c for c in ("sex", "age", "set") if c in meta.columns),
            )
        )

    sdf = pd.DataFrame(summary)
    sdf.to_csv(PREP / "study_summary.tsv", sep="\t", index=False)
    print("\n" + sdf.to_string(index=False))


if __name__ == "__main__":
    main()
