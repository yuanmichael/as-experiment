"""Comprehensive sweep for spondyloarthritis datasets, across the whole spectrum.

Earlier sweeps queried "ankylosing spondylitis", "axial spondyloarthritis" and
"spondyloarthritis". That is the right spelling set for AS itself but misses the
rest of the SpA family, which shares the HLA-B27 association and in several
cases the same public TCR clonotypes: reactive arthritis (where this project's
clone was first described), enthesitis-related juvenile arthritis, acute
anterior uveitis, enteropathic arthritis, and the abbreviations axSpA/nr-axSpA
that recent deposits use in place of "ankylosing".

Two assay passes per disease term, because the two things we want are indexed
differently: dedicated TCR sequencing, and single-cell transcriptomes (which
often carry a V(D)J library that the summary text never mentions).

    python scripts/spa_inventory.py --out data/airr/spa_inventory.tsv
"""
from __future__ import annotations

import argparse
import json
import os
import time
import urllib.parse
import urllib.request

import pandas as pd

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
PAUSE = 0.4

DISEASES = {
    "ankylosing spondylitis": '"ankylosing spondylitis"[All Fields]',
    "axial spondyloarthritis": '"axial spondyloarthritis"[All Fields] OR axSpA[All Fields]',
    "spondyloarthritis (any)": ('"spondyloarthritis"[All Fields] OR '
                                '"spondyloarthropathy"[All Fields] OR '
                                '"spondylarthritis"[All Fields]'),
    "psoriatic arthritis": '"psoriatic arthritis"[All Fields]',
    "reactive arthritis": '"reactive arthritis"[All Fields] OR "Reiter"[All Fields]',
    "enthesitis-related arthritis": ('"enthesitis"[All Fields] OR '
                                     '"enthesitis-related arthritis"[All Fields]'),
    "acute anterior uveitis": ('"anterior uveitis"[All Fields] OR '
                               '"HLA-B27 uveitis"[All Fields]'),
    "enteropathic arthritis": ('"enteropathic arthritis"[All Fields] OR '
                               '"IBD-associated arthritis"[All Fields]'),
    "HLA-B27": '"HLA-B27"[All Fields] OR "HLA-B*27"[All Fields]',
    "sacroiliitis": '"sacroiliitis"[All Fields] OR "sacroiliac"[All Fields]',
}

ASSAYS = {
    "TCR": ('"T cell receptor"[All Fields] OR TCR[All Fields] OR '
            '"repertoire"[All Fields] OR "immunosequencing"[All Fields] OR '
            '"immunoSEQ"[All Fields] OR "VDJ"[All Fields] OR "TCR-seq"[All Fields]'),
    "single-cell": ('"single cell"[All Fields] OR "single-cell"[All Fields] OR '
                    'scRNA[All Fields] OR "10x Genomics"[All Fields] OR '
                    'CITE-seq[All Fields]'),
}


def get(url: str) -> dict:
    with urllib.request.urlopen(url, timeout=60) as fh:
        return json.loads(fh.read().decode())


def esearch(db: str, term: str, retmax: int = 500) -> list[str]:
    url = (f"{EUTILS}/esearch.fcgi?db={db}&retmax={retmax}&retmode=json"
           f"&term={urllib.parse.quote(term)}")
    try:
        return get(url)["esearchresult"].get("idlist", []) or []
    except Exception as exc:  # noqa: BLE001
        print(f"    search failed: {exc}")
        return []


def esummary(db: str, uids: list[str]) -> dict:
    out = {}
    for i in range(0, len(uids), 200):
        chunk = uids[i:i + 200]
        url = (f"{EUTILS}/esummary.fcgi?db={db}&retmode=json"
               f"&id={','.join(chunk)}")
        try:
            res = get(url).get("result", {})
        except Exception as exc:  # noqa: BLE001
            print(f"    summary failed: {exc}")
            continue
        for uid in res.get("uids", []):
            out[uid] = res[uid]
        time.sleep(PAUSE)
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default="data/airr/spa_inventory.tsv")
    ap.add_argument("--retmax", type=int, default=500)
    args = ap.parse_args()

    found: dict[str, dict] = {}
    for disease, dq in DISEASES.items():
        for assay, aq in ASSAYS.items():
            term = f"({dq}) AND ({aq}) AND (\"Homo sapiens\"[Organism])"
            uids = esearch("gds", term, args.retmax)
            time.sleep(PAUSE)
            print(f"  {disease:<30}{assay:<13}{len(uids):>5} hits")
            if not uids:
                continue
            for uid, rec in esummary("gds", uids).items():
                acc = rec.get("accession", uid)
                e = found.setdefault(acc, dict(
                    accession=acc, title=rec.get("title", ""),
                    summary=rec.get("summary", "")[:700],
                    n_samples=rec.get("n_samples", 0),
                    gdstype=rec.get("gdstype", ""), pdat=rec.get("pdat", ""),
                    taxon=rec.get("taxon", ""), diseases=set(), assays=set()))
                e["diseases"].add(disease)
                e["assays"].add(assay)

    df = pd.DataFrame([
        dict(r, diseases="; ".join(sorted(r["diseases"])),
             assays="; ".join(sorted(r["assays"]))) for r in found.values()])
    if df.empty:
        print("nothing returned")
        return
    df = df[df["taxon"].str.contains("Homo sapiens", na=False)]
    txt = (df["title"] + " " + df["summary"] + " " + df["gdstype"]).str.lower()
    df["likely_vdj"] = txt.str.contains(
        r"\bvdj\b|v\(d\)j|tcr|t cell receptor|repertoire|immunoseq", regex=True)
    df["likely_sc"] = txt.str.contains(
        r"single.cell|scrna|10x|cite-seq", regex=True)
    df = df.sort_values(["likely_vdj", "n_samples"], ascending=[False, False])

    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    df.to_csv(args.out, sep="\t", index=False)

    print(f"\n{len(df)} distinct human GEO series across the spectrum")
    print(df.groupby("diseases").size().sort_values(ascending=False).head(20).to_string())
    print(f"\nlikely TCR/VDJ: {int(df['likely_vdj'].sum())};  "
          f"likely single-cell: {int(df['likely_sc'].sum())}")
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
