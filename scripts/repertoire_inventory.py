"""Where does real TCR repertoire data live, including the places the ADC sweep could not see?

`airr_inventory.py` queried the AIRR Data Commons and found 5,356 repertoires,
zero of them AS. That number was wrong in a specific way: the ADC indexes
CURATED repertoire studies, so Komech's 115-run 5'RACE TCRb dataset was
invisible to it -- it sits in ArrayExpress/ENA, deposited as ordinary
sequencing.

This sweep looks in the three places that blind spot hides:

  1. GEO -- especially 10x 5' libraries with VDJ enrichment, which get indexed
     as "single-cell RNA-seq" and so never appear in a repertoire inventory
     despite being dedicated, deep TCR sequencing WITH a paired transcriptome.
  2. ENA -- the ground truth for raw reads, where ArrayExpress submissions land.
  3. BioStudies/ArrayExpress -- Komech's own home.

The criterion that matters is clonotype recovery per sample: >=10^4 separates
real repertoire data from TCR reads scraped out of bulk RNA-seq (our TRUST4 run
on GSE181364 recovered 87-332, i.e. 0.0009% of reads). Metadata cannot report
recovery, so assay type is the proxy -- AMPLICON / targeted / 5'RACE / VDJ
enrichment / immunoSEQ mean every read is a TCR read; RNA-Seq means they are
incidental.

Downloads nothing. Metadata only.

Usage:  python scripts/repertoire_inventory.py
"""

from __future__ import annotations

import json
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "airr"

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
ENA = "https://www.ebi.ac.uk/ena/portal/api"
BIOSTUDIES = "https://www.ebi.ac.uk/biostudies/api/v1"

UA = {"User-Agent": "as-experiment/1.0 (repertoire inventory; metadata only)"}

# Assays where every read is a TCR read. These are the ones worth having.
DEDICATED = [
    "immunoseq", "immunoseq", "tcr repertoire", "bcr repertoire",
    "t cell receptor repertoire", "tcrseq", "tcr-seq", "tcr seq",
    "5'race", "5' race", "race-", "amplicon", "targeted",
    "vdj", "v(d)j", "airr", "adaptive biotech", "multiplex pcr",
    "repertoire sequencing", "rep-seq", "repseq",
]

# 10x 5' + VDJ: dedicated TCR sequencing WITH a paired transcriptome. The
# whole point of the single-cell sequence->state project.
PAIRED = ["vdj", "v(d)j", "5'", "5 prime", "five prime", "tcr enrich",
          "immune profiling", "chromium"]

DISEASES = {
    "ankylosing spondylitis": ["ankylosing spondylitis", "axial spondyloarthritis",
                               "spondyloarthritis", "spondylarthritis"],
    "psoriatic arthritis": ["psoriatic arthritis"],
    "rheumatoid arthritis": ["rheumatoid arthritis"],
    "IBD": ["crohn", "ulcerative colitis", "inflammatory bowel"],
    "multiple sclerosis": ["multiple sclerosis"],
    "SLE / lupus": ["lupus", "systemic lupus"],
    "type 1 diabetes": ["type 1 diabetes", "type i diabetes"],
    "celiac": ["celiac", "coeliac"],
    "psoriasis": ["psoriasis"],
    "uveitis": ["uveitis"],
    "sjogren": ["sjogren", "sjögren"],
    "vasculitis": ["vasculitis", "behcet", "behçet"],
    "sarcoidosis": ["sarcoid"],
    "myositis": ["myositis", "dermatomyositis"],
}

HEALTHY = ["healthy", "control", "normal donor", "hc ", "volunteer"]


def get(url: str, timeout: int = 90) -> bytes | None:
    try:
        with urllib.request.urlopen(urllib.request.Request(url, headers=UA),
                                    timeout=timeout) as fh:
            return fh.read()
    except Exception as exc:
        print(f"    FAILED {type(exc).__name__}: {str(exc)[:110]}")
        return None


def any_of(text: str, keys: list[str]) -> bool:
    t = text.lower()
    return any(k in t for k in keys)


def which_disease(text: str) -> str | None:
    t = text.lower()
    for name, keys in DISEASES.items():
        if any(k in t for k in keys):
            return name
    return None


# ----------------------------------------------------------------- GEO

def geo_search(term: str, retmax: int = 400) -> list[str]:
    url = (f"{EUTILS}/esearch.fcgi?db=gds&retmax={retmax}&retmode=json"
           f"&term={urllib.parse.quote(term)}")
    raw = get(url)
    if not raw:
        return []
    try:
        return json.loads(raw)["esearchresult"].get("idlist", []) or []
    except Exception:
        return []


def geo_summaries(uids: list[str]) -> list[dict]:
    rows = []
    for i in range(0, len(uids), 150):
        chunk = uids[i:i + 150]
        url = (f"{EUTILS}/esummary.fcgi?db=gds&retmode=json"
               f"&id={','.join(chunk)}")
        raw = get(url)
        time.sleep(0.4)
        if not raw:
            continue
        try:
            res = json.loads(raw).get("result", {})
        except Exception:
            continue
        for uid in res.get("uids", []):
            r = res.get(uid, {})
            rows.append(dict(
                source="GEO",
                acc=r.get("accession", ""),
                title=(r.get("title") or "").strip(),
                summary=(r.get("summary") or "").strip(),
                n_samples=int(r.get("n_samples") or 0),
                taxon=r.get("taxon", ""),
                gdstype=r.get("gdstype", ""),
                pdat=r.get("PDAT", ""),
            ))
    return rows


def sweep_geo() -> pd.DataFrame:
    print("\n" + "=" * 92)
    print("GEO")
    print("=" * 92)

    queries = {
        "repertoire (any)": (
            '("T cell receptor repertoire"[All Fields] OR "TCR repertoire"[All Fields] '
            'OR "immunoSEQ"[All Fields] OR "TCR-seq"[All Fields] OR "Rep-seq"[All Fields] '
            'OR "immune repertoire"[All Fields]) AND "Homo sapiens"[Organism]'),
        "10x VDJ / paired": (
            '("VDJ"[All Fields] OR "V(D)J"[All Fields] OR "TCR enrichment"[All Fields] '
            'OR "immune profiling"[All Fields]) AND ("single cell"[All Fields]) '
            'AND "Homo sapiens"[Organism]'),
    }
    # One targeted query per inflammatory disease -- the whole point of the sweep
    for name, keys in DISEASES.items():
        dis = " OR ".join(f'"{k}"[All Fields]' for k in keys)
        queries[f"disease: {name}"] = (
            f'({dis}) AND ("TCR"[All Fields] OR "repertoire"[All Fields] '
            f'OR "VDJ"[All Fields] OR "T cell receptor"[All Fields]) '
            f'AND "Homo sapiens"[Organism]')

    seen, rows = set(), []
    for label, term in queries.items():
        uids = geo_search(term)
        time.sleep(0.4)
        new = [u for u in uids if u not in seen]
        seen.update(new)
        got = geo_summaries(new)
        for r in got:
            r["found_by"] = label
        rows += got
        print(f"  {label:34s} {len(uids):4d} hits, {len(got):4d} new")

    df = pd.DataFrame(rows)
    if df.empty:
        return df
    df = df[df["acc"].str.startswith("GSE")]
    df["text"] = df["title"] + " || " + df["summary"]
    df["disease"] = df["text"].apply(which_disease)
    df["dedicated"] = df["text"].apply(lambda t: any_of(t, DEDICATED))
    df["paired_sc"] = df["text"].apply(lambda t: any_of(t, PAIRED))
    df["has_control"] = df["text"].apply(lambda t: any_of(t, HEALTHY))
    return df


# ----------------------------------------------------------------- ENA

def sweep_ena() -> pd.DataFrame:
    print("\n" + "=" * 92)
    print("ENA  (raw reads -- where ArrayExpress submissions land, e.g. Komech)")
    print("=" * 92)

    terms = ["repertoire", "T cell receptor", "TCR", "immunoSEQ", "RACE TCR"]
    for name, keys in DISEASES.items():
        terms.append(keys[0])

    fields = ("study_accession,study_title,tax_id,library_strategy,"
              "library_selection,instrument_platform")
    rows, seen = [], set()
    for t in terms:
        q = f'study_title="*{t}*" AND tax_id=9606'
        url = (f"{ENA}/search?result=read_study&limit=0&format=json"
               f"&fields={urllib.parse.quote(fields)}"
               f"&query={urllib.parse.quote(q)}")
        raw = get(url)
        time.sleep(0.3)
        if not raw:
            print(f"  {t:28s} query failed")
            continue
        try:
            recs = json.loads(raw)
        except Exception:
            print(f"  {t:28s} unparseable")
            continue
        new = 0
        for r in recs:
            acc = r.get("study_accession", "")
            if not acc or acc in seen:
                continue
            seen.add(acc)
            new += 1
            rows.append(dict(source="ENA", acc=acc,
                             title=(r.get("study_title") or "").strip(),
                             summary="", n_samples=0, taxon="Homo sapiens",
                             gdstype=r.get("library_strategy", ""),
                             pdat="", found_by=f"ena: {t}"))
        print(f"  {t:28s} {len(recs):5d} studies, {new:4d} new")

    df = pd.DataFrame(rows)
    if df.empty:
        return df
    df["text"] = df["title"] + " || " + df["gdstype"]
    df["disease"] = df["text"].apply(which_disease)
    df["dedicated"] = df["text"].apply(lambda t: any_of(t, DEDICATED))
    df["paired_sc"] = df["text"].apply(lambda t: any_of(t, PAIRED))
    df["has_control"] = df["text"].apply(lambda t: any_of(t, HEALTHY))
    return df


def ena_run_stats(study: str) -> dict:
    """Run count and total bytes for one study -- the cost of using it."""
    fields = "run_accession,read_count,base_count,fastq_bytes,library_strategy"
    url = (f"{ENA}/filereport?accession={study}&result=read_run&format=json"
           f"&fields={urllib.parse.quote(fields)}")
    raw = get(url, timeout=120)
    if not raw:
        return {}
    try:
        recs = json.loads(raw)
    except Exception:
        return {}
    gb = 0.0
    for r in recs:
        for b in str(r.get("fastq_bytes") or "").split(";"):
            if b.strip().isdigit():
                gb += int(b) / 1e9
    reads = sum(int(r.get("read_count") or 0) for r in recs)
    strat = {r.get("library_strategy", "") for r in recs}
    return dict(n_runs=len(recs), total_gb=round(gb, 1),
                total_reads_M=round(reads / 1e6, 1),
                strategy=";".join(sorted(s for s in strat if s)))


# ----------------------------------------------------------- ArrayExpress

def sweep_arrayexpress() -> pd.DataFrame:
    print("\n" + "=" * 92)
    print("BioStudies / ArrayExpress  (Komech's own home)")
    print("=" * 92)

    terms = ["TCR repertoire", "T cell receptor repertoire", "immune repertoire"]
    for keys in DISEASES.values():
        terms.append(f"{keys[0]} TCR")

    rows, seen = [], set()
    for t in terms:
        url = (f"{BIOSTUDIES}/search?query={urllib.parse.quote(t)}"
               f"&collection=arrayexpress&pageSize=100")
        raw = get(url)
        time.sleep(0.3)
        if not raw:
            print(f"  {t:34s} query failed")
            continue
        try:
            hits = json.loads(raw).get("hits", []) or []
        except Exception:
            print(f"  {t:34s} unparseable")
            continue
        new = 0
        for h in hits:
            acc = h.get("accession", "")
            if not acc or acc in seen:
                continue
            seen.add(acc)
            new += 1
            rows.append(dict(source="ArrayExpress", acc=acc,
                             title=(h.get("title") or "").strip(),
                             summary="", n_samples=0, taxon="",
                             gdstype="", pdat=str(h.get("release_date") or ""),
                             found_by=f"ae: {t}"))
        print(f"  {t:34s} {len(hits):4d} hits, {new:3d} new")

    df = pd.DataFrame(rows)
    if df.empty:
        return df
    df["text"] = df["title"]
    df["disease"] = df["text"].apply(which_disease)
    df["dedicated"] = df["text"].apply(lambda t: any_of(t, DEDICATED))
    df["paired_sc"] = df["text"].apply(lambda t: any_of(t, PAIRED))
    df["has_control"] = df["text"].apply(lambda t: any_of(t, HEALTHY))
    return df


# ----------------------------------------------------------------- report

def main() -> None:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    OUT.mkdir(parents=True, exist_ok=True)

    parts = [sweep_geo(), sweep_ena(), sweep_arrayexpress()]
    parts = [p for p in parts if not p.empty]
    if not parts:
        print("\nnothing returned from any source")
        return
    df = pd.concat(parts, ignore_index=True)
    df = df.drop_duplicates(subset=["source", "acc"])
    df.to_csv(OUT / "repertoire_inventory_all.tsv", sep="\t", index=False)

    print("\n" + "=" * 92)
    print(f"TOTAL {len(df):,} candidate studies "
          f"({', '.join(f'{k}={v}' for k, v in df['source'].value_counts().items())})")
    print("=" * 92)

    # The question: inflammatory disease WITH dedicated repertoire sequencing.
    hit = df[df["disease"].notna()].copy()
    print(f"\nmentioning an inflammatory/autoimmune disease: {len(hit):,}")

    print("\n" + "-" * 92)
    print("BY DISEASE  (dedicated = assay where every read is a TCR read)")
    print("-" * 92)
    tab = (hit.groupby("disease")
              .agg(studies=("acc", "size"),
                   dedicated=("dedicated", "sum"),
                   paired_single_cell=("paired_sc", "sum"),
                   with_controls=("has_control", "sum"))
              .sort_values("dedicated", ascending=False))
    print(tab.to_string())

    print("\n" + "=" * 92)
    print("THE SHORTLIST -- inflammatory disease + dedicated repertoire assay")
    print("=" * 92)
    short = hit[hit["dedicated"] | hit["paired_sc"]].copy()
    short = short.sort_values(["disease", "source", "acc"])
    if short.empty:
        print("  none")
    else:
        for _, r in short.iterrows():
            flags = []
            if r["dedicated"]:
                flags.append("DEDICATED")
            if r["paired_sc"]:
                flags.append("PAIRED-SC")
            if r["has_control"]:
                flags.append("controls?")
            print(f"\n  {r['source']:13s} {r['acc']:16s} {r['disease']}  "
                  f"[{', '.join(flags)}]  n={r['n_samples'] or '?'}")
            print(f"    {r['title'][:150]}")

    short.to_csv(OUT / "repertoire_shortlist.tsv", sep="\t", index=False)

    # AS specifically, every hit, no assay filter -- small enough to eyeball
    print("\n" + "=" * 92)
    print("EVERY AS / SpA HIT (no assay filter)")
    print("=" * 92)
    as_hits = df[df["text"].str.lower().str.contains(
        "ankylosing|spondylo|spondylar", regex=True, na=False)]
    for _, r in as_hits.iterrows():
        print(f"  {r['source']:13s} {r['acc']:16s} n={r['n_samples'] or '?':>4} "
              f"ded={int(bool(r['dedicated']))} sc={int(bool(r['paired_sc']))}")
        print(f"    {r['title'][:150]}")

    # Cost of the ENA studies that survived -- runs and bytes, no download
    ena_short = short[short["source"] == "ENA"]["acc"].tolist()
    if ena_short:
        print("\n" + "=" * 92)
        print("SIZE OF THE ENA SHORTLIST (metadata only -- nothing downloaded)")
        print("=" * 92)
        for acc in ena_short[:25]:
            st = ena_run_stats(acc)
            time.sleep(0.3)
            if st:
                print(f"  {acc:14s} {st['n_runs']:4d} runs  {st['total_gb']:8.1f} GB  "
                      f"{st['total_reads_M']:9.1f}M reads  {st['strategy'][:40]}")

    print(f"\nwrote {OUT / 'repertoire_inventory_all.tsv'}")
    print(f"wrote {OUT / 'repertoire_shortlist.tsv'}")
    print("\nNOTE: 'dedicated' is inferred from title/summary text, which is a weak")
    print("proxy. Clonotype recovery per sample is the real criterion and can only")
    print("be measured from the data. Treat the shortlist as candidates to check.")


if __name__ == "__main__":
    main()
