"""What immune repertoire data actually exists publicly, and for which diseases?

Same exercise as the GEO sweep at the start of the project, pointed at the AIRR
Data Commons -- the federated set of repositories (iReceptor, VDJServer and
others) that index TCR/BCR repertoire studies behind one API.

The question that decides whether a repertoire-based model is feasible today:
how many repertoires exist for inflammatory disease, and how many of those
studies contain controls?

Usage:  python scripts/airr_inventory.py
"""

from __future__ import annotations

import json
import sys
import urllib.error
import urllib.request
from collections import Counter
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "airr"

# Public ADC (AIRR Data Commons) endpoints. The commons is federated, so each
# repository must be queried separately and the results pooled.
REPOS = {
    "ireceptor-ipa1": "https://ipa1.ireceptor.org/airr/v1",
    "ireceptor-ipa2": "https://ipa2.ireceptor.org/airr/v1",
    "ireceptor-ipa3": "https://ipa3.ireceptor.org/airr/v1",
    "ireceptor-ipa4": "https://ipa4.ireceptor.org/airr/v1",
    "ireceptor-ipa5": "https://ipa5.ireceptor.org/airr/v1",
    "vdjserver": "https://vdjserver.org/airr/v1",
    "covid19-1": "https://covid19-1.ireceptor.org/airr/v1",
    "covid19-2": "https://covid19-2.ireceptor.org/airr/v1",
    "covid19-3": "https://covid19-3.ireceptor.org/airr/v1",
    "covid19-4": "https://covid19-4.ireceptor.org/airr/v1",
}

INFLAMMATORY = [
    "spondyl", "ankylos", "arthrit", "crohn", "colitis", "inflammatory bowel",
    "lupus", "sclerosis", "psoria", "diabetes", "celiac", "coeliac", "sjogren",
    "vasculit", "uveitis", "myositis", "autoimmun", "ibd", "sle", "ra ",
    "thyroid", "graves", "hashimoto", "behcet", "sarcoid", "asthma", "allerg",
]


def post(url: str, payload: dict, timeout: int = 120) -> dict | None:
    req = urllib.request.Request(
        url, data=json.dumps(payload).encode(),
        headers={"Content-Type": "application/json", "Accept": "application/json",
                 "User-Agent": "as-experiment/1.0"},
    )
    try:
        with urllib.request.urlopen(req, timeout=timeout) as fh:
            return json.loads(fh.read().decode("utf-8", "replace"))
    except Exception as exc:
        print(f"    query failed: {type(exc).__name__}: {str(exc)[:120]}")
        return None


def facet(base: str, field: str) -> list[dict]:
    """Ask the repository to count repertoires grouped by one field."""
    res = post(f"{base}/repertoire", {"facets": field})
    if not res:
        return []
    return res.get("Facet", []) or []


def main() -> None:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    OUT.mkdir(parents=True, exist_ok=True)

    rows, study_rows = [], []
    for name, base in REPOS.items():
        print(f"\n{name}")
        info = None
        try:
            with urllib.request.urlopen(f"{base}/info", timeout=60) as fh:
                info = json.loads(fh.read().decode("utf-8", "replace"))
        except Exception as exc:
            print(f"  unreachable: {type(exc).__name__}")
            continue
        print(f"  online, API {info.get('version', '?')}")

        dis = facet(base, "subject.diagnosis.disease_diagnosis.label")
        if not dis:
            dis = facet(base, "subject.diagnosis.disease_diagnosis.id")
        total = sum(int(d.get("count") or 0) for d in dis)
        print(f"  {len(dis)} distinct diagnoses, {total:,} repertoires")
        for d in dis:
            label = str(d.get("subject.diagnosis.disease_diagnosis.label")
                        or d.get("subject.diagnosis.disease_diagnosis.id") or "").strip()
            rows.append(dict(repo=name, diagnosis=label or "(blank)",
                             n_repertoires=int(d.get("count") or 0)))

        st = facet(base, "study.study_id")
        for s in st:
            study_rows.append(dict(repo=name, study=str(s.get("study.study_id") or ""),
                                   n_repertoires=int(s.get("count") or 0)))
        print(f"  {len(st)} studies")

    if not rows:
        print("\nno data returned from any repository")
        return

    df = pd.DataFrame(rows)
    df.to_csv(OUT / "airr_diagnoses.tsv", sep="\t", index=False)
    pd.DataFrame(study_rows).to_csv(OUT / "airr_studies.tsv", sep="\t", index=False)

    agg = (df.groupby("diagnosis", as_index=False)["n_repertoires"].sum()
             .sort_values("n_repertoires", ascending=False))

    print("\n" + "=" * 88)
    print(f"TOTAL: {agg['n_repertoires'].sum():,} repertoires across "
          f"{df['repo'].nunique()} repositories, {len(agg)} diagnoses")
    print("=" * 88)
    print("\nTop 30 diagnoses by repertoire count:")
    print(agg.head(30).to_string(index=False))

    low = agg["diagnosis"].str.lower()
    mask = low.apply(lambda s: any(k in s for k in INFLAMMATORY))
    infl = agg[mask]
    print("\n" + "=" * 88)
    print("INFLAMMATORY / AUTOIMMUNE DIAGNOSES")
    print("=" * 88)
    print(infl.to_string(index=False) if len(infl) else "  none matched")
    print(f"\n  total inflammatory repertoires: {infl['n_repertoires'].sum():,}")

    healthy = agg[low.str.contains("healthy|normal|control|no disease|none",
                                   na=False)]
    print(f"  total healthy/control repertoires: {healthy['n_repertoires'].sum():,}")
    print(f"\nwrote {OUT / 'airr_diagnoses.tsv'}")


if __name__ == "__main__":
    main()
