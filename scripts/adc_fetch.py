"""Download already-processed repertoires from the AIRR Data Commons.

Why this exists: the synchronous ADC `/rearrangement` endpoint caps a page at
1000 rows. Emerson's 786 repertoires hold ~70M rearrangements, which is ~70,000
requests -- roughly 29 hours and a rude amount of load on a public academic
service. VDJServer also runs the **ADC async API**: submit one query, it builds
the whole result server-side and hands back a single file. One repertoire
(~90,000 rows, 7.8 MB) comes back in under 10 seconds.

Storage: the raw TSV for the full study is ~6 GB, almost all of it the
repertoire_id repeated on every row and allele-level V/J calls we do not use.
Each repertoire is reduced to four columns and written as parquet, which is
~1.5 MB. Nothing raw is kept.

Resumable: a repertoire already on disk is skipped, so an interrupted run
continues where it stopped.

Usage:
    python scripts/adc_fetch.py --study emerson-2017-natgen
    python scripts/adc_fetch.py --study emerson-2017-natgen --limit 20
"""

from __future__ import annotations

import argparse
import gzip
import io
import json
import sys
import time
import urllib.error
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]

# Repositories that serve the ADC async API. Only VDJServer is confirmed.
SYNC = "https://vdjserver.org/airr/v1"
ASYNC = "https://vdjserver.org/airr/async/v1"

FIELDS = ["repertoire_id", "v_call", "j_call", "junction_aa", "duplicate_count"]

QC_COLUMNS = ["repertoire_id", "subject_id", "n_raw", "n_kept", "templates",
              "n_unique_cdr3", "error"]

# Rearrangements per async query. Larger means fewer requests but a bigger
# server-side job and a longer retry if one fails.
REPS_PER_QUERY = 10
MAX_PARALLEL = 3
POLL_SECONDS = 6
POLL_TIMEOUT = 1800


def post(url: str, payload: dict, timeout: int = 180) -> dict:
    req = urllib.request.Request(
        url, data=json.dumps(payload).encode(),
        headers={"Content-Type": "application/json", "Accept": "application/json",
                 "User-Agent": "as-experiment/1.0"},
    )
    with urllib.request.urlopen(req, timeout=timeout) as fh:
        return json.loads(fh.read().decode("utf-8", "replace"))


def get(url: str, timeout: int = 120) -> dict:
    req = urllib.request.Request(url, headers={"User-Agent": "as-experiment/1.0"})
    with urllib.request.urlopen(req, timeout=timeout) as fh:
        return json.loads(fh.read().decode("utf-8", "replace"))


# --------------------------------------------------------------------------
# metadata


def fetch_metadata(study: str) -> pd.DataFrame:
    """One row per repertoire: ids, cohort, HLA alleles, tissue."""
    # The sync endpoint caps a page at 1000 for repertoires too, so page it.
    reps, start = [], 0
    while True:
        res = post(f"{SYNC}/repertoire", {
            "filters": {"op": "=", "content": {"field": "study.study_id", "value": study}},
            "size": 1000, "from": start,
        })
        page = res.get("Repertoire", [])
        reps.extend(page)
        if len(page) < 1000:
            break
        start += 1000

    rows = []
    for r in reps:
        subj = r.get("subject") or {}
        diag = (subj.get("diagnosis") or [{}])[0]
        samp = (r.get("sample") or [{}])[0]
        mhc = [a for a in (subj.get("mhc") or []) if a]
        rows.append(dict(
            repertoire_id=r.get("repertoire_id"),
            subject_id=subj.get("subject_id"),
            study_group=diag.get("study_group_description"),
            tissue=((samp.get("tissue") or {}).get("label")),
            template_class=samp.get("template_class"),
            # Documented alleles are the typed ones; "Inferred HLA-*" entries are
            # Emerson's own model output and must never be used as labels, since
            # that is the thing we are trying to predict.
            hla_documented=";".join(sorted(a for a in mhc if not a.startswith("Inferred"))),
            hla_inferred=";".join(sorted(a for a in mhc if a.startswith("Inferred"))),
        ))
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------
# rearrangements


def submit(rep_ids: list[str]) -> str:
    payload = {
        "filters": {"op": "in", "content": {"field": "repertoire_id", "value": rep_ids}},
        "fields": FIELDS,
        "format": "tsv",
    }
    return post(f"{ASYNC}/rearrangement", payload)["query_id"]


def wait(qid: str) -> dict:
    """Poll until the query leaves the running states."""
    deadline = time.time() + POLL_TIMEOUT
    last = None
    while time.time() < deadline:
        time.sleep(POLL_SECONDS)
        try:
            st = get(f"{ASYNC}/status/{qid}")
        except (urllib.error.URLError, TimeoutError) as exc:
            print(f"    status poll failed ({type(exc).__name__}), retrying")
            continue
        last = st
        status = (st.get("status") or "").upper()
        if status in ("FINISHED", "ERROR", "EXPIRED", "ABORTED"):
            return st
    return last or {"status": "TIMEOUT"}


def download_tsv(url: str) -> pd.DataFrame:
    req = urllib.request.Request(url, headers={"User-Agent": "as-experiment/1.0",
                                               "Accept-Encoding": "gzip"})
    with urllib.request.urlopen(req, timeout=900) as fh:
        raw = fh.read()
        if fh.headers.get("Content-Encoding") == "gzip":
            raw = gzip.decompress(raw)
    return pd.read_csv(io.BytesIO(raw), sep="\t", dtype=str)


def gene(call: str) -> str:
    """Allele-level call -> gene. 'TRBV3-1*01,TRBV3-2*01' -> 'TRBV3-1'.

    Ambiguous calls list several alleles; the first is taken, which is the
    convention IgBLAST output is normally collapsed with. Keeping allele
    resolution would split the same gene across columns for no gain -- the
    cell-level work already established that gene usage, not allele, is the
    usable signal.
    """
    if not isinstance(call, str) or not call:
        return ""
    return call.split(",")[0].split("*")[0]


def write_repertoire(df: pd.DataFrame, out: Path) -> dict:
    """Reduce one repertoire to four columns and persist. Returns a QC row."""
    n_raw = len(df)
    df = df.assign(
        v_gene=df["v_call"].map(gene),
        j_gene=df["j_call"].map(gene),
        junction_aa=df["junction_aa"].fillna(""),
        count=pd.to_numeric(df["duplicate_count"], errors="coerce").fillna(0).astype("int64"),
    )
    # A productive rearrangement has an in-frame CDR3 that starts with C, ends
    # with F or W, and carries no stop codon. The ADC copy is already filtered,
    # but this is cheap and the count is worth recording.
    prod = (df["junction_aa"].str.len() >= 6) & ~df["junction_aa"].str.contains(r"[*_]", regex=True)
    keep = df.loc[prod, ["v_gene", "j_gene", "junction_aa", "count"]]
    keep = keep[(keep["v_gene"] != "") & (keep["count"] > 0)]
    keep.to_parquet(out, index=False, compression="zstd")
    return dict(n_raw=n_raw, n_kept=len(keep), templates=int(keep["count"].sum()),
                n_unique_cdr3=keep["junction_aa"].nunique())


def fetch_batch(rep_ids: list[str], id2subject: dict[str, str], outdir: Path,
                attempts: int = 3) -> list[dict]:
    """Fetch one batch, retrying transient failures.

    Every failure mode here is transient and worth a retry: the async job can
    report FINISHED and still hand back a file that is empty or an error page.
    The first run of this lost 238 repertoires because one such file raised
    inside a thread pool and took the whole run down with it, so nothing below
    is allowed to raise -- a batch that cannot be fetched comes back as an
    error row and the other 77 batches carry on.
    """
    last = "not attempted"
    for attempt in range(1, attempts + 1):
        try:
            st = wait(submit(rep_ids))
            status = (st.get("status") or "").upper()
            if status != "FINISHED" or not st.get("download_url"):
                last = f"{status}: {st.get('message')}"
            else:
                df = download_tsv(st["download_url"])
                missing = [c for c in FIELDS if c not in df.columns]
                if missing:
                    # Usually a zero-row result, which has no header to parse.
                    last = f"result missing columns {missing}, {len(df)} rows"
                else:
                    return _write_batch(df, rep_ids, id2subject, outdir)
        except Exception as exc:                       # noqa: BLE001 - see docstring
            last = f"{type(exc).__name__}: {str(exc)[:160]}"
        print(f"    batch attempt {attempt}/{attempts} failed: {last}", flush=True)
        if attempt < attempts:
            time.sleep(10 * attempt)
    return [dict(repertoire_id=r, error=last) for r in rep_ids]


def _write_batch(df: pd.DataFrame, rep_ids: list[str],
                 id2subject: dict[str, str], outdir: Path) -> list[dict]:
    rows = []
    for rid, sub in df.groupby("repertoire_id"):
        name = id2subject.get(rid, rid)
        qc = write_repertoire(sub, outdir / f"{name}.parquet")
        rows.append(dict(repertoire_id=rid, subject_id=name, **qc))
    missing = set(rep_ids) - set(df["repertoire_id"].unique())
    for rid in missing:
        rows.append(dict(repertoire_id=rid, error="no rows returned"))
    return rows


def main() -> None:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    ap = argparse.ArgumentParser()
    ap.add_argument("--study", default="emerson-2017-natgen")
    ap.add_argument("--out", default=None, help="output dir (default data/airr/<study>)")
    ap.add_argument("--limit", type=int, default=0, help="only N repertoires, for a smoke test")
    ap.add_argument("--metadata-only", action="store_true")
    args = ap.parse_args()

    tag = args.study.replace("/", "_").replace(":", "_")
    outdir = Path(args.out) if args.out else ROOT / "data" / "airr" / tag
    reps_dir = outdir / "reps"
    reps_dir.mkdir(parents=True, exist_ok=True)

    print(f"study: {args.study}")
    meta = fetch_metadata(args.study)
    if meta.empty:
        print("no repertoires found")
        return
    meta.to_csv(outdir / "metadata.tsv", sep="\t", index=False)
    print(f"  {len(meta)} repertoires, {meta['subject_id'].nunique()} subjects")
    print(f"  with documented HLA: {(meta['hla_documented'] != '').sum()}")
    print(f"  wrote {outdir / 'metadata.tsv'}")
    if args.metadata_only:
        return

    id2subject = dict(zip(meta["repertoire_id"], meta["subject_id"]))
    todo = [r for r in meta["repertoire_id"]
            if not (reps_dir / f"{id2subject[r]}.parquet").exists()]
    print(f"  {len(meta) - len(todo)} already on disk, {len(todo)} missing")
    if args.limit:
        todo = todo[: args.limit]
        print(f"  --limit: fetching {len(todo)} of them")
    if not todo:
        return

    batches = [todo[i:i + REPS_PER_QUERY] for i in range(0, len(todo), REPS_PER_QUERY)]
    qc_path = outdir / "fetch_qc.tsv"
    done, t0 = 0, time.time()
    with ThreadPoolExecutor(max_workers=MAX_PARALLEL) as pool:
        for rows in pool.map(lambda b: fetch_batch(b, id2subject, reps_dir), batches):
            # Append after every batch so an interrupted run keeps its QC.
            # Reindexed to a fixed column set: appending frames whose columns
            # differ (an `error` column only when something failed) puts a
            # second header mid-file and nothing can read it back.
            pd.DataFrame(rows).reindex(columns=QC_COLUMNS).to_csv(
                qc_path, sep="\t", index=False, mode="a", header=not qc_path.exists())
            done += len(rows)
            el = time.time() - t0
            rate = done / el if el else 0
            eta = (len(todo) - done) / rate / 60 if rate else 0
            ok = sum(1 for r in rows if "error" not in r)
            print(f"  {done}/{len(todo)}  (+{ok} ok)  {el/60:.1f} min elapsed, "
                  f"~{eta:.0f} min left", flush=True)

    # The summary is a convenience, not the result. The repertoires are already
    # on disk and committed; a summary that cannot be computed must report that
    # and exit 0, never discard a finished run.
    n_on_disk = len(list(reps_dir.glob("*.parquet")))
    print(f"\ndone: {n_on_disk} repertoires on disk at {reps_dir}")
    try:
        qc = pd.read_csv(qc_path, sep="\t", on_bad_lines="warn")
        bad = qc[qc["error"].notna()] if "error" in qc else qc.iloc[:0]
        print(f"  QC rows: {len(qc) - len(bad)} ok, {len(bad)} failed")
        g = qc["n_kept"].dropna()
        print(f"  clones per repertoire: median {g.median():,.0f} "
              f"(min {g.min():,.0f}, max {g.max():,.0f})")
    except Exception as exc:                           # noqa: BLE001
        print(f"  QC summary unavailable ({type(exc).__name__}: {exc}); "
              f"repertoires are unaffected")


if __name__ == "__main__":
    main()
