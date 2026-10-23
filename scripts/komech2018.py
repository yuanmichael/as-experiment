"""Scan Komech 2018 (PRJNA393498) for the AS clonotype.

This deposit was written off as unlabelled because ENA's mirror serves an empty
`disease` field. The labels are there -- NCBI's BioSample XML carries
submitter-supplied `disease`, `health_state`, `Patient No`, `cell_type`,
`tissue`, `age` and `sex` for all 78 runs. So the cohort is:

    52 AS runs / 26 patients, 26 healthy runs
    unsorted PBMC only: 27 AS runs / 25 patients, 12 healthy runs

Only unsorted PBMC is used, and that matters more than it looks:

  * The 14 TRBV9-sorted runs are *sorted for the V gene our target uses*
    (TRBV9) and are all healthy donors. Including them would inflate healthy
    carriage and read as a refutation of the finding.
  * The 24 CD4/CD8-sorted AS runs sample a different slice of the repertoire
    than whole PBMC, so they are a different denominator, not extra patients.

Carriage is aggregated per *subject*, never per run: 27 runs come from 25
patients, and GE-style repeats and the TwHM triplicates would otherwise be
counted as separate people.

The scan logic is deliberately duplicated from as_trust4.py rather than
imported: that module builds its own modal.App at import time, and two apps in
one process is a worse problem than forty lines of repetition.

    python scripts/komech2018.py                  # list subjects, no compute
    modal run --detach scripts/komech2018.py::scan_all
    modal run scripts/komech2018.py::report
"""
import modal

app = modal.App("komech2018")
vol = modal.Volume.from_name("as-data", create_if_missing=True)

PROJECT = "PRJNA393498"
OUT = f"/data/tcr/{PROJECT}"

MOTIF = "SVGLYSTDTQ"
FULL_CLONOTYPE = "CASSVGLYSTDTQYF"
CONTROL_PEPTIDE = "CASS"

# Unsorted whole PBMC only -- see the module docstring.
KEEP_CELL_TYPES = ("PBMC",)

CODONS = {
    "A": "GC[ACGT]", "C": "TG[CT]", "D": "GA[CT]", "E": "GA[AG]",
    "F": "TT[CT]", "G": "GG[ACGT]", "H": "CA[CT]", "I": "AT[ACT]",
    "K": "AA[AG]", "L": "(?:CT[ACGT]|TT[AG])", "M": "ATG", "N": "AA[CT]",
    "P": "CC[ACGT]", "Q": "CA[AG]", "R": "(?:CG[ACGT]|AG[AG])",
    "S": "(?:TC[ACGT]|AG[CT])", "T": "AC[ACGT]", "V": "GT[ACGT]",
    "W": "TGG", "Y": "TA[CT]",
}

image = (
    modal.Image.debian_slim(python_version="3.12")
    .apt_install("wget")
    .pip_install("pandas==2.2.3", "pyarrow==17.0.0")
)


def peptide_regex(peptide: str) -> str:
    return "".join(CODONS[a] for a in peptide)


def ena_runs(project: str = PROJECT) -> list[dict]:
    """ENA for the files, NCBI for the labels, joined on BioSample accession."""
    import json
    import re
    import time
    import urllib.parse
    import urllib.request

    url = "https://www.ebi.ac.uk/ena/portal/api/filereport?" + urllib.parse.urlencode({
        "accession": project, "result": "read_run", "format": "json", "limit": "0",
        "fields": "run_accession,sample_accession,run_alias,read_count,fastq_ftp",
    })
    with urllib.request.urlopen(url, timeout=180) as fh:
        rows = json.load(fh)

    samples = sorted({r["sample_accession"] for r in rows})
    attrs: dict[str, dict] = {}
    for i in range(0, len(samples), 20):
        ids = ",".join(samples[i:i + 20])
        e = ("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?"
             f"db=biosample&id={ids}&rettype=full&retmode=xml")
        with urllib.request.urlopen(e, timeout=180) as fh:
            xml = fh.read().decode("utf-8", "replace")
        for blk in re.split(r"(?=<BioSample )", xml):
            acc = re.search(r'accession="(SAMN\d+)"', blk)
            if not acc:
                continue
            attrs[acc.group(1)] = {
                m.group(1): m.group(2) for m in re.finditer(
                    r'<Attribute attribute_name="([^"]+)"[^>]*>([^<]*)</Attribute>', blk)}
        time.sleep(0.4)

    out = []
    for r in rows:
        a = attrs.get(r["sample_accession"], {})
        cell = a.get("cell_type", "")
        if KEEP_CELL_TYPES and cell not in KEEP_CELL_TYPES:
            continue
        health = a.get("health_state", "") or a.get("disease", "")
        disease = "as" if health.lower().startswith("ankylos") else (
            "healthy" if health.lower().startswith("healthy") else "unknown")
        # AS runs carry `Patient No`. Healthy ones do not, so the alias has to
        # serve -- but its first token is a *sequencing plate*, not a donor
        # (AS patients P9, P15 and P16 all came off plate `P91`). Plate plus
        # index is unique per library; plate alone would merge two donors who
        # happened to share a plate.
        alias = r.get("run_alias") or r["run_accession"]
        subject = a.get("Patient No", "").replace("Patient ", "P") or \
            "_".join(re.split(r"[._-]", alias)[:2])
        out.append({**r, "cell_type": cell, "disease": disease,
                    "subject": subject, "age": a.get("age", ""),
                    "sex": a.get("sex", "")})
    return out


def _meta(run: dict) -> dict:
    return {k: run.get(k, "") for k in ("run_alias", "disease", "subject",
                                        "cell_type", "age", "sex")}


@app.function(image=image, volumes={"/data": vol}, cpu=2.0, memory=2048,
              timeout=3600, max_containers=8)
def scan_one(run: dict) -> dict:
    import gzip
    import json
    import os
    import re
    import subprocess
    import time

    acc = run["run_accession"]
    os.makedirs(OUT, exist_ok=True)
    out = f"{OUT}/{acc}_scan.json"
    if os.path.exists(out):
        with open(out) as fh:
            return json.load(fh)

    pat = re.compile(peptide_regex(MOTIF).encode())
    ctl = re.compile(peptide_regex(CONTROL_PEPTIDE).encode())
    rc = bytes.maketrans(b"ACGTN", b"TGCAN")

    work = f"/tmp/{acc}"
    os.makedirs(work, exist_ok=True)
    t0 = time.time()
    res = {"run": acc, "status": "ok", "n_reads": 0, "motif_reads": 0,
           "control_reads": 0, "hits": [], **_meta(run)}
    try:
        urls = [u for u in (run.get("fastq_ftp") or "").split(";") if u]
        if not urls:
            res["status"] = "no fastq url"
            return res
        for u in urls:
            dest = f"{work}/{os.path.basename(u)}"
            subprocess.run(["wget", "-q", "-O", dest, f"https://{u}"], check=True)
            with gzip.open(dest, "rb") as fh:
                for i, line in enumerate(fh):
                    if i % 4 != 1:
                        continue
                    res["n_reads"] += 1
                    seq = line.strip()
                    if ctl.search(seq):
                        res["control_reads"] += 1
                    m = pat.search(seq)
                    if m is None:
                        back = seq.translate(rc)[::-1]
                        m = pat.search(back)
                        seq = back
                    if m is not None:
                        res["motif_reads"] += 1
                        if len(res["hits"]) < 5:
                            lo = max(m.start() - 30, 0)
                            res["hits"].append(seq[lo:m.end() + 30].decode())
            os.remove(dest)
    except Exception as exc:                                # noqa: BLE001
        res["status"] = f"{type(exc).__name__}: {str(exc)[:160]}"
        return res

    res["minutes"] = round((time.time() - t0) / 60, 1)
    with open(out, "w") as fh:
        json.dump(res, fh)
    vol.commit()
    print(f"{acc} {run['disease']}/{run['subject']}: {res['n_reads']:,} reads, "
          f"motif={res['motif_reads']}, CASS={res['control_reads']:,}, "
          f"{res['minutes']} min", flush=True)
    return res


@app.function(image=image, volumes={"/data": vol}, timeout=86400)
def scan_all():
    runs = ena_runs()
    print(f"scanning {len(runs)} unsorted-PBMC runs")
    results = list(scan_one.map(runs, return_exceptions=True))
    rows = [r for r in results if isinstance(r, dict)]
    for r in results:
        if not isinstance(r, dict):
            print(f"  worker died: {type(r).__name__}: {str(r)[:160]}")
    for r in rows:
        if r.get("status") != "ok":
            print(f"  FAILED {r['run']}: {r['status']}")
    print(f"\n{sum(r.get('status') == 'ok' for r in rows)}/{len(results)} ok")
    return rows


def summarise(rows: list[dict], limit: float = 4.5e-6) -> None:
    """One row per subject. Carriage is called on frequency, not raw reads,
    because depth here spans 27-fold (0.72M to 19.7M reads)."""
    import collections

    per = collections.defaultdict(lambda: [0, 0, 0, []])
    for r in rows:
        if r.get("status") != "ok":
            continue
        p = per[(r["disease"], r["subject"])]
        p[0] += r["motif_reads"]
        p[1] += r["n_reads"]
        p[2] += r["control_reads"]
        p[3].append(r["run"])
    print(f"\n{'grp':<8}{'subject':<10}{'motif':>7}{'CASS':>12}{'freq':>11}  runs")
    cnt = collections.defaultdict(lambda: [0, 0])
    for (d, s), (mot, n, ctl, runs_) in sorted(per.items()):
        f = mot / ctl if ctl else 0.0
        cnt[d][0] += 1
        cnt[d][1] += f >= limit
        print(f"{d:<8}{s:<10}{mot:>7}{ctl:>12,}{f:>11.2e}  {len(runs_)}")
    print()
    for d, (n, k) in sorted(cnt.items()):
        print(f"  {d}: {k}/{n} carriers at freq >= {limit:.1e}")


@app.function(image=image, volumes={"/data": vol}, timeout=1800)
def report():
    import glob
    import json

    rows = []
    for p in sorted(glob.glob(f"{OUT}/*_scan.json")):
        with open(p) as fh:
            rows.append(json.load(fh))
    summarise(rows)
    return rows


if __name__ == "__main__":
    import collections

    runs = ena_runs()
    print(f"{len(runs)} unsorted-PBMC runs after filtering\n")
    print(f"{'run':<12}{'alias':<22}{'grp':<9}{'subject':<8}{'age':<5}{'sex':<8}{'reads':>12}")
    for r in sorted(runs, key=lambda x: (x["disease"], x["subject"])):
        print(f"{r['run_accession']:<12}{(r.get('run_alias') or '')[:21]:<22}"
              f"{r['disease']:<9}{r['subject']:<8}{r['age']:<5}{r['sex'][:7]:<8}"
              f"{int(r.get('read_count') or 0):>12,}")
    by = collections.defaultdict(set)
    for r in runs:
        by[r["disease"]].add(r["subject"])
    print()
    for d, s in sorted(by.items()):
        print(f"  {d}: {len(s)} subjects, "
              f"{sum(1 for r in runs if r['disease'] == d)} runs")
    miss = [r["run_accession"] for r in runs if not r.get("fastq_ftp")]
    print(f"  runs without a fastq url: {len(miss)} {miss[:5]}")
