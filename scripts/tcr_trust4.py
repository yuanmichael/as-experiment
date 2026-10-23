"""Pull TCR CDR3 sequences out of ordinary bulk RNA-seq, and look for the
AS-associated TRBV9 motif.

Standard Illumina RNA-seq fragments the whole transcript, so reads DO cross the
CDR3 junction -- just at ~1/1000 the depth of a dedicated TCR library. TRUST4
assembles those reads into clonotypes.

Target: the motif Faham reported in 30-39% of AS cases and 0/227 controls --
CDR3 amino acid "SVGLYSTDTQ", classically paired with TRBV9 (Komech's full
clonotype: CASSVGLYSTDTQYF).

Expect low sensitivity. Recovering ~10^2-10^3 clonotypes instead of ~10^6 means
only clones expanded above roughly 0.1-1% of the repertoire are visible, so a
negative result here does NOT mean the motif is absent.

    modal run --detach scripts/tcr_trust4.py::run_all
    modal run scripts/tcr_trust4.py::collect
"""
import modal

app = modal.App("as-experiment")
vol = modal.Volume.from_name("as-data", create_if_missing=True)
OUT = "/data/tcr"

image = (
    modal.Image.debian_slim(python_version="3.12")
    .apt_install("build-essential", "zlib1g-dev", "git", "wget", "curl", "unzip")
    .run_commands(
        "git clone https://github.com/liulab-dfci/TRUST4.git /opt/TRUST4",
        "cd /opt/TRUST4 && make -j4",
    )
    .pip_install("pandas", "pyarrow")
)

# Cohorts with raw reads in ENA. GSE299639 is Oxford Nanopore single-end; TRUST4
# is built for short paired reads and ONT error rates corrupt CDR3 calling, so
# it is excluded rather than run and quietly misinterpreted.
# Scoped deliberately: one small cohort, to validate the method end to end
# before spending anything larger. The expensive step is downloading reads, not
# analysing them, so the right design is to pull one cohort once and extract
# every sequence layer from it -- not one layer across many cohorts.
COHORTS = {
    "GSE181364": "PRJNA751846",   # 5 AS vs 3 HC, whole blood, 29 GB total
}
# Held back pending results: GSE221786 (PRJNA916206, 20 AS vs 8 HC, 222 GB)

MOTIF = "SVGLYSTDTQ"


def ena_runs(project: str) -> list[dict]:
    import json
    import urllib.parse
    import urllib.request

    url = "https://www.ebi.ac.uk/ena/portal/api/filereport?" + urllib.parse.urlencode({
        "accession": project, "result": "read_run", "format": "json",
        "fields": "run_accession,sample_title,library_layout,read_count,fastq_ftp",
    })
    with urllib.request.urlopen(url, timeout=180) as fh:
        return json.load(fh)


@app.function(image=image, volumes={"/data": vol}, timeout=600)
def list_runs():
    rows = []
    for gse, proj in COHORTS.items():
        for r in ena_runs(proj):
            rows.append({"gse": gse, **r})
    for r in rows:
        print(f"{r['gse']}  {r['run_accession']}  {r.get('sample_title','')[:30]:30s} "
              f"{r.get('library_layout')}  {int(r.get('read_count') or 0)/1e6:.1f}M")
    print(f"\n{len(rows)} runs total")
    return rows


@app.function(image=image, volumes={"/data": vol}, cpu=4.0, memory=16384,
              timeout=10800, max_containers=8)
def trust4_one(run: dict) -> dict:
    """Download one run's FASTQs from ENA and assemble its TCR repertoire."""
    import os
    import subprocess
    import time

    acc = run["run_accession"]
    gse = run["gse"]
    os.makedirs(f"{OUT}/{gse}", exist_ok=True)
    report = f"{OUT}/{gse}/{acc}_report.tsv"
    if os.path.exists(report):
        print(f"{acc}: already done")
        return {"run": acc, "gse": gse, "status": "cached"}

    work = f"/tmp/{acc}"
    os.makedirs(work, exist_ok=True)
    t0 = time.time()

    urls = [u for u in (run.get("fastq_ftp") or "").split(";") if u]
    if not urls:
        return {"run": acc, "gse": gse, "status": "no fastq url"}
    local = []
    for u in urls:
        dest = f"{work}/{os.path.basename(u)}"
        subprocess.run(["wget", "-q", "-O", dest, f"https://{u}"], check=True)
        local.append(dest)
    dl = time.time() - t0
    size = sum(os.path.getsize(f) for f in local) / 1e9
    print(f"{acc}: downloaded {size:.1f} GB in {dl/60:.1f} min")

    # TRUST4 wants R1/R2 for paired data; single files go through -u
    cmd = ["/opt/TRUST4/run-trust4",
           "-f", "/opt/TRUST4/hg38_bcrtcr.fa",
           "--ref", "/opt/TRUST4/human_IMGT+C.fa",
           "-t", "8", "-o", acc, "--od", work]
    if len(local) >= 2:
        cmd += ["-1", local[0], "-2", local[1]]
    else:
        cmd += ["-u", local[0]]

    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode != 0:
        print(f"{acc}: TRUST4 FAILED\n{res.stderr[-2000:]}")
        return {"run": acc, "gse": gse, "status": "trust4 failed"}

    src = f"{work}/{acc}_report.tsv"
    if not os.path.exists(src):
        return {"run": acc, "gse": gse, "status": "no report"}
    with open(src) as fh:
        text = fh.read()
    with open(report, "w") as fh:
        fh.write(text)
    vol.commit()

    n = max(len(text.splitlines()) - 1, 0)
    print(f"{acc}: {n} clonotypes, total {(time.time()-t0)/60:.1f} min")
    for f in local:
        os.remove(f)
    return {"run": acc, "gse": gse, "status": "ok", "n_clonotypes": n}


@app.function(image=image, volumes={"/data": vol}, timeout=3600)
def run_all():
    runs = []
    for gse, proj in COHORTS.items():
        for r in ena_runs(proj):
            runs.append({"gse": gse, **r})
    print(f"dispatching {len(runs)} runs")
    results = list(trust4_one.map(runs))
    ok = sum(1 for r in results if r.get("status") in ("ok", "cached"))
    print(f"\n{ok}/{len(results)} succeeded")
    for r in results:
        if r.get("status") not in ("ok", "cached"):
            print(f"  FAILED {r['run']}: {r.get('status')}")
    return results


@app.function(image=image, volumes={"/data": vol}, timeout=3600)
def collect():
    """Search every sample's repertoire for the AS motif."""
    import glob
    import os

    import pandas as pd

    rows = []
    meta = {}
    for gse, proj in COHORTS.items():
        for r in ena_runs(proj):
            meta[r["run_accession"]] = (gse, r.get("sample_title", ""))

    for path in sorted(glob.glob(f"{OUT}/*/*_report.tsv")):
        acc = os.path.basename(path).replace("_report.tsv", "")
        gse, title = meta.get(acc, ("?", "?"))
        try:
            df = pd.read_csv(path, sep="\t")
        except Exception as exc:
            print(f"{acc}: unreadable ({exc})")
            continue
        if df.empty:
            rows.append(dict(gse=gse, run=acc, sample=title, n_clono=0,
                             n_trb=0, trbv9=0, motif=0, motif_trbv9=0))
            continue

        df.columns = [c.strip("#") for c in df.columns]
        aa = df.get("CDR3aa", pd.Series(dtype=str)).astype(str)
        v = df.get("V", pd.Series(dtype=str)).astype(str)
        is_trb = v.str.startswith("TRBV")
        is_v9 = v.str.startswith("TRBV9")
        has_motif = aa.str.contains(MOTIF, na=False)

        rows.append(dict(
            gse=gse, run=acc, sample=title,
            n_clono=len(df), n_trb=int(is_trb.sum()),
            trbv9=int(is_v9.sum()),
            motif=int(has_motif.sum()),
            motif_trbv9=int((has_motif & is_v9).sum()),
        ))
        for _, row in df[has_motif].iterrows():
            print(f"  MOTIF HIT  {gse} {title} {acc}: {row.get('CDR3aa')} "
                  f"{row.get('V')} count={row.get('count')}")

    out = pd.DataFrame(rows)
    out["is_AS"] = out["sample"].str.contains(r"^AS|patient", case=False, regex=True)
    out.to_csv(f"{OUT}/motif_summary.tsv", sep="\t", index=False)
    vol.commit()
    print("\n" + out.to_string(index=False))

    print("\n" + "=" * 70)
    for gse, g in out.groupby("gse"):
        a, c = g[g["is_AS"]], g[~g["is_AS"]]
        print(f"{gse}: AS n={len(a)} median {a['n_clono'].median():.0f} clonotypes, "
              f"{int((a['motif'] > 0).sum())} with motif  |  "
              f"HC n={len(c)} median {c['n_clono'].median():.0f} clonotypes, "
              f"{int((c['motif'] > 0).sum())} with motif")
    return out.to_dict("records")
