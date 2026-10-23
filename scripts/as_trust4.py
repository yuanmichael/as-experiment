"""Assemble TCR repertoires for AS and PsA patients from Komech 2022.

PRJEB51315 (= E-MTAB-11498). Unlike the 2018 deposit, whose disease labels
survive only inside FASTQ filenames, every sample here is titled
`<disease>_<patient>_<tissue>_<subset>`, e.g. `as_Bel_SF_F`, so cases and
controls are unambiguous:

    10 AS patients, 17 PsA patients, 115 runs
    blood only: 7 AS + 16 PsA patients, 51 runs, ~16 GB paired

The handoff doc records this study as "29 AS + 86 PsA". Those are *run* counts.
Each patient contributes several runs (blood and synovial fluid x CD4, CD8 and
unsorted), so the real subject counts are 10 and 17.

Why Modal rather than the laptop: the work is embarrassingly parallel and the
dominant cost is pulling reads from ENA, which runs at datacentre bandwidth
here instead of over a home connection. 51 runs in parallel is ~10-15 minutes
against ~2.5 hours locally, and nothing lands on the laptop.

What it is looking for: CASSVGLYSTDTQYF (TRBV9/TRBJ2-3), the clonotype Komech
reported in AS synovium and the motif Faham reported in 30-39% of AS cases.
In GSE288581 it appeared in 3 of 4 AS patients against 15 of 786 healthy
subjects, and crucially against 6 of 44 *HLA-B27-positive* healthy subjects
(p = 0.017) -- the comparison that separates disease from HLA type. Four
patients cannot estimate a rate; ten can narrow it.

A caveat on the assembly: these are 5'RACE libraries that carry UMIs, and
TRUST4 does not use them, so clone *abundances* here are PCR-biased. That is
acceptable because the test is presence/absence per patient, not frequency.

    modal run scripts/as_trust4.py::list_runs
    modal run --detach scripts/as_trust4.py::run_all
    modal run scripts/as_trust4.py::collect
"""
import modal

app = modal.App("as-experiment")
vol = modal.Volume.from_name("as-data", create_if_missing=True)

PROJECT = "PRJEB51315"
OUT = f"/data/tcr/{PROJECT}"

# Blood and synovial fluid. SF was held back at first because it is a different
# compartment from the healthy blood controls, and it still is: SF carriage is
# never pooled into the blood denominator. It earns its place on its own terms.
#
#     AS SF   7 patients / 17 runs /  86M reads
#     PsA SF 17 patients / 47 runs / 245M reads
#
# Three things make the SF arm worth scanning rather than skipping:
#
#   * 3 AS patients (Lap, Non, Shor) appear only in SF, so they are invisible
#     to the blood scan. Lap and Shor are CD8-sorted, which is the compartment
#     the target clone occupies.
#   * 4 AS patients (Bel, Kos, Pon, Zar) have both PB and SF, giving a
#     within-person blood-versus-joint comparison rather than a group one.
#   * There is no healthy synovial fluid anywhere -- healthy joints have no
#     effusion to sample -- so the only possible control for AS synovium is
#     synovium from another arthritis. PsA SF, same lab and same 5'RACE assay,
#     is that control, and it makes the SF numbers interpretable at all.
#
# Blood runs already on the volume return their cached result, so widening this
# tuple re-scans nothing.
TISSUES = ("PB", "SF")

# AS only. The first dispatch tried all 51 blood runs and the real throughput
# turned out to be ~10 container-minutes per million reads, so 300M reads is
# ~50 container-hours (~$16) rather than the ~15 minutes first costed. The AS
# arm alone is 12 runs / 7 patients / 86M reads, and it is the whole headline
# test: 7 AS patients against the 6/44 healthy HLA-B27+ carriers already on
# disk. PsA answers a different question -- AS-specific versus
# spondyloarthritis-general -- and costs the same later as now.
# PsA added once the AS scan came back 4-5/7: it tests whether the clonotype
# marks AS specifically or spondyloarthritis generally. Already-scanned AS runs
# return their cached result, so re-including them costs seconds, not minutes.
DISEASES = ("as", "psa")

# TRUST4's thread count must match the container's CPU allocation. The first
# dispatch ran `-t 8` inside cpu=4.0, oversubscribing four cores. Eight cores
# with eight threads costs the same per core-hour but halves wall time, which
# also halves the memory-seconds bill.
CPUS = 8

MOTIF = "SVGLYSTDTQ"
FULL_CLONOTYPE = "CASSVGLYSTDTQYF"

image = (
    modal.Image.debian_slim(python_version="3.12")
    .apt_install("build-essential", "zlib1g-dev", "git", "wget", "curl", "unzip")
    .run_commands(
        "git clone https://github.com/liulab-dfci/TRUST4.git /opt/TRUST4",
        "cd /opt/TRUST4 && make -j4",
    )
    # Pinned: an unpinned pandas once broke a finished run at its last step.
    .pip_install("pandas==2.2.3", "pyarrow==17.0.0")
)


def parse_title(title: str) -> dict:
    """`as_Bel_SF_F` -> disease, patient, tissue, subset.

    The patient id is what matters: a clone seen in three libraries from one
    person is one person carrying it. Counting libraries instead would be the
    same pseudo-replication that halved this project's cell-level result.
    """
    parts = str(title).split("_")
    disease = parts[0] if parts else ""
    patient = parts[1] if len(parts) > 1 else ""
    tissue = next((p for p in parts if p in ("PB", "SF")), "")
    subset = ""
    if tissue in parts:
        i = parts.index(tissue)
        if len(parts) > i + 1 and parts[i + 1] in ("4", "8", "F"):
            subset = parts[i + 1]
    return dict(disease=disease, patient=patient, tissue=tissue, subset=subset)


def ena_runs(project: str) -> list[dict]:
    import json
    import urllib.parse
    import urllib.request

    url = "https://www.ebi.ac.uk/ena/portal/api/filereport?" + urllib.parse.urlencode({
        "accession": project, "result": "read_run", "format": "json", "limit": "0",
        "fields": "run_accession,sample_title,library_layout,read_count,fastq_ftp",
    })
    with urllib.request.urlopen(url, timeout=180) as fh:
        rows = json.load(fh)
    out = []
    for r in rows:
        meta = parse_title(r.get("sample_title", ""))
        if TISSUES and meta["tissue"] not in TISSUES:
            continue
        if DISEASES and meta["disease"] not in DISEASES:
            continue
        out.append({**r, **meta})
    return out


@app.function(image=image, timeout=600)
def list_runs():
    runs = ena_runs(PROJECT)
    for r in runs:
        print(f"{r['run_accession']}  {r['sample_title']:<22} {r['disease']:<4} "
              f"{r['patient']:<6} {int(r.get('read_count') or 0)/1e6:>6.1f}M")
    import collections
    pat = collections.defaultdict(set)
    for r in runs:
        pat[r["disease"]].add(r["patient"])
    print(f"\n{len(runs)} runs, " + ", ".join(
        f"{d}: {len(p)} patients" for d, p in sorted(pat.items())))
    return runs


@app.function(image=image, volumes={"/data": vol}, cpu=float(CPUS), memory=8192,
              timeout=21600, max_containers=8)
def trust4_one(run: dict) -> dict:
    """Download one run from ENA and assemble its repertoire.

    Resumable: a run whose report is already on the volume is skipped, so an
    interrupted dispatch costs nothing to repeat. Nothing here raises -- a
    failed run is reported as a row so the other 50 finish.
    """
    import os
    import subprocess
    import time

    acc = run["run_accession"]
    os.makedirs(OUT, exist_ok=True)
    report = f"{OUT}/{acc}_report.tsv"
    if os.path.exists(report):
        return {"run": acc, "status": "cached", **_meta(run)}

    work = f"/tmp/{acc}"
    os.makedirs(work, exist_ok=True)
    t0 = time.time()
    try:
        urls = [u for u in (run.get("fastq_ftp") or "").split(";") if u]
        if not urls:
            return {"run": acc, "status": "no fastq url", **_meta(run)}
        local = []
        for u in urls:
            dest = f"{work}/{os.path.basename(u)}"
            subprocess.run(["wget", "-q", "-O", dest, f"https://{u}"], check=True)
            local.append(dest)
        size = sum(os.path.getsize(f) for f in local) / 1e9
        print(f"{acc}: {size:.2f} GB in {(time.time() - t0) / 60:.1f} min", flush=True)

        base = ["/opt/TRUST4/run-trust4",
                "-f", "/opt/TRUST4/hg38_bcrtcr.fa",
                "--ref", "/opt/TRUST4/human_IMGT+C.fa",
                "-t", str(CPUS), "-o", acc, "--od", work]
        reads = ([["-1", local[0], "-2", local[1]], ["-u", local[0]]]
                 if len(local) >= 2 else [["-u", local[0]]])
        # ERR9030005 died on "mate-pair read files have different number of
        # reads" -- a truncated deposit at ENA. Mate 1 alone still yields CDR3s
        # from a 5'RACE library, so a paired failure falls back to single-end
        # rather than losing the patient.
        for attempt, extra in enumerate(reads):
            res = subprocess.run(base + extra, capture_output=True, text=True)
            if res.returncode == 0:
                break
            mode = "paired" if "-1" in extra else "single-end"
            print(f"{acc}: TRUST4 {mode} failed\n{res.stderr[-800:]}", flush=True)
            if attempt == len(reads) - 1:
                return {"run": acc, "status": "trust4 failed", **_meta(run)}

        src = f"{work}/{acc}_report.tsv"
        if not os.path.exists(src):
            return {"run": acc, "status": "no report", **_meta(run)}
        # Persist the expensive artefact before doing anything optional.
        with open(src) as fh:
            text = fh.read()
        with open(report, "w") as fh:
            fh.write(text)
        vol.commit()
    except Exception as exc:                            # noqa: BLE001
        return {"run": acc, "status": f"{type(exc).__name__}: {str(exc)[:120]}",
                **_meta(run)}

    n = max(len(text.splitlines()) - 1, 0)
    print(f"{acc}: {n:,} clonotypes, {(time.time() - t0) / 60:.1f} min", flush=True)
    return {"run": acc, "status": "ok", "n_clonotypes": n, **_meta(run)}


# --------------------------------------------------------------------------
# Read-level scan. Assembly is the wrong tool here: these are 5'RACE amplicon
# libraries, so nearly every read is TCR and TRUST4's assembly step -- built to
# run on the ~1% of reads RNA-seq leaves behind -- blew past a 6 h timeout.
#
# The question is presence/absence of one known sequence, which needs no
# assembly. A nucleotide regex built from the peptide's codon alternatives
# finds it in any frame at C speed, and is strictly more sensitive than
# assembly because no contig has to survive TRUST4's filters.

CODONS = {
    "A": "GC[ACGT]", "C": "TG[CT]", "D": "GA[CT]", "E": "GA[AG]",
    "F": "TT[CT]", "G": "GG[ACGT]", "H": "CA[CT]", "I": "AT[ACT]",
    "K": "AA[AG]", "L": "(?:CT[ACGT]|TT[AG])", "M": "ATG", "N": "AA[CT]",
    "P": "CC[ACGT]", "Q": "CA[AG]", "R": "(?:CG[ACGT]|AG[AG])",
    "S": "(?:TC[ACGT]|AG[CT])", "T": "AC[ACGT]", "V": "GT[ACGT]",
    "W": "TGG", "Y": "TA[CT]",
}

# A control pattern: nearly every TRB CDR3 starts CASS, so a high CASS count
# with a zero motif count means the scan worked and the sample is genuinely
# negative -- without it, zero is ambiguous between absence and a broken scan.
CONTROL_PEPTIDE = "CASS"


def peptide_regex(peptide: str) -> str:
    return "".join(CODONS[a] for a in peptide)


@app.function(image=image, volumes={"/data": vol}, cpu=2.0, memory=2048,
              timeout=3600, max_containers=8)
def scan_one(run: dict) -> dict:
    """Count reads containing the AS motif, at full depth, without assembling."""
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
        res["status"] = f"{type(exc).__name__}: {str(exc)[:120]}"
        return res

    res["minutes"] = round((time.time() - t0) / 60, 1)
    with open(out, "w") as fh:
        json.dump(res, fh)
    vol.commit()
    print(f"{acc} {run.get('sample_title')}: {res['n_reads']:,} reads, "
          f"motif={res['motif_reads']}, CASS={res['control_reads']:,}, "
          f"{res['minutes']} min", flush=True)
    return res


@app.function(image=image, volumes={"/data": vol}, timeout=86400)
def scan_all():
    """Dispatch the scan. `return_exceptions` is the fix the last two failures
    earned: a container timeout or OOM is not a Python exception inside the
    worker, so it propagated through `.map` and killed every sibling. Returned
    as a value instead, one bad run costs one run."""
    import collections

    runs = ena_runs(PROJECT)
    print(f"scanning {len(runs)} runs at full depth")
    results = list(scan_one.map(runs, return_exceptions=True))

    rows = [r for r in results if isinstance(r, dict)]
    for r in results:
        if not isinstance(r, dict):
            print(f"  worker died: {type(r).__name__}: {str(r)[:160]}")
    for r in rows:
        if r.get("status") != "ok":
            print(f"  FAILED {r['run']}: {r['status']}")

    per = collections.defaultdict(lambda: [0, 0, 0, []])
    for r in rows:
        if r.get("status") != "ok":
            continue
        p = per[(r["disease"], r["patient"])]
        p[0] += r["motif_reads"]
        p[1] += r["n_reads"]
        p[2] += r["control_reads"]
        p[3].append(f"{r['subset'] or '?'}:{r['motif_reads']}")
    print(f"\nper-patient carriage of {FULL_CLONOTYPE} ({MOTIF})")
    carriers = 0
    for (d, p), (mot, n, ctl, det) in sorted(per.items()):
        carriers += mot > 0
        print(f"  {d:<4} {p:<6} motif_reads={mot:<5} reads={n:>12,} "
              f"CASS={ctl:>10,}  [{' '.join(det)}]")
    print(f"\n{carriers}/{len(per)} patients carry the motif")
    return rows


# --------------------------------------------------------------------------
# CDR3 extraction. The motif scan answers "is this exact sequence present?"
# but T cells recognising one antigen usually share a *motif*, not a sequence,
# so a patient carrying CASSVGLYSTDTQF is biologically a carrier and the scan
# calls them negative.
#
# Still no assembly. The CDR3 sits between two conserved anchors -- the Cys at
# the end of the V gene and the FGxG at the start of the J gene -- both of
# which are findable in a translated read. Extracting what lies between gives
# the whole repertoire from the same single pass over the reads.

GENETIC_CODE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L", "CTT": "L", "CTC": "L",
    "CTA": "L", "CTG": "L", "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V", "TCT": "S", "TCC": "S",
    "TCA": "S", "TCG": "S", "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T", "GCT": "A", "GCC": "A",
    "GCA": "A", "GCG": "A", "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q", "AAT": "N", "AAC": "N",
    "AAA": "K", "AAG": "K", "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W", "CGT": "R", "CGC": "R",
    "CGA": "R", "CGG": "R", "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

# IMGT junction convention: starts at the conserved Cys, ends at the Phe of
# the J gene's FGxG. The anchor before the Cys ([YF][FYLIVH]) is FR3's
# conserved tail and keeps the match from firing on an incidental Cys.
JUNCTION_RE = r"[YF][FYLIVH]C([ACDEFGHIKLMNPQRSTVWY]{4,30}?[FW])G[ACDEFGHIKLMNPQRSTVWY]G"

# Keep every CDR3 within this edit distance of the target, however rare, and
# otherwise keep CDR3s seen twice or more. Singletons are dominated by PCR and
# sequencing error -- without UMIs there is no way to tell a real rare clone
# from a misread common one, so they inflate apparent diversity.
NEIGHBOUR_DIST = 2
MIN_READS = 2


def edit_distance(a: str, b: str, cap: int) -> int:
    """Levenshtein, abandoned once it exceeds `cap`."""
    if abs(len(a) - len(b)) > cap:
        return cap + 1
    prev = list(range(len(b) + 1))
    for i, ca in enumerate(a, 1):
        cur = [i] + [0] * len(b)
        for j, cb in enumerate(b, 1):
            cur[j] = min(prev[j] + 1, cur[j - 1] + 1, prev[j - 1] + (ca != cb))
        if min(cur) > cap:
            return cap + 1
        prev = cur
    return prev[-1]


@app.function(image=image, volumes={"/data": vol}, cpu=2.0, memory=8192,
              timeout=5400, max_containers=8)
def cdr3_one(run: dict) -> dict:
    """Extract every CDR3 in one run by translating reads, without assembling.

    Translation is vectorised in numpy: a batch of reads is joined, mapped to
    2-bit codes, reshaped into codons and looked up in one pass. Doing this a
    read at a time in Python would be slower than the download.
    """
    import gzip
    import json
    import os
    import re
    import subprocess
    import time
    from collections import Counter

    import numpy as np

    acc = run["run_accession"]
    os.makedirs(OUT, exist_ok=True)
    out = f"{OUT}/{acc}_cdr3.tsv.gz"
    meta_out = f"{OUT}/{acc}_cdr3.json"
    if os.path.exists(meta_out):
        with open(meta_out) as fh:
            return json.load(fh)

    bases = "ACGT"
    codon_aa = np.array(
        [GENETIC_CODE[bases[i] + bases[j] + bases[k]]
         for i in range(4) for j in range(4) for k in range(4)] + ["X"],
        dtype="U1")
    code = np.full(256, 4, dtype=np.uint8)
    for n, v in zip(b"ACGT", range(4)):
        code[n] = v
    comp = bytes.maketrans(b"ACGTN", b"TGCAN")
    rx = re.compile(JUNCTION_RE)

    def translated_frames(buf: bytes):
        arr = code[np.frombuffer(buf, dtype=np.uint8)]
        for frame in (0, 1, 2):
            a = arr[frame:]
            a = a[: len(a) - len(a) % 3].reshape(-1, 3).astype(np.int32)
            bad = (a == 4).any(axis=1)
            idx = a[:, 0] * 16 + a[:, 1] * 4 + a[:, 2]
            idx[bad] = 64
            yield "".join(codon_aa[idx])

    counts: Counter = Counter()
    n_reads = 0
    t0 = time.time()
    work = f"/tmp/{acc}"
    os.makedirs(work, exist_ok=True)
    res = {"run": acc, "status": "ok", **_meta(run)}
    try:
        urls = [u for u in (run.get("fastq_ftp") or "").split(";") if u]
        if not urls:
            res["status"] = "no fastq url"
            return res
        for u in urls:
            dest = f"{work}/{os.path.basename(u)}"
            subprocess.run(["wget", "-q", "-O", dest, f"https://{u}"], check=True)
            batch: list[bytes] = []
            with gzip.open(dest, "rb") as fh:
                for i, line in enumerate(fh):
                    if i % 4 != 1:
                        continue
                    n_reads += 1
                    batch.append(line.strip())
                    if len(batch) >= 200_000:
                        _harvest(batch, comp, translated_frames, rx, counts)
                        batch = []
            if batch:
                _harvest(batch, comp, translated_frames, rx, counts)
            os.remove(dest)
    except Exception as exc:                                # noqa: BLE001
        res["status"] = f"{type(exc).__name__}: {str(exc)[:160]}"
        return res

    keep = {c: n for c, n in counts.items() if n >= MIN_READS}
    near = {c: n for c, n in counts.items()
            if c not in keep
            and edit_distance(c, FULL_CLONOTYPE, NEIGHBOUR_DIST) <= NEIGHBOUR_DIST}
    keep.update(near)
    with gzip.open(out, "wt") as fh:
        fh.write("cdr3\treads\n")
        for c, n in sorted(keep.items(), key=lambda kv: -kv[1]):
            fh.write(f"{c}\t{n}\n")

    res.update(n_reads=n_reads, n_cdr3_all=len(counts), n_cdr3_kept=len(keep),
               n_singleton_neighbours=len(near),
               junction_reads=int(sum(counts.values())),
               minutes=round((time.time() - t0) / 60, 1))
    with open(meta_out, "w") as fh:
        json.dump(res, fh)
    vol.commit()
    print(f"{acc} {run.get('sample_title')}: {n_reads:,} reads, "
          f"{len(counts):,} CDR3s ({len(keep):,} kept), {res['minutes']} min",
          flush=True)
    return res


def _harvest(batch, comp, translated_frames, rx, counts) -> None:
    """Both orientations: 5'RACE mates are not consistently stranded."""
    joined = b"NNN".join(batch)
    for buf in (joined, joined.translate(comp)[::-1]):
        for pep in translated_frames(buf):
            for m in rx.finditer(pep):
                counts["C" + m.group(1)] += 1


@app.function(image=image, volumes={"/data": vol}, timeout=86400)
def cdr3_all():
    runs = ena_runs(PROJECT)
    print(f"extracting CDR3s from {len(runs)} runs")
    results = list(cdr3_one.map(runs, return_exceptions=True))
    rows = [r for r in results if isinstance(r, dict)]
    for r in results:
        if not isinstance(r, dict):
            print(f"  worker died: {type(r).__name__}: {str(r)[:160]}")
    for r in rows:
        if r.get("status") != "ok":
            print(f"  FAILED {r['run']}: {r['status']}")
    ok = [r for r in rows if r.get("status") == "ok"]
    print(f"\n{len(ok)}/{len(results)} runs extracted")
    for r in sorted(ok, key=lambda x: x.get("sample_title", "")):
        print(f"  {r['sample_title']:<22} {r['n_reads']:>12,} reads "
              f"{r['n_cdr3_kept']:>9,} CDR3s")
    return rows


def _meta(run: dict) -> dict:
    return {k: run.get(k, "") for k in ("sample_title", "disease", "patient",
                                        "tissue", "subset")}


# 24 h, not 2 h. The orchestrator blocks on `.map` until every worker returns,
# so its own timeout is a ceiling on the whole dispatch. At 7200 s it expired
# two hours in and Modal cancelled all eight workers mid-assembly, losing the
# run for 49 of 51 runs.
@app.function(image=image, volumes={"/data": vol}, timeout=86400)
def run_all():
    runs = ena_runs(PROJECT)
    print(f"dispatching {len(runs)} runs")
    results = list(trust4_one.map(runs))
    ok = [r for r in results if r.get("status") in ("ok", "cached")]
    print(f"\n{len(ok)}/{len(results)} succeeded")
    for r in results:
        if r.get("status") not in ("ok", "cached"):
            print(f"  FAILED {r['run']} ({r.get('sample_title')}): {r.get('status')}")
    return results


@app.function(image=image, volumes={"/data": vol}, timeout=7200)
def collect():
    """Per-patient carriage of the AS clonotype, plus the full clonotype table.

    Aggregated by patient, not by library: a patient counts as a carrier if any
    of their libraries contains the sequence.
    """
    import glob
    import os

    import pandas as pd

    meta = {r["run_accession"]: r for r in ena_runs(PROJECT)}
    rows, clones = [], []
    for path in sorted(glob.glob(f"{OUT}/*_report.tsv")):
        acc = os.path.basename(path).replace("_report.tsv", "")
        m = meta.get(acc)
        # The volume also holds two PsA reports from the first dispatch, which
        # DISEASES now excludes. Skipping them keeps blank-labelled rows out of
        # the per-patient table rather than counting them as a nameless patient.
        if m is None:
            continue
        try:
            df = pd.read_csv(path, sep="\t")
        except Exception as exc:                        # noqa: BLE001
            print(f"{acc}: unreadable ({exc})")
            continue
        df.columns = [c.strip("#") for c in df.columns]
        aa = df.get("CDR3aa", pd.Series(dtype=str)).astype(str)
        v = df.get("V", pd.Series(dtype=str)).astype(str)
        # TRUST4 writes "out_of_frame" / "partial" for unproductive contigs.
        good = aa.str.match(r"^[ACDEFGHIKLMNPQRSTVWY]+$", na=False) & (aa.str.len() >= 6)
        trb = v.str.startswith("TRBV") & good
        rows.append(dict(
            run=acc, sample=m.get("sample_title", ""), disease=m.get("disease", ""),
            patient=m.get("patient", ""), tissue=m.get("tissue", ""),
            subset=m.get("subset", ""),
            n_clonotypes=int(trb.sum()),
            trbv9=int((v.str.startswith("TRBV9") & trb).sum()),
            has_motif=int((aa.str.contains(MOTIF, na=False) & trb).sum()),
            has_full=int((aa == FULL_CLONOTYPE).sum()),
        ))
        sub = pd.DataFrame({"cdr3": aa[trb], "v": v[trb]})
        sub["run"] = acc
        sub["patient"] = m.get("patient", "")
        sub["disease"] = m.get("disease", "")
        clones.append(sub)

    summary = pd.DataFrame(rows)
    summary.to_csv(f"{OUT}/run_summary.tsv", sep="\t", index=False)
    if clones:
        allc = pd.concat(clones, ignore_index=True)
        allc.to_parquet(f"{OUT}/clonotypes.parquet", index=False)
    vol.commit()

    print(summary.to_string(index=False))
    print("\nper-patient carriage of", FULL_CLONOTYPE)
    if not summary.empty:
        per = (summary.groupby(["disease", "patient"])
                      .agg(libraries=("run", "size"),
                           clonotypes=("n_clonotypes", "sum"),
                           carries=("has_full", "max"))
                      .reset_index())
        print(per.to_string(index=False))
        for d, g in per.groupby("disease"):
            print(f"  {d}: {int(g['carries'].sum())}/{len(g)} patients "
                  f"= {100 * g['carries'].mean():.0f}%")
    return summary.to_dict("records")
