"""Measure how much TCR a bulk RNA-seq deposit actually contains, on 2 runs.

The question this answers before any full run: can a transcriptome dataset
support a clonotype-frequency test at all? The limit is not the analysis, it is
the assay. Whole-blood bulk RNA-seq yields ~146-498 TCRb reads per sample
(0.0009% of reads), which puts the detection floor near 3e-3 -- thirty times
coarser than the 1e-4 threshold the AS result rests on.

GSE221786 might be the exception: PBMCs stimulated with CytoStim for 4 h and
enriched for IL-17-producing cells. Activated T cells transcribe far more TCR
mRNA than resting whole blood, so the TCR fraction could be orders of magnitude
higher. Worth measuring rather than assuming, and two runs settle it.

Read the ANCHOR count, not the CASS count. The first version of this pilot used
the degenerate `CASS` peptide as its denominator and reported 0.50% of reads as
TCR, which would have been a 500-fold lift over resting blood and a green light
to spend real money. It was noise. `CASS` expands to 12 nucleotides matching one
position in ~58,000, a 150bp read offers 139 positions, and both strands are
scanned: 278 x 1.7e-5 = 0.48% of reads hit *with no TCR present at all*. The two
pilot libraries came in at 0.5055% and 0.5036% -- on the chance line, and
agreeing far too precisely for biology. This never showed up in the 5'RACE
amplicon libraries because CASS ran at 20-80% there, a hundredfold above
background, so the noise was invisible.

The denominator here is instead an exact 25-mer from the TRB constant region,
inside a 173-nt stretch identical in TRBC1 and TRBC2 (Ensembl ENST00000633705 /
ENST00000466254), so one pattern counts both alleles. Chance expectation is
2.5e-13 per read: every hit is a real TCRb transcript. CASS is still counted,
now purely as the background check -- its ratio to chance is the diagnostic, and
near 1.0 means "no measurable TCR", not "lots of TCR".

Measured on a subsample, not the whole library. The first attempt scanned every
read of two 72M-read libraries and died at its 90-minute cap having learned
nothing. But the quantity wanted is a *fraction* -- what share of reads carry a
TCR -- and a fraction does not need the whole file. A million reads pins it to
the order of magnitude, which is all the verdict turns on, and the floor at full
depth follows by multiplication. Streaming also means the multi-GB download
never happens: the pipe is closed once the sample is in hand.

A constant-region count is not yet the number wanted: to call a clonotype a read
must span the whole CDR3, and the anchor only shows the transcript was there. The
two are related by how many read start positions satisfy each, so for 150bp reads
a CDR3-spanning read needs its start in a (150-45+1)=106 window against the
anchor's (150-25+1)=126, giving CDR3 reads ~= 0.84 x anchor reads. This assumes
coverage is even along the transcript; where it is not, RNA-seq leans 3', and the
constant region sits 3' of CDR3, so the assumption errs towards caution.

Further deliberately conservative choices, so the error runs towards NOT USABLE
and never towards a false green light:

  * Only mate 1 is read. Both mates of one fragment come from one cDNA
    molecule, so counting hits in both would double-count the molecule and
    flatter the floor.
  * The projection uses ENA's `read_count` (fragments), not reads.

    modal run scripts/rnaseq_pilot.py::pilot
"""
import modal

app = modal.App("rnaseq-pilot")
vol = modal.Volume.from_name("as-data", create_if_missing=True)

PROJECT = "PRJNA916206"          # GSE221786: 20 AS + 8 control, IL-17 enriched
N_RUNS = 2
N_SAMPLE = 4_000_000             # fragments sampled per library
MOTIF = "SVGLYSTDTQ"
CONTROL_PEPTIDE = "CASS"         # background check only -- see the docstring
TARGET_FLOOR = 1e-4              # the threshold the AS result is called at

# 173 nt identical in TRBC1 and TRBC2; the 25-mer is its middle.
TRBC_SHARED = ("GACCCGCAGCCCCTCAAGGAGCAGCCCGCCCTCAATGACTCCAGATACTGCCTGAGCAGC"
               "CGCCTGAGGGTCTCGGCCACCTTCTGGCAGAACCCCCGCAACCACTTCCGCTGTCAAGTC"
               "CAGTTCTACGGGCTCTCGGAGAATGACGAGTGGACCCAGGATAGGGCCAAACC")
ANCHOR = TRBC_SHARED[74:99]
CDR3_NT = 45                     # typical CDR3 nucleotide length

CODONS = {
    "A": "GC[ACGT]", "C": "TG[CT]", "D": "GA[CT]", "E": "GA[AG]",
    "F": "TT[CT]", "G": "GG[ACGT]", "H": "CA[CT]", "I": "AT[ACT]",
    "K": "AA[AG]", "L": "(?:CT[ACGT]|TT[AG])", "M": "ATG", "N": "AA[CT]",
    "P": "CC[ACGT]", "Q": "CA[AG]", "R": "(?:CG[ACGT]|AG[AG])",
    "S": "(?:TC[ACGT]|AG[CT])", "T": "AC[ACGT]", "V": "GT[ACGT]",
    "W": "TGG", "Y": "TA[CT]",
}

image = (modal.Image.debian_slim(python_version="3.12")
         .apt_install("wget")
         .pip_install("pandas==2.2.3"))


def peptide_regex(p: str) -> str:
    return "".join(CODONS[a] for a in p)


def ena_runs(project: str = PROJECT) -> list[dict]:
    import json
    import urllib.parse
    import urllib.request

    url = "https://www.ebi.ac.uk/ena/portal/api/filereport?" + urllib.parse.urlencode({
        "accession": project, "result": "read_run", "format": "json", "limit": "0",
        "fields": "run_accession,sample_title,read_count,fastq_ftp,library_strategy",
    })
    with urllib.request.urlopen(url, timeout=180) as fh:
        return json.load(fh)


@app.function(image=image, volumes={"/data": vol}, cpu=2.0, memory=2048,
              timeout=900, max_containers=4)
def scan(run: dict) -> dict:
    """Stream the first N_SAMPLE fragments of mate 1 and count TCR reads."""
    import gzip
    import json
    import os
    import re
    import subprocess
    import time

    acc = run["run_accession"]
    pat = re.compile(peptide_regex(MOTIF).encode())
    ctl = re.compile(peptide_regex(CONTROL_PEPTIDE).encode())
    anc = re.compile(ANCHOR.encode())
    rc = bytes.maketrans(b"ACGTN", b"TGCAN")
    total = int(run.get("read_count") or 0)
    t0 = time.time()
    res = {"run": acc, "title": run.get("sample_title", ""), "status": "ok",
           "sampled": 0, "cass": 0, "motif": 0, "anchor": 0, "readlen": 0,
           "read_count": total}
    proc = None
    try:
        urls = [u for u in (run.get("fastq_ftp") or "").split(";") if u]
        if not urls:
            res["status"] = "no fastq"
            return res
        # Pipe rather than download: the sample is in hand long before the
        # multi-GB file would be, and closing the pipe ends the transfer.
        proc = subprocess.Popen(["wget", "-q", "-O", "-", f"https://{urls[0]}"],
                                stdout=subprocess.PIPE)
        with gzip.GzipFile(fileobj=proc.stdout) as fh:
            for i, line in enumerate(fh):
                if i % 4 != 1:
                    continue
                res["sampled"] += 1
                s = line.strip()
                if not res["readlen"]:
                    res["readlen"] = len(s)
                back = s.translate(rc)[::-1]
                if anc.search(s) or anc.search(back):
                    res["anchor"] += 1
                if ctl.search(s) or ctl.search(back):
                    res["cass"] += 1
                if pat.search(s) or pat.search(back):
                    res["motif"] += 1
                if res["sampled"] >= N_SAMPLE:
                    break
    except Exception as exc:                                # noqa: BLE001
        res["status"] = f"{type(exc).__name__}: {str(exc)[:140]}"
        return res
    finally:
        if proc is not None:
            proc.kill()
            if proc.stdout:
                proc.stdout.close()
            proc.wait()

    res["minutes"] = round((time.time() - t0) / 60, 1)
    n = max(res["sampled"], 1)
    L = res["readlen"] or 150

    # The real measurement: TCRb transcripts per read, then CDR3-spanning reads
    # at full depth, which is what a clonotype frequency is actually called on.
    res["tcr_frac"] = res["anchor"] / n
    res["proj_anchor"] = res["tcr_frac"] * total
    res["proj_cdr3"] = res["proj_anchor"] * max(L - CDR3_NT + 1, 1) / max(L - len(ANCHOR) + 1, 1)
    res["floor"] = 1.0 / res["proj_cdr3"] if res["proj_cdr3"] else float("inf")

    # Background check. Near 1.0 means the CASS count is pure chance.
    # 1.717e-05 is P(CASS matches at one frame) = (2/64)(4/64)(6/64)(6/64).
    res["cass_frac"] = res["cass"] / n
    res["cass_expected"] = 1.0 - (1.0 - 1.717e-05) ** (2 * max(L - 12 + 1, 1))
    res["cass_ratio"] = res["cass_frac"] / res["cass_expected"]

    os.makedirs("/data/pilot", exist_ok=True)
    with open(f"/data/pilot/{acc}_pilot.json", "w") as fh:
        json.dump(res, fh)
    vol.commit()
    print(f"{acc} {res['title'][:24]}: {res['sampled']:,} sampled x {L}bp\n"
          f"    TRBC anchor={res['anchor']:,} ({100 * res['tcr_frac']:.5f}%) -> "
          f"{res['proj_anchor']:,.0f} TRB reads, {res['proj_cdr3']:,.0f} CDR3-spanning "
          f"at {total:,} fragments\n"
          f"    floor={res['floor']:.1e}  motif={res['motif']}  "
          f"CASS background={res['cass_ratio']:.2f}x chance  {res['minutes']} min",
          flush=True)
    return res


@app.function(image=image, volumes={"/data": vol}, timeout=1800)
def pilot():
    runs = ena_runs()
    print(f"{len(runs)} runs in {PROJECT}; "
          f"strategies: {sorted({r.get('library_strategy') for r in runs})}")
    print(f"sampling {N_SAMPLE:,} fragments per library, mate 1 only")
    pick = sorted(runs, key=lambda r: -int(r.get("read_count") or 0))[:N_RUNS]
    for r in pick:
        print(f"  piloting {r['run_accession']} {r.get('sample_title','')} "
              f"{int(r.get('read_count') or 0):,} reads")
    out = list(scan.map(pick, return_exceptions=True))
    rows = [r for r in out if isinstance(r, dict) and r.get("status") == "ok"]
    print("\nverdict")
    if not rows:
        print("  pilot failed; nothing to judge")
        return out
    for r in rows:
        print(f"  {r['run']}: CASS ran at {r['cass_ratio']:.2f}x its chance rate"
              + ("  <- background only, as expected" if r["cass_ratio"] < 2 else ""))
    thin = [r for r in rows if r["anchor"] < 20]
    best = max(r["floor"] for r in rows)          # worst case across the pilot
    print(f"\n  detection floor across pilot runs: {best:.1e}")
    if thin:
        print(f"  INCONCLUSIVE: {len(thin)}/{len(rows)} libraries gave under 20 "
              f"anchor reads in {N_SAMPLE:,}, so the estimate is too noisy to "
              f"project; lowest was {min(r['anchor'] for r in thin)}")
    elif best <= TARGET_FLOOR:
        print(f"  USABLE: floor is at or below {TARGET_FLOOR:.0e}, so the "
              f"frequency test can run on this deposit")
    else:
        print(f"  NOT USABLE: floor is {best / TARGET_FLOOR:.0f}x coarser than "
              f"{TARGET_FLOOR:.0e}; a zero here would mean nothing")
    return out
