"""Can we recover the V gene these reads already contain? A two-library pilot.

Read-level CDR3 extraction finds the junction by pattern matching and never
compares the read to anything, so it reports the seam and discards the germline
V call. That is the wrong half to keep for an HLA-driven disease: CDR3 contacts
the peptide, but CDR1 and CDR2 -- which come entirely from the V segment --
contact the MHC molecule itself. HLA-B27 is MHC class I, and our own strongest
mechanistic result is that this clone sits 35x higher in CD8 than CD4 cells,
i.e. that it is class-I restricted. The V gene is the part of the receptor that
reads class I, and the published TRBV9 association in AS is a V-gene claim we
have never been able to test.

Measured first, on 200k reads of one library: 95% of junction-spanning reads
carry >= 20 aa (~60 nt) of V-segment sequence upstream of the conserved Y-x-C
anchor. That is ample for assignment, so this is a tooling gap rather than a
data gap.

Assignment here is deliberately simple: take the nucleotides immediately 5' of
the conserved cysteine, and compare them to the same region of every germline
TRBV allele. No aligner, no indel handling -- V segments are not spliced in
this region, so a fixed-offset comparison is well posed. MiXCR would do this
more thoroughly, but it needs a paid licence for commercial use; the germline
reference here comes from OLGA (GPL, IMGT-derived), so this pilot has no
licence encumbrance.

What the pilot answers: what fraction of reads get an unambiguous V call, how
concentrated the V repertoire is, and which V gene carries CASSVGLYSTDTQYF.

    python scripts/vgene_pilot.py --runs SRR5812652 SRR5812641
"""
from __future__ import annotations

import argparse
import collections
import gzip
import itertools
import json
import re
import sys
import urllib.parse
import urllib.request

CODON_TABLE = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
GENETIC_CODE = {a + b + c: CODON_TABLE[i] for i, (a, b, c)
                in enumerate(itertools.product("TCAG", repeat=3))}
AA = "[ACDEFGHIKLMNPQRSTVWY]"
JUNCTION_RE = re.compile(rf"[YF][FYLIVH]C({AA}{{4,30}}?[FW])G{AA}G")
RC = bytes.maketrans(b"ACGTN", b"TGCAN")

OLGA = ("https://raw.githubusercontent.com/statbiophys/OLGA/master/olga/"
        "default_models/human_T_beta/")
PROBE_NT = 25          # nucleotides 5' of the conserved C used for assignment

# 25, not 45, because read length varies across the deposit and sets the ceiling.
# At 150 bp the conserved C sits a median 78 nt into the read and 99.9% of
# junction-spanning reads clear a 45 nt probe; at 100 bp the median is 31 nt and
# only 0.2% do. Shortening costs almost nothing in specificity -- across 59 V
# genes, a 25 nt probe leaves 6 sequences shared by two genes against 4 at 45 nt
# -- and it is the difference between assigning 93% of a 150 bp library and 0.1%
# of a 100 bp one.
MAX_MISMATCH = 4       # beyond this the read is left unassigned
MIN_MARGIN = 2         # best match must beat the runner-up GENE by this much

# An absolute distance cutoff alone is not safe here, because germline V probes
# are not evenly spaced: TRBV9's 25 nt probe sits only 3 mismatches from TRBV13
# and TRBV7-4, inside a cutoff of 4. The margin rule is what makes a call
# trustworthy -- two independent sequencing errors would be needed to flip it.
# In practice this costs little: of 39 reads carrying CASSVGLYSTDTQYF, 35 match
# TRBV9 at distance 0 with a median margin of 3, and the other 4 sit 9 away and
# are correctly refused.
TARGET = "CASSVGLYSTDTQYF"


def germline_v(probe_nt: int = PROBE_NT) -> dict[str, str]:
    """{allele: the PROBE_NT nucleotides immediately 5' of the conserved C}.

    OLGA's anchor_index is the position of the conserved cysteine codon within
    each V allele's genomic sequence, which is exactly the landmark needed.
    """
    anchors = {}
    with urllib.request.urlopen(OLGA + "V_gene_CDR3_anchors.csv", timeout=120) as fh:
        for line in fh.read().decode().splitlines()[1:]:
            parts = line.split(",")
            if len(parts) >= 2 and parts[1].strip().isdigit():
                anchors[parts[0]] = int(parts[1])

    seqs = {}
    with urllib.request.urlopen(OLGA + "model_params.txt", timeout=120) as fh:
        section = None
        for line in fh.read().decode().splitlines():
            if line.startswith("#"):
                section = line
            elif line.startswith("%") and section and "v_choice" in section:
                name, _, seq = line[1:].partition(";")
                seqs[name] = seq.strip().upper()

    out = {}
    for allele, idx in anchors.items():
        s = seqs.get(allele)
        if s and idx >= probe_nt:
            out[allele] = s[idx - probe_nt:idx]
    return out


def fastq_url(run: str) -> list[str]:
    q = urllib.parse.urlencode({"accession": run, "result": "read_run",
                                "format": "json", "limit": "0",
                                "fields": "run_accession,fastq_ftp"})
    with urllib.request.urlopen(
            "https://www.ebi.ac.uk/ena/portal/api/filereport?" + q, timeout=180) as fh:
        rows = json.load(fh)
    return [u for u in (rows[0].get("fastq_ftp") or "").split(";") if u]


def assign(probe: str, ref: dict[str, str]) -> tuple[str | None, int]:
    """Nearest germline gene by Hamming distance, refusing ambiguous calls.

    Returns (gene, distance), or (None, distance) when the read is either too
    far from everything or too close to a second gene to be called.
    """
    best, bestd = None, 10**9
    runner = 10**9
    for allele, r in ref.items():
        d = sum(1 for a, b in zip(probe, r) if a != b)
        gene = allele.split("*")[0]
        if d < bestd:
            if best is not None and gene != best:
                runner = bestd
            best, bestd = gene, d
        elif gene != best and d < runner:
            runner = d
    if best is None or bestd > MAX_MISMATCH or runner - bestd < MIN_MARGIN:
        return None, bestd
    return best, bestd


def scan(run: str, ref: dict[str, str], n_reads: int, probe_nt: int,
         truncate: int = 0) -> dict:
    urls = fastq_url(run)
    if not urls:
        return {"run": run, "status": "no fastq url"}
    url = "https://" + (urls[1] if len(urls) > 1 else urls[0])   # mate 2 carries the CDR3

    v_counts = collections.Counter()
    pair_counts = collections.Counter()
    target_v = collections.Counter()
    n = junctions = assigned = 0
    with urllib.request.urlopen(url, timeout=600) as resp:
        gz = gzip.GzipFile(fileobj=resp)
        for i, line in enumerate(gz):
            if i % 4 != 1:
                continue
            n += 1
            raw = line.strip()
            if truncate:
                raw = raw[:truncate]
            nt = raw.translate(RC)[::-1].decode()
            for off in (0, 1, 2):
                s = nt[off:]
                s = s[:len(s) // 3 * 3]
                pep = "".join(GENETIC_CODE.get(s[j:j + 3], "X")
                              for j in range(0, len(s), 3))
                m = JUNCTION_RE.search(pep)
                if not m:
                    continue
                junctions += 1
                cdr3 = "C" + m.group(1)
                # nucleotide index of the conserved C, then the probe 5' of it
                c_nt = off + (m.start() + 2) * 3
                if c_nt < probe_nt:
                    break
                probe = nt[c_nt - probe_nt:c_nt]
                gene, d = assign(probe, ref)
                if gene:
                    assigned += 1
                    v_counts[gene] += 1
                    pair_counts[(gene, len(cdr3))] += 1
                    if cdr3 == TARGET:
                        target_v[gene] += 1
                break
            if n >= n_reads:
                break
    return {"run": run, "status": "ok", "n_reads": n, "junctions": junctions,
            "assigned": assigned, "v_counts": v_counts, "target_v": target_v}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--runs", nargs="+", required=True)
    ap.add_argument("--n-reads", type=int, default=400_000)
    ap.add_argument("--truncate", type=int, default=100,
                    help="trim every read to this many nt before assignment. "
                         "Required for any cross-group V-usage comparison: 11%% "
                         "of AS libraries here are long-read against 33%% of "
                         "healthy ones (p = 0.091), longer reads raise the "
                         "assignment rate gene-specifically, so an untruncated "
                         "comparison measures read length as well as biology.")
    ap.add_argument("--probe-nt", type=int, default=PROBE_NT,
                    help="nucleotides 5' of the conserved C used for assignment; "
                         "lower it for short-read libraries")
    args = ap.parse_args()

    ref = germline_v(args.probe_nt)
    genes = {a.split("*")[0] for a in ref}
    print(f"germline reference: {len(ref)} TRBV alleles, {len(genes)} genes, "
          f"{args.probe_nt} nt probe")
    if "TRBV9" not in genes:
        print("  WARNING: TRBV9 absent from the reference")

    for run in args.runs:
        r = scan(run, ref, args.n_reads, args.probe_nt, args.truncate)
        if r["status"] != "ok":
            print(f"\n{run}: {r['status']}")
            continue
        n, j, a = r["n_reads"], r["junctions"], r["assigned"]
        print(f"\n{run}: {n:,} reads -> {j:,} junctions ({j / n:.1%}) -> "
              f"{a:,} V-assigned ({a / max(j, 1):.1%} of junctions)")
        print("  top V genes:")
        for g, c in r["v_counts"].most_common(8):
            print(f"    {g:<12}{c:>8,}  {c / max(a, 1):6.2%}")
        tv = r["target_v"]
        if tv:
            print(f"  {TARGET} V gene: " +
                  ", ".join(f"{g} x{c}" for g, c in tv.most_common()))
        else:
            print(f"  {TARGET}: not seen in this subsample")


if __name__ == "__main__":
    main()
